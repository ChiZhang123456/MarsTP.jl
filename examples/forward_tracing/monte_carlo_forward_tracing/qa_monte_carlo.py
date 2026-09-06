"""Validate completed exports and summarize the targeted timestep refinement."""
from pathlib import Path
import json
import sys
import tomllib
import numpy as np
import h5py
from analyze_monte_carlo import load_csv

out=Path(sys.argv[1])
meta=tomllib.loads((out/"metadata.toml").read_text())
particles=load_csv(out/"particles.csv")
ref=load_csv(out/"timestep_refinement.csv")
lookup={int(p["particle_id"]):p for p in particles}
groups=0;checked=0;size=0
for path in sorted(out.glob("trajectories_*.jld2")):
    size+=path.stat().st_size
    with h5py.File(path) as f:
        names=sorted((k for k in f if k.startswith("p") and k[1:].isdigit()),key=lambda s:int(s[1:]))
        groups+=len(names)
        for name in set((names[0],names[-1])):
            state=f[f"{name}/state"][:]
            p=lookup[int(name[1:])]
            assert state.ndim==2 and state.shape[1]==7
            assert np.isfinite(state).all()
            assert np.all(np.diff(state[:,0])>0)
            assert np.all(np.diff(state[:,0])<=meta["dt_s"]+1e-10)
            radius=np.linalg.norm(state[:,1:4],axis=1)
            assert np.all(radius>=meta["Rm_m"]+meta["inner_altitude_km"]*1e3-1e-6)
            assert np.all(radius<=meta["outer_radius_Rm"]*meta["Rm_m"]+1e-6)
            weight_key="rate_weight_s1" if "rate_weight_s1" in f[name] else "weight_s1"
            assert f[f"{name}/{weight_key}"][()]==p["weight_s1"]
            np.testing.assert_allclose(state[0,1:],[p[k] for k in ("x0_m","y0_m","z0_m","vx0_ms","vy0_ms","vz0_ms")],rtol=0,atol=1e-9)
            np.testing.assert_allclose(state[-1,1:],[p[k] for k in ("xend_m","yend_m","zend_m","vxend_ms","vyend_ms","vzend_ms")],rtol=0,atol=1e-9)
            checked+=1
assert groups==meta["n_particles"]
summary={"saved_trajectory_groups":groups,"full_histories_checked":checked,"trajectory_bytes":size,
    "check_scope":"Count every particle group; validate first and last complete trajectory in every batch against CSV, domain bounds, times, finite values and weights",
    "refinement_scope":"Up to 200 largest coarse probe contributors plus positive-rate spaced source controls. Not a full-ensemble timestep convergence test."}
rates={}
for dt in (.1,.05,.025):
    a=ref[ref["dt_s"]==dt]
    rates[str(dt)]=float(np.sum(a["weight_s1"]*a["probe_residence_s"])/meta["cube_volume_m3"])
summary["refinement_subset_density_m3"]=rates
summary["density_change_01_to_005_relative"]=(rates["0.05"]-rates["0.1"])/rates["0.1"]
summary["density_change_005_to_0025_relative"]=(rates["0.025"]-rates["0.05"])/rates["0.05"]
summary["refinement_status_changes_01_to_005"]=int(np.sum(ref[ref["dt_s"]==.1]["status"]!=ref[ref["dt_s"]==.05]["status"]))
if (out/"energy_refinement.csv").exists():
    e=np.genfromtxt(out/"energy_refinement.csv",delimiter=",",names=True,dtype=None,encoding="utf8")
    summary["worst_energy_subset_max_abs_residual_eV"]={str(dt):float(np.max(abs(e["energy_residual_eV"][e["dt_s"]==dt]))) for dt in (.1,.05,.025)}
closure=abs(particles["residual_eV"])/np.maximum.reduce([abs(particles["deltaK_eV"]),abs(particles["work_eV"]),np.ones(len(particles))])
summary["particles_energy_closure_over_1_percent_with_1eV_floor"]=int(np.sum(closure>.01))
summary["probe_particles_max_energy_closure_relative"]=float(closure[particles["probe_residence_s"]>0].max())
summary["max_absolute_energy_residual_eV"]=float(np.max(abs(particles["residual_eV"])))
summary["limitations"]="Original dt=0.1 results retained. Probe sampling and flight-age convergence remain unresolved."
(out/"qa_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf8")
print(json.dumps(summary,indent=2))
