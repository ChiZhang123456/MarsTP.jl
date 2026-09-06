"""Re-bin finite-volume O2+ Monte Carlo residence records in SI units.

Usage: python analyze_monte_carlo.py RUN_DIR [--dv-kms 5] [--vmax-kms 300]
The two primary panels integrate over the omitted velocity. A separate figure
shows central-bin slab averages, which retain the 3D PSD units.
"""
from pathlib import Path
import argparse
import json
import math
import tomllib
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont


def load_csv(path):
    a=np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf8", ndmin=1)
    # Internal compatibility alias only. New files explicitly store rate_weight_s1.
    a.dtype.names=tuple("weight_s1" if n=="rate_weight_s1" else n for n in a.dtype.names)
    return a


def velocity_segments(v0, v1, edges):
    """Exact bin residence fractions for linearly interpolated velocity."""
    cuts = [0., 1.]
    for k in range(3):
        if v1[k] != v0[k]:
            f = (edges[k] - v0[k]) / (v1[k] - v0[k])
            cuts.extend(f[(f > 0) & (f < 1)])
    cuts = np.unique(cuts)
    for a, b in zip(cuts[:-1], cuts[1:]):
        vmid = v0 + .5 * (a+b) * (v1-v0)
        idx = tuple(int(np.searchsorted(edges[k], vmid[k], side="right")-1) for k in range(3))
        if any(i < 0 or i >= len(edges[k])-1 for k, i in enumerate(idx)):
            raise ValueError("Velocity outside requested bins; enlarge --vmax-kms")
        yield idx, b-a


def bin_residence(records, edges, volume):
    shape = tuple(len(e)-1 for e in edges)
    numbers = np.zeros(shape)
    # Each particle is one independent draw; multiple time samples are correlated.
    per_particle_bin = {}
    total_number = 0.
    current = np.zeros(3)
    for r in records:
        dt = r["t1_s"]-r["t0_s"]
        if dt <= 0 or r["weight_s1"] < 0:
            raise ValueError("Invalid dwell time/weight")
        v0 = np.array([r[f"v{k}0_ms"] for k in "xyz"])
        v1 = np.array([r[f"v{k}1_ms"] for k in "xyz"])
        number = r["weight_s1"] * dt
        total_number += number
        current += number*(v0+v1)/2/volume
        for idx, fraction in velocity_segments(v0, v1, edges):
            amount = number*fraction
            numbers[idx] += amount
            key = (int(r["particle_id"]), *idx)
            per_particle_bin[key] = per_particle_bin.get(key, 0.)+amount
    sumsq = np.zeros(shape)
    unique = np.zeros(shape, dtype=np.int32)
    for key, amount in per_particle_bin.items():
        idx = key[1:]
        sumsq[idx] += amount**2
        unique[idx] += 1
    neff = np.divide(numbers**2, sumsq, out=np.zeros(shape), where=sumsq > 0)
    dv3 = np.diff(edges[0])[:,None,None]*np.diff(edges[1])[None,:,None]*np.diff(edges[2])[None,None,:]
    f3d = numbers/volume/dv3
    assert np.isclose(numbers.sum(), total_number, rtol=1e-12, atol=1e-100)
    return f3d, numbers, unique, neff, current


def draw_panels(path, edges, a, b, label, subtitle, hits, neff, meta, zoom=False):
    mpl.rcParams.update({"font.family":"Arial", "font.size":10, "pdf.fonttype":42, "svg.fonttype":"none"})
    findfont("Arial", fallback_to_default=False)
    positive = np.r_[a[a>0], b[b>0]]
    norm = None
    if positive.size:
        vmin, vmax = positive.min(), positive.max()
        if vmax <= vmin:
            vmin, vmax = vmin/2, vmax*2
        norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, axes = plt.subplots(1,2,figsize=(10.4,4.8),layout="constrained")
    cmap=plt.get_cmap("turbo").copy()
    cmap.set_bad("#f0f0f0")
    supports=[]
    for matrix,ydim in zip((a,b),(1,2)):
        ix,iy=np.nonzero(matrix>0)
        supports.append(None if len(ix)==0 else (edges[0][ix.min()],edges[0][ix.max()+1],edges[ydim][iy.min()],edges[ydim][iy.max()+1]))
    available=[s for s in supports if s is not None]
    if zoom and available:
        xlo=min(s[0] for s in available);xhi=max(s[1] for s in available)
        span=max(xhi-xlo,max(s[3]-s[2] for s in available))+4*np.diff(edges[0])[0]
        xmid=(xlo+xhi)/2
    for ax, matrix, ydim, panel in zip(axes,(a,b),(1,2),("a","b")):
        im=ax.pcolormesh(edges[0]/1e3,edges[ydim]/1e3,np.ma.masked_less_equal(matrix.T,0),cmap=cmap,norm=norm,rasterized=True)
        ax.set(xlabel=r"$v_x$ (km/s)",ylabel=rf"$v_{'y' if ydim==1 else 'z'}$ (km/s)",aspect="equal")
        ax.set_title(f"({panel}) "+("Vx-Vy" if ydim==1 else "Vx-Vz"),loc="left")
        ax.axhline(0,lw=.5,c="gray",alpha=.4);ax.axvline(0,lw=.5,c="gray",alpha=.4)
        ax.set_facecolor("#f0f0f0")
        if not np.any(matrix>0):
            ax.text(.5,.5,"No sampled contribution in this view",ha="center",transform=ax.transAxes,fontsize=9)
        elif zoom:
            s=supports[ydim-1];ymid=(s[2]+s[3])/2
            ax.set_xlim((xmid-span/2)/1e3,(xmid+span/2)/1e3)
            ax.set_ylim((ymid-span/2)/1e3,(ymid+span/2)/1e3)
        if meta.get("plot_limit_kms") is not None:
            lim=meta["plot_limit_kms"]
            ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
            mx=(edges[0][:-1]<lim*1e3)&(edges[0][1:]>-lim*1e3)
            my=(edges[ydim][:-1]<lim*1e3)&(edges[ydim][1:]>-lim*1e3)
            if np.any(matrix>0) and not np.any(matrix[np.ix_(mx,my)]>0):
                ax.text(.5,.5,f"Sampled signal lies outside\nthe displayed ±{lim:g} km/s limits",ha="center",va="center",transform=ax.transAxes,fontsize=9)
    if positive.size:
        fig.colorbar(im,ax=axes,label=label,shrink=.78)
    probe=", ".join(f"{v:g}" for v in meta["detector_Rm"])
    fig.suptitle(f"O$_2^+$, probe ({probe}) $R_M$, cube side {meta['cube_side_Rm']:g} $R_M$\n"+subtitle)
    fig.savefig(path.with_suffix(".png"),dpi=190,bbox_inches="tight",pad_inches=.12)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    parser.add_argument("--dv-kms",type=float,default=5.)
    parser.add_argument("--vmax-kms",type=float,default=300.)
    parser.add_argument("--plot-limit-kms",type=float,default=300.,help="Axis display only; retain all velocity bins in saved PSD")
    parser.add_argument("--output-dir",type=Path,help="New directory for derived outputs; read original run without overwriting its analysis")
    args=parser.parse_args()
    if args.dv_kms <= 0:
        raise ValueError("dv must be positive")
    out=args.run
    meta=tomllib.loads((out/"metadata.toml").read_text())
    meta["plot_limit_kms"]=args.plot_limit_kms
    complete=tomllib.loads((out/"completion.toml").read_text())
    assert complete["complete"]
    particles=load_csv(out/"particles.csv")
    cells=load_csv(out/"source_cells.csv")
    records=load_csv(out/"probe_residence.csv")
    events=load_csv(out/"probe_crossings.csv")
    if args.output_dir is not None:
        out=args.output_dir
        out.mkdir(parents=True,exist_ok=False)
    assert len(particles)==meta["n_particles"]
    starts=np.searchsorted(particles["cell_id"],cells["cell_id"],side="left")
    ends=np.searchsorted(particles["cell_id"],cells["cell_id"],side="right")
    analytic_rate=0.;rate_variance=0.
    for r,start,end in zip(cells,starts,ends):
        p=particles[start:end]
        w=p["weight_s1"]
        assert len(w)==meta["per_cell"]
        if meta["flux_model"]=="reservoir_maxwellian_rate":
            x=np.column_stack([p[f"{k}0_m"] for k in "xyz"])
            v=np.column_stack([p[f"v{k}0_ms"] for k in "xyz"])
            vr=np.sum(v*x,axis=1)/np.linalg.norm(x,axis=1)
            expected=r["n_m3"]*r["area_m2"]*np.maximum(vr,0)*np.exp(p["log_importance"])/meta["per_cell"]
            # Near-tangent v dot er suffers cancellation. Bound the difference
            # from recomputing the normalized direction after CSV roundtrip.
            rounding=64*np.finfo(float).eps*r["n_m3"]*r["area_m2"]*np.exp(p["log_importance"])/meta["per_cell"]*np.linalg.norm(v,axis=1)
            assert np.all(abs(w-expected)<=1e-12*abs(expected)+rounding)
            assert np.isclose(p["source_density_weight_m3"].sum(),r["n_m3"],rtol=1e-12)
            normal=x[0]/np.linalg.norm(x[0])
            ur=np.dot(normal,[r["ux_ms"],r["uy_ms"],r["uz_ms"]])
            sigma=math.sqrt(1.380649e-23*r["Ti_K"]/meta["particle_mass_kg"])
            a=ur/sigma
            analytic_rate+=r["n_m3"]*r["area_m2"]*(sigma*math.exp(-a*a/2)/math.sqrt(2*math.pi)+ur*.5*math.erfc(-a/math.sqrt(2)))
            N=meta["per_cell"]
            if N>1:
                rate_variance+=N/(N-1)*np.sum((w-w.mean())**2)
        else:
            assert np.isclose(w.sum(),r["flux_m2_s"]*r["area_m2"],rtol=1e-12,atol=1e-100)
    vmax=100e3
    if len(records):
        vmax=max(vmax,max(float(np.max(abs(records[f"v{k}{i}_ms"]))) for k in "xyz" for i in (0,1)))
    dv=args.dv_kms*1e3
    extent=(np.ceil(vmax/dv)+.5)*dv if args.vmax_kms is None else args.vmax_kms*1e3
    if extent <= 0:
        raise ValueError("vmax must be positive")
    if not np.isclose(2*extent/dv,round(2*extent/dv)):
        raise ValueError("The full velocity span 2*vmax must be an integer multiple of dv")
    # Centers at zero for central-slab views; custom extent must align to dv.
    edge=np.arange(-extent,extent+dv*.1,dv)
    edges=(edge,edge.copy(),edge.copy())
    f3d,number,unique,neff,current=bin_residence(records,edges,meta["cube_volume_m3"])
    fxy=np.sum(f3d*np.diff(edge)[None,None,:],axis=2)
    fxz=np.sum(f3d*np.diff(edge)[None,:,None],axis=1)
    center=int(np.searchsorted(edge,0,side="right")-1)
    slice_xy=f3d[:,:,center];slice_xz=f3d[:,center,:]
    density=number.sum()/meta["cube_volume_m3"]
    dwell_number=particles["weight_s1"]*particles["probe_residence_s"]
    assert np.isclose(density,dwell_number.sum()/meta["cube_volume_m3"],rtol=1e-10,atol=1e-100)
    hit_ids=particles["particle_id"][particles["probe_residence_s"]>0]
    effective=dwell_number.sum()**2/np.sum(dwell_number**2) if np.any(dwell_number) else 0.
    cell_neff=[]
    for start,end in zip(starts,ends):
        w=particles["weight_s1"][start:end]
        if np.any(w):cell_neff.append(w.sum()**2/np.sum(w*w))
    faces={}
    for face in (-1,1,-2,2,-3,3):
        e=events[events["face"]==face]
        faces[str(face)]={"incoming_m2_s":float(e["face_flux_m2_s"][e["direction"]==1].sum()),
            "outgoing_m2_s":float(e["face_flux_m2_s"][e["direction"]==-1].sum())}
    residual=abs(particles["residual_eV"])
    closure_scale=np.maximum.reduce([abs(particles["deltaK_eV"]),abs(particles["work_eV"]),np.ones(len(particles))])
    summary={"number_density_m3":float(density),"number_density_cm3":float(density/1e6),
        "unique_probe_particles":len(hit_ids),"probe_effective_particles":float(effective),
        "residence_segments":len(records),"crossing_events":len(events),
        "volume_averaged_number_flux_vector_m2_s":current.tolist(),"face_fluxes":faces,
        "source_total_rate_s1":float(particles["weight_s1"].sum()),
        "source_cell_neff_min_median_max":np.quantile(cell_neff,[0,.5,1]).tolist(),
        "velocity_bin_width_kms":args.dv_kms,"velocity_edges_kms":[float(edge[0]/1e3),float(edge[-1]/1e3)],
        "central_slab_bounds_kms":[float(edge[center]/1e3),float(edge[center+1]/1e3)],
        "max_energy_closure_error_eV":float(residual.max()),
        "energy_closure_relative_p50_p95_max":np.quantile(residual/closure_scale,[.5,.95,1]).tolist(),
        "max_gyro_angle_rad":float(particles["max_gyro_angle_rad"].max()),
        "time_limit_particles":int(np.sum(particles["status"]=="time_limit")),
        "time_limit_weight_fraction":float(particles["weight_s1"][particles["status"]=="time_limit"].sum()/particles["weight_s1"].sum()),
        "status_counts":complete["status_counts"],
        "psd_definition":"sum W_i * dwell_time_i_in_bin / (cube_volume * dvx*dvy*dvz)",
        "projection_definition":"fxy = integral f3d dvz; fxz = integral f3d dvy",
        "flux_model":meta["flux_model"],"proposal_draws_per_cell":meta["per_cell"],
        "positive_rate_particles":int(np.sum(particles["weight_s1"]>0)),
        "plot_limit_kms":args.plot_limit_kms,
        "source_run_directory":str(args.run.resolve()),
        "limitations":"Finite source samples, finite flight-age window, finite cube average. Sampling zeros are not physical zeros. No steady-state or MC convergence claim."}
    if meta["flux_model"]=="reservoir_maxwellian_rate":
        summary["analytic_source_outward_rate_s1"]=float(analytic_rate)
        summary["source_rate_mc_standard_error_s1"]=float(math.sqrt(rate_variance))
        summary["source_rate_mc_minus_analytic_sigma"]=float((particles["weight_s1"].sum()-analytic_rate)/math.sqrt(rate_variance)) if rate_variance>0 else None
    if args.plot_limit_kms is not None:
        lim=args.plot_limit_kms*1e3
        overlap=np.maximum(0,np.minimum(edge[1:],lim)-np.maximum(edge[:-1],-lim))
        summary["displayed_density_fraction_xy"]=float(np.sum(fxy*overlap[:,None]*overlap[None,:])/density) if density>0 else 0.
        summary["displayed_density_fraction_xz"]=float(np.sum(fxz*overlap[:,None]*overlap[None,:])/density) if density>0 else 0.
    np.savez_compressed(out/"probe_psd.npz",vx_edges_ms=edge,vy_edges_ms=edge,vz_edges_ms=edge,
        f3d_s3_m6=f3d,fxy_s2_m5=fxy,fxz_s2_m5=fxz,number_per_bin=number,
        unique_particles_per_bin=unique,effective_particles_per_bin=neff,
        slice_xy_s3_m6=slice_xy,slice_xz_s3_m6=slice_xz)
    (out/"analysis_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf8")
    draw_panels(out/"probe_psd_projections",edges,fxy,fxz,r"Reduced PSD (s$^2$ m$^{-5}$)",
        f"Velocity-integrated projections; {args.dv_kms:g} km/s bins; maximum flight age {meta['max_flight_time_s']:g} s",len(hit_ids),effective,meta)
    draw_panels(out/"probe_psd_slices",edges,slice_xy,slice_xz,r"3D PSD (s$^3$ m$^{-6}$)",
        f"Central velocity slabs [{edge[center]/1e3:g}, {edge[center+1]/1e3:g}] km/s in omitted component",len(hit_ids),effective,meta)
    if args.plot_limit_kms is None:
        draw_panels(out/"probe_psd_projections_zoom",edges,fxy,fxz,r"Reduced PSD (s$^2$ m$^{-5}$)",
            f"Velocity-integrated projections; {args.dv_kms:g} km/s bins; all nonzero bins shown",len(hit_ids),effective,meta,zoom=True)
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
