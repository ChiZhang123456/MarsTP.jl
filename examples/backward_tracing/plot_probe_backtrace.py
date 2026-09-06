"""Run the reproducible probe backtrace and save data, diagnostics and figures."""
import json
import argparse
import os
import subprocess
from pathlib import Path
from datetime import datetime
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from py_space_zc.maven import bs_mpb, plot_mars

ROOT = Path(__file__).resolve().parents[2]
N = int(os.environ.get('PARTICLE_COUNT', '5000'))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data', type=Path, help='Existing trajectories.jsonl to replot without integration')
args = parser.parse_args()
OUT = ROOT / 'outputs' / ('probe_backtrace_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
IMAGE = ROOT / 'examples/backward_tracing/images/probe_backtrace_xz_5000'
IMAGE.parent.mkdir(exist_ok=True)
OUT.mkdir(parents=True, exist_ok=False)
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9,
    'axes.labelsize': 9, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6})
records = []
if args.data:
    records = [json.loads(line) for line in args.data.open()]
else:
    with (OUT/'run.log').open('w') as log, (OUT/'trajectories.jsonl').open('w') as data:
        proc = subprocess.Popen(['julia', '--startup-file=no', f'--project={ROOT}',
            str(ROOT/'examples/backward_tracing/random_probe_backtrace.jl')], stdout=subprocess.PIPE,
            stderr=log, text=True, cwd=ROOT)
        for line in proc.stdout:
            if line.startswith('{'):
                records.append(json.loads(line))
                data.write(line)
            else:
                log.write(line)
        code = proc.wait()
    if code:
        raise RuntimeError((OUT/'run.log').read_text())
assert len(records) == N
counts = Counter(r['status'] for r in records)
for r in records:
    p = np.asarray(r['points'])
    assert np.isfinite(p).all() and np.allclose(p[0], [0,0,2])
    assert 10 <= r['speed_kms'] <= 200 and r['v0_kms'][1] == 0
    assert len(p) == len(r['times_s']) and np.all(np.diff(r['times_s']) < 0)
    radius = np.linalg.norm(p, axis=1)
    assert radius.min() >= 3590/3390-1e-9 and radius.max() <= 4+1e-9
    if r['status'] != 'time_limit':
        assert abs(radius[-1]-(3590/3390 if r['status']=='inner' else 4)) < 1e-8
summary = dict(particle_count=N, counts=counts, seed=20260905, species='O2+', Rm_km=3390,
    detector_Rm=[0,0,2], speed_sampling='uniform 10 to 200 km/s',
    direction_sampling='uniform angle 0 to 2pi in XZ; initial Vy=0',
    electric_field='static total E from VTS', coordinate_system='native MHD Cartesian axes',
    coordinate_note='MSO/MSE designation not established by repository metadata',
    integrator='nonrelativistic Boris, negative time; original positive ion charge',
    dt_s=-0.05, comparison_dt_s=-0.1, limit_s=500, saved_cadence_s=0.5,
    boundary_endpoint='linear segment-sphere intersection',
    max_common_time_separation_Rm=max(r['max_separation_Rm'] for r in records),
    max_endpoint_difference_Rm=max(r['endpoint_difference_Rm'] for r in records),
    max_termination_time_difference_s=max(r['time_difference_s'] for r in records),
    changed_termination_count=sum(r['status']!=r['coarse_status'] for r in records),
    git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    input_file=str(ROOT/'data/mars_fields_spherical_from_dat.vts'))
(OUT/'metadata.json').write_text(json.dumps(summary, indent=2))
np.savetxt(OUT/'initial_velocities_kms.csv',
    [[r['id'],*r['v0_kms'],r['speed_kms']] for r in records],
    delimiter=',',header='id,vx_kms,vy_kms,vz_kms,speed_kms',comments='')

norm = mpl.colors.Normalize(10,200)
cmap = plt.get_cmap('turbo')
colors = [cmap(norm(r['speed_kms'])) for r in records]
fig, ax = plt.subplots(figsize=(183/25.4, 172/25.4), layout='constrained')
bs_mpb(ax=ax, draw_mpb=False, sphere=False, boundary_color='#444444',
    boundary_ls='--', boundary_lw=.8, mars_lw=0, mars_ls='-')
bs_mpb(ax=ax, draw_bs=False, sphere=False, boundary_color='#444444',
    boundary_ls=':', boundary_lw=.9, mars_lw=0, mars_ls='-')
plot_mars(ax=ax, radius=1.0, alpha=.85, zorder=3)
for radius in [3590/3390,4]:
    ax.add_patch(Circle((0,0),radius,fill=False,ls='--',lw=.8,color='#999999'))
paths = [np.asarray(r['points'])[:,[0,2]] for r in records]
ax.add_collection(LineCollection(paths,colors=colors,linewidths=.35,alpha=.28,rasterized=True))
ax.scatter(0,2,marker='*',s=190,c='crimson',edgecolors='white',zorder=5)
ax.set(xlim=(-4.2,4.2),ylim=(-4.2,4.2),aspect='equal',
    xlabel=r'X / $R_M$',ylabel=r'Z / $R_M$')
ax.grid(alpha=.15)
fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,
    shrink=.8,label='Initial speed (km/s)',pad=.025)
ax.set_title(f'{N:,} O$_2^+$ backtraces from (0, 0, 2 $R_M$)\n'
    f'Inner: {counts["inner"]}   |   Outer: {counts["outer"]}   |   '
    f'500 s limit: {counts["time_limit"]}',fontsize=10,pad=9)
fig.savefig(IMAGE.with_suffix('.png'),dpi=400)
print(json.dumps(summary,indent=2))
print('OUTPUT:',OUT,flush=True)
