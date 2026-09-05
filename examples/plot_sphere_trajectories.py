"""Run Julia, keep trajectories in memory, and preview their XZ projection."""
import json
import os
from pathlib import Path
import subprocess
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from py_space_zc.maven import bs_mpb, plot_mars

ROOT = Path(__file__).resolve().parents[1]
PREVIEW = Path(os.environ.get('TRAJECTORY_PREVIEW', ROOT / 'examples/images/trajectories_800km.png'))
ALTITUDE = float(os.environ.get('RELEASE_ALTITUDE_KM', '800'))
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 11})
process = subprocess.Popen(
    ['julia', '--startup-file=no', '--threads=4', f'--project={ROOT}',
     str(ROOT / 'examples' / 'sphere_trajectories.jl')],
    stdout=subprocess.PIPE, text=True, encoding='utf-8', cwd=ROOT)
records = [json.loads(line) for line in process.stdout if line.startswith('{')]
if process.wait():
    raise RuntimeError('Julia trajectory calculation failed')
assert len(records) == int(os.environ.get('PARTICLE_COUNT', '1000'))
counts = Counter(r['status'] for r in records)
print('Termination counts:', dict(counts), flush=True)
print('Maximum flight time (s):', max(r['time'] for r in records), flush=True)
for side, test in [('X>0', lambda x: x > 0), ('X<=0', lambda x: x <= 0)]:
    print(side, dict(Counter(r['status'] for r in records if test(r['points'][0][0]))), flush=True)
assert not counts['nonfinite'], 'Nonfinite trajectory detected'
for r in records:
    p = np.array(r['points'])
    assert np.isfinite(p).all()
    if r['status'] in ('inner', 'outer'):
        target = (3390 + 200) / 3390 if r['status'] == 'inner' else 4
        assert abs(np.linalg.norm(p[-1]) - target) < 1e-8

fig, axes = plt.subplots(1, 3, figsize=(18, 6.5), layout='constrained')
colors = {'inner': '#b46924', 'outer': '#16749a', 'time_limit': '#9b418e'}
groups = [records,
          [r for r in records if r['points'][0][0] > 0],
          [r for r in records if r['points'][0][0] <= 0]]
assert len(groups[1]) + len(groups[2]) == len(records)
for ax, selected, title in zip(axes, groups,
        ['All O2+', 'Dayside O2+', 'Nightside O2+']):
    lim = 4.15
    bs_mpb(ax=ax, draw_mpb=False, sphere=False, boundary_color='#454545',
           boundary_ls='--', boundary_lw=1, mars_lw=0, mars_ls='-')
    bs_mpb(ax=ax, draw_bs=False, sphere=False, boundary_color='#454545',
           boundary_ls=':', boundary_lw=1.2, mars_lw=0, mars_ls='-')
    plot_mars(ax=ax, zorder=1, alpha=.6)
    for status, color in colors.items():
        paths = [np.array(r['points'])[:, [0, 2]] for r in selected if r['status'] == status]
        ax.add_collection(LineCollection(paths, colors=color, linewidths=.45, alpha=.32, zorder=3))
    starts = np.array([r['points'][0] for r in selected])
    ax.scatter(starts[:, 0], starts[:, 2], s=3, c='#333333', alpha=.25, zorder=4)
    angle = np.linspace(0, 2*np.pi, 500)
    radius = (3390 + ALTITUDE) / 3390
    ax.plot(radius*np.cos(angle), radius*np.sin(angle), color='#999999', lw=.7)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect='equal',
           xlabel=r'X / $R_M$', ylabel=r'Z / $R_M$', title=title)
    ax.grid(alpha=.12)
handles = [Line2D([], [], color=c, label=f'{s.replace("_", " ")}: {counts[s]}')
           for s, c in colors.items() if counts[s]]
handles += [Line2D([], [], color='#454545', ls='--', label='BS'),
            Line2D([], [], color='#454545', ls=':', label='MPB')]
fig.legend(handles=handles, loc='outside lower center', ncol=len(handles), frameon=False)
fig.savefig(PREVIEW, dpi=160)
print('Preview:', PREVIEW, flush=True)
