"""Replot saved probe trajectories with the current publication style."""
import json
import sys
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
SOURCE = Path(sys.argv[1]).resolve()
OUT = SOURCE / 'turbo_nature'
OUT.mkdir(exist_ok=True)
records = [json.loads(line) for line in (SOURCE/'trajectories.jsonl').open()]
N = len(records)
counts = Counter(r['status'] for r in records)
assert N == 5000
# Contract: show the spatial reach of 5000 fixed-detector backtraces.
# Single quantitative panel; original coordinates and all trajectories retained.
# Turbo is the user-selected initial-speed encoding; no inference statistics.
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9,
    'axes.labelsize': 9, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
norm = mpl.colors.Normalize(10,200)
cmap = plt.get_cmap('turbo')
colors = [cmap(norm(r['speed_kms'])) for r in records]
fig, ax = plt.subplots(figsize=(183/25.4, 172/25.4), layout='constrained')
ax.add_patch(Circle((0,0),1,facecolor='#d9c2a5',edgecolor='#8d7963',zorder=2))
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
fig.savefig(OUT/'trajectories_xz.png',dpi=400)
fig.savefig(OUT/'trajectories_xz.pdf',dpi=600)
fig.savefig(OUT/'trajectories_xz.svg',dpi=600)

print(OUT)
