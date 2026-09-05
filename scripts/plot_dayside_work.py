"""Signed work on the same dayside trajectories; data are kept in memory."""
import os, json, subprocess
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.collections import LineCollection
from py_space_zc.maven import bs_mpb, plot_mars
ROOT=Path(__file__).resolve().parents[1]
process=subprocess.Popen(['julia','--startup-file=no','--threads=4',f'--project={ROOT}',
    str(ROOT/'scripts/dayside_work.jl')],stdout=subprocess.PIPE,text=True)
records=[json.loads(line) for line in process.stdout if line.startswith('{')]
if process.wait(): raise RuntimeError('Work-profile calculation failed')
assert len(records)==500
summaries=np.array([r['summary'] for r in records])
assert np.isfinite(summaries).all()
print('Counts:',dict(Counter(r['status'] for r in records)),flush=True)
print('Timesteps:',dict(Counter(r['dt'] for r in records)),flush=True)
print('Final eV min/median/max (conv, hall, kinetic, residual):',np.quantile(summaries,[0,.5,1],axis=0).tolist(),flush=True)
print('Negative net work count (conv, hall):',np.sum(summaries[:,:2]<0,axis=0).tolist(),flush=True)
relative=np.abs(summaries[:,3])/np.maximum(np.abs(summaries[:,2]),1)
print('Relative closure error median/max:',np.median(relative),np.max(relative),flush=True)
print('Hall > convection count:',int(np.sum(summaries[:,1]>summaries[:,0])),flush=True)
mpl.rcParams.update({'font.family':'Arial','font.size':11})
fig,axes=plt.subplots(1,3,figsize=(17,6),layout='constrained')
limit=max(np.max(np.abs(np.array(r['points'])[:,2:]/1000)) for r in records)
norm=SymLogNorm(linthresh=.1,vmin=-limit,vmax=limit,base=10)
for ax,column,title in zip(axes[:2],[2,3],['Convection work','Hall work']):
    bs_mpb(ax=ax,sphere=False,mars_lw=0,mars_ls='-',boundary_color='0.4')
    plot_mars(ax=ax,alpha=.5,zorder=1)
    segments,values=[],[]
    for r in records:
        p=np.array(r['points'])
        segments.extend(np.stack([p[:-1,:2],p[1:,:2]],axis=1))
        values.extend((p[:-1,column]+p[1:,column])/2000)
    lines=LineCollection(segments,cmap='RdBu_r',norm=norm,linewidths=.6,alpha=.85,zorder=3)
    lines.set_array(np.array(values)); ax.add_collection(lines)
    ax.set(xlim=(-4.15,4.15),ylim=(-4.15,4.15),aspect='equal',
           xlabel=r'X / $R_M$',ylabel=r'Z / $R_M$',title=title)
fig.colorbar(lines,ax=axes[:2],location='bottom',label='Cumulative signed work (keV), symmetric log scale',shrink=.9,pad=.03)
ax=axes[2]
for status,marker,color,label in [('outer','o','#16749a','Outer boundary'),('inner','^','#b46924','Inner boundary'),('time_limit','s','#9b418e','Time limit')]:
    m=np.array([r['status']==status for r in records])
    if m.any(): ax.scatter(summaries[m,0]/1000,summaries[m,1]/1000,s=18,alpha=.6,marker=marker,c=color,label=label)
ax.axhline(0,c='0.6',lw=.7); ax.axvline(0,c='0.6',lw=.7)
ax.set_xscale('symlog',linthresh=.1); ax.set_yscale('symlog',linthresh=.1)
ax.margins(x=.1,y=.1)
ax.set(xlabel='Final convection work (keV)',ylabel='Final Hall work (keV)',title='Per-particle net work')
ax.legend(frameon=False,fontsize=9)
for ax in axes: ax.grid(alpha=.12)
fig.savefig(os.environ['DAYSIDE_WORK_PREVIEW'],dpi=170)
print('Preview:',os.environ['DAYSIDE_WORK_PREVIEW'],flush=True)
