"""Day/night work maps: shared scales within rows, independent scales between rows."""
import os, json, subprocess
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.collections import LineCollection
import sys
_VIS_ROOT = next(p for p in Path(__file__).resolve().parents if (p / "src" / "MarsTP.jl").exists())
sys.path.insert(0, str(_VIS_ROOT / "src"))
from visualization.mars import bs_mpb, plot_mars
ROOT=Path(__file__).resolve().parents[2]
PREVIEW=Path(os.environ.get('HEMISPHERE_WORK_PREVIEW', ROOT/'examples/forward_tracing/images/electric_work_800km.png'))
env=dict(os.environ,WORK_HEMISPHERE='all')
process=subprocess.Popen(['julia','--startup-file=no','--threads=4',f'--project={ROOT}',
    str(ROOT/'examples/forward_tracing/hemisphere_work.jl')],stdout=subprocess.PIPE,text=True,env=env)
records=[json.loads(line) for line in process.stdout if line.startswith('{')]
if process.wait(): raise RuntimeError('Particle work calculation failed')
assert len(records)==1000
groups=[[r for r in records if r['points'][0][0]>0],
        [r for r in records if r['points'][0][0]<=0]]
assert [len(g) for g in groups]==[500,500]
mpl.rcParams.update({'font.family':'Arial','font.size':12})
fig,axes=plt.subplots(2,3,figsize=(16,10),layout='constrained')
for row,(group,name) in enumerate(zip(groups,['Dayside','Nightside'])):
    arrays=[np.array(r['points']) for r in group]
    assert all(np.isfinite(a).all() for a in arrays)
    limit=max(np.max(np.abs(a[:,2:5])) for a in arrays)/1000
    norm=SymLogNorm(linthresh=.1,vmin=-limit,vmax=limit,base=10)
    summary=np.array([r['summary'] for r in group])
    closure=np.abs(summary[:,3])/np.maximum(np.abs(summary[:,2]),1)
    print(name,'clim (keV):',(-limit,limit),'counts:',dict(Counter(r['status'] for r in group)),
          'max relative closure:',closure.max(),flush=True)
    assert closure.max()<=1e-3
    segments=np.concatenate([np.stack([a[:-1,:2],a[1:,:2]],axis=1) for a in arrays])
    for col,title in enumerate(['Convection work','Hall work','Total work']):
        ax=axes[row,col]
        bs_mpb(ax=ax,sphere=False,mars_lw=0,mars_ls='-',boundary_color='0.4')
        plot_mars(ax=ax,alpha=.5,zorder=1)
        values=np.concatenate([(a[:-1,col+2]+a[1:,col+2])/2000 for a in arrays])
        lines=LineCollection(segments,cmap='RdBu_r',norm=norm,linewidths=.6,alpha=.85,zorder=3)
        lines.set_array(values); ax.add_collection(lines)
        ax.set(xlim=(-4.15,4.15),ylim=(-4.15,4.15),aspect='equal',
               xlabel=r'X / $R_M$',ylabel=(name+'\n' if col==0 else '')+r'Z / $R_M$',
               title=title if row==0 else '')
        ax.grid(alpha=.12)
    bar=fig.colorbar(lines,ax=axes[row,:].tolist(),location='right',fraction=.025,pad=.02,
        label='Cumulative work (keV)')
    ticks=sorted(set([-limit,-10.,-1.,-.1,0.,.1,1.,10.,limit]))
    ticks=[v for v in ticks if -limit<=v<=limit]
    bar.set_ticks(ticks)
    bar.set_ticklabels([f'{v:.3g}' for v in ticks])
fig.savefig(PREVIEW,dpi=170)
print('Preview:',PREVIEW,flush=True)
