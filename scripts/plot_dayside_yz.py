"""YZ projection of full trajectories selected by dayside initial position."""
import os, json, subprocess
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from py_space_zc.maven import plot_mars
ROOT=Path(__file__).resolve().parents[1]
env=dict(os.environ,RELEASE_ALTITUDE_KM='800',PARTICLE_COUNT='1000',TRACE_DT='0.1')
process=subprocess.Popen(['julia','--startup-file=no','--threads=4',f'--project={ROOT}',
    str(ROOT/'examples/sphere_trajectories.jl')],stdout=subprocess.PIPE,text=True,env=env)
records=[]
total=0
for line in process.stdout:
    if not line.startswith('{'): continue
    record=json.loads(line)
    total+=1
    if record['points'][0][0]>0: records.append(record)
if process.wait(): raise RuntimeError('Trajectory calculation failed')
assert total==1000 and len(records)==500
paths=[np.array(r['points'])[:,[1,2]] for r in records]
assert all(np.isfinite(p).all() for p in paths)
print('Dayside termination counts:',dict(Counter(r['status'] for r in records)),flush=True)
mpl.rcParams.update({'font.family':'Arial','font.size':13})
fig,ax=plt.subplots(figsize=(7.3,7),layout='constrained')
# A geometric disk avoids assigning an unverified surface-map orientation in YZ.
plot_mars(ax=ax,texture=False,facecolor='#dddddd',edgecolor='#888888',lw=.8,alpha=.55,zorder=1)
ax.add_collection(LineCollection(paths,colors='#16749a',linewidths=.65,alpha=.45,zorder=3))
starts=np.array([p[0] for p in paths])
ax.scatter(starts[:,0],starts[:,1],c='#333333',s=4,alpha=.35,zorder=4)
ax.set(xlim=(-4.15,4.15),ylim=(-4.15,4.15),aspect='equal',
       xlabel=r'Y / $R_M$',ylabel=r'Z / $R_M$',title='Dayside O2+')
ax.grid(alpha=.12)
fig.savefig(os.environ['YZ_PREVIEW'],dpi=170)
print('Preview:',os.environ['YZ_PREVIEW'],flush=True)
