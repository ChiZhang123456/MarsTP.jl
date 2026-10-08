"""Display trajectory and signed field-work diagnostics, without saving data."""
import io
import os
import subprocess
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from py_space_zc.maven import bs_mpb, plot_mars
ROOT = Path(__file__).resolve().parents[1]
run = subprocess.run(['julia','--startup-file=no',f'--project={ROOT}',
    str(ROOT/'scripts/single_particle_work.jl')],stdout=subprocess.PIPE,text=True,check=True)
d = np.genfromtxt(io.StringIO(run.stdout),delimiter=',',names=True)
mpl.rcParams.update({'font.family':'Arial','font.size':11})
fig,axes = plt.subplots(1,3,figsize=(16,5),layout='constrained')
ax = axes[0]
bs_mpb(ax=ax,sphere=False,mars_lw=0,mars_ls='-',boundary_color='0.4')
plot_mars(ax=ax,alpha=.55,zorder=1)
ax.plot(d['x'],d['z'],color='#16749a',lw=1.7,zorder=3)
ax.scatter(d['x'][0],d['z'][0],s=40,c='#339354',label='Start',zorder=4)
ax.scatter(d['x'][-1],d['z'][-1],s=40,c='#bd432c',label='End',zorder=4)
ax.set(xlim=(-2,2),ylim=(-.5,4.2),aspect='equal',xlabel=r'X / $R_M$',ylabel=r'Z / $R_M$',title='Dayside O2+')
ax.legend(frameon=False,fontsize=9)
ax = axes[1]
ax.plot(d['t'],d['Wconv'],label='Convection work',color='#16749a')
ax.plot(d['t'],d['Whall'],label='Hall work',color='#b46924')
ax.plot(d['t'],d['Wtotal'],label='Total work',color='#9560a8',lw=2.5,alpha=.6)
ax.plot(d['t'],d['K']-d['K'][0],label='Kinetic energy change',color='black',ls='--',lw=1)
ax.axhline(0,color='0.6',lw=.6)
ax.set(xlabel='Time (s)',ylabel='Energy (eV)',title='Cumulative work')
ax.legend(frameon=False,fontsize=9)
ax=axes[2]
ax.plot(d['t'],d['Pconv'],label='Convection',color='#16749a')
ax.plot(d['t'],d['Phall'],label='Hall',color='#b46924')
ax.axhline(0,color='0.6',lw=.6)
ax.set(xlabel='Time (s)',ylabel='Power (eV/s)',title='Instantaneous power')
ax.legend(frameon=False,fontsize=9)
for ax in axes: ax.grid(alpha=.15)
fig.savefig(os.environ['WORK_PREVIEW'],dpi=170)
print('Final values:', {n:float(d[n][-1]) for n in d.dtype.names})
print('Power extrema:', {n:(float(d[n].min()),float(d[n].max())) for n in ['Pconv','Phall']})
