"""XZ trajectory segments colored by physical forward-time electric power."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.font_manager import findfont
from py_space_zc.maven import bs_mpb,plot_mars
out=Path(sys.argv[1])
case=json.loads((out/'case.json').read_text()) if (out/'case.json').exists() else {'mode':'original','detector_Rm':[0,0,2],'dt_s':-.05}
probe=case['detector_Rm']
tail=case['mode']!='original'
image_dir=Path(__file__).resolve().parent/'images'
with (out/'segments.csv').open() as stream:
    columns=stream.readline().strip().split(',')
a=np.loadtxt(out/'segments.csv',delimiter=',',skiprows=1,dtype=[(name,'f8') for name in columns],ndmin=1)
b=np.genfromtxt(out/'particles.csv',delimiter=',',names=True)
assert np.isfinite(a.view(float)).all()
if tail:
    closure_scale=np.maximum.reduce([abs(b['total_eV']),abs(b['deltaK_eV']),np.ones(len(b))])
    flagged=abs(b['residual_eV'])/closure_scale>.01
    np.savetxt(out/'refine_points.csv',np.column_stack([b[k][flagged] for k in ('vx_kms','vy_kms','vz_kms')]),delimiter=',',header='vx_kms,vy_kms,vz_kms',comments='')
segments=np.stack([np.column_stack([a['x0_Rm'],a['z0_Rm']]),np.column_stack([a['x1_Rm'],a['z1_Rm']])],axis=1)
ids=a['id'].astype(int)-1
assert np.array_equal(b['id'].astype(int),np.arange(1,len(b)+1))
values=[a[k] for k in ('conv_eV_s','hall_eV_s','total_eV_s')]
# C(t)=integral from the past endpoint to t. Saved lookback sums run from probe.
# Segment color is the average of its cumulative endpoint values.
cumulative=[]
max_endpoint_error=0.
for component,power in zip(('conv','hall','total'),values):
    after=a[f'lookback_{component}_eV'];dw=power*a['elapsed_s']
    cumulative.append(b[f'{component}_eV'][ids]-after+.5*dw)
    ends=np.r_[ids[1:]!=ids[:-1],True]
    max_endpoint_error=max(max_endpoint_error,float(np.max(abs(after[ends]-b[f'{component}_eV']))))
    assert np.allclose(after[ends],b[f'{component}_eV'],rtol=1e-11,atol=1e-8)
    starts=np.r_[True,ids[1:]!=ids[:-1]]
    assert np.allclose((after-dw)[starts],0,atol=1e-8)
limits=[10.**np.ceil(np.log10(max(np.max(abs(v)) for v in row))) for row in (values,cumulative)]
cmap=plt.get_cmap('coolwarm')
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.6})
fig,axs=plt.subplots(2,3,figsize=(183/25.4,177/25.4),layout='constrained',sharey=True)
for row,rowvalues in enumerate((values,cumulative)):
    norm=SymLogNorm(linthresh=1.,vmin=-limits[row],vmax=limits[row],base=10)
    for col,(ax,v,title) in enumerate(zip(axs[row],rowvalues,('Convection','Hall','Total'))):
        bs_mpb(ax=ax,draw_mpb=False,sphere=False,boundary_color='#444444',boundary_ls='--',boundary_lw=.6,mars_lw=0,mars_ls='-')
        bs_mpb(ax=ax,draw_bs=False,sphere=False,boundary_color='#444444',boundary_ls=':',boundary_lw=.6,mars_lw=0,mars_ls='-')
        ax.add_collection(LineCollection(segments,array=v,cmap=cmap,norm=norm,linewidths=.20,alpha=.7,rasterized=True))
        plot_mars(ax=ax,radius=1.,alpha=.95,zorder=3)
        for radius in (3590/3390,4): ax.add_patch(Circle((0,0),radius,fill=False,ls='--',lw=.5,color='#888888'))
        ax.scatter(probe[0],probe[2],marker='*',s=45,c='black',edgecolors='white',linewidths=.7,zorder=5)
        ax.set(xlim=(-4.2,4.2),ylim=(-4.2,4.2),aspect='equal',xlabel=r'X / $R_M$',title=title)
        ax.text(0,1.06,'abcdef'[3*row+col],transform=ax.transAxes,fontweight='bold')
    axs[row,0].set_ylabel(r'Z / $R_M$')
    label=r'Local power (eV/s)' if row==0 else 'Cumulative work from past endpoint (eV)'
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axs[row],location='bottom',shrink=.85,pad=.06,label=label)
fig.suptitle(f'{len(b):,} O$_2^+$ backtraces at ({probe[0]:g}, {probe[1]:g}, {probe[2]:g}) $R_M$'+'\nLocal power and cumulative work',fontsize=10)
image=image_dir/('tail_probe_path_power_xz_5000.png' if tail else 'probe_path_power_xz_5000.png')
fig.savefig(image,dpi=400)
qa={'particles':len(b),'segments':len(a),'seed':20260905,'case':case,'dt_s':case['dt_s'],'display_cadence_s':.5,
    'sign':'Positive forward-time power means local gain, negative means loss; physical velocity is not negated.',
    'display':'XZ projection; segment-averaged power, no spatial binning; overlap and transparency affect appearance.',
    'color_scale':'coolwarm SymLogNorm; linear within +/-1 eV/s (power) or +/-1 eV (work)',
    'symmetric_limit_eV_s':limits[0],'symmetric_limit_eV':limits[1],
    'cumulative_definition':'Integral from past endpoint to current location; zero at past endpoint, full net gain at probe. Segment color averages cumulative endpoints.',
    'cumulative_endpoint_check_max_eV':max_endpoint_error,
    'status_counts':{str(int(k)):int(sum(b['status']==k)) for k in np.unique(b['status'])},
    'max_energy_closure_residual_eV':float(max(abs(b['residual_eV']))),
    'max_power_component_sum_residual_eV_s':float(max(abs(a['total_eV_s']-a['conv_eV_s']-a['hall_eV_s'])))}
(out/'path_power_qa.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))

if tail:
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(183/25.4,172/25.4),layout='constrained')
    bs_mpb(ax=ax,draw_mpb=False,sphere=False,boundary_color='#444444',boundary_ls='--',boundary_lw=.8,mars_lw=0,mars_ls='-')
    bs_mpb(ax=ax,draw_bs=False,sphere=False,boundary_color='#444444',boundary_ls=':',boundary_lw=.8,mars_lw=0,mars_ls='-')
    speed=np.sqrt(sum(b[k]**2 for k in ('vx_kms','vy_kms','vz_kms')))
    snorm=mpl.colors.Normalize(0,100)
    ax.add_collection(LineCollection(segments,array=speed[ids],cmap='turbo',norm=snorm,linewidths=.3,alpha=.35))
    plot_mars(ax=ax,radius=1.,alpha=.95,zorder=3)
    for radius in (3590/3390,4): ax.add_patch(Circle((0,0),radius,fill=False,ls='--',lw=.7,color='#888888'))
    ax.scatter(probe[0],probe[2],marker='*',s=180,c='crimson',edgecolors='white',zorder=5)
    ax.set(xlim=(-4.2,4.2),ylim=(-4.2,4.2),aspect='equal',xlabel=r'X / $R_M$',ylabel=r'Z / $R_M$')
    ax.set_title(f'{len(b):,} O$_2^+$ backtraces from (-2, 0, 0) $R_M$',fontsize=11)
    fig.colorbar(mpl.cm.ScalarMappable(norm=snorm,cmap='turbo'),ax=ax,shrink=.8,label='Initial speed (km/s)')
    fig.savefig(image_dir/'tail_probe_trajectories_xz_5000.png',dpi=400)
