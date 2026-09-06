"""Full-path forward-time energy gains, PSD-weighted over Vy; PNG only."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.font_manager import findfont
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent/'data'
a=np.genfromtxt(out/'energy_xz.csv',delimiter=',',names=True)
b=np.genfromtxt(out/'energy_3d.csv',delimiter=',',names=True)
vx=np.unique(a['vx_kms']);vz=np.unique(a['vz_kms'])
shape=(len(vz),len(vx));f=a['f_s2_m5'].reshape(shape)
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.6})
values=[a[k].reshape(shape)/1000 for k in ('conv_eV','hall_eV','total_eV')]
limit=max(np.nanmax(abs(v)) for v in values)
fig,axs=plt.subplots(1,3,figsize=(183/25.4,82/25.4),layout='constrained',sharey=True)
for ax,v,title,letter in zip(axs,values,('Convection','Hall','Total'),'abc'):
    im=ax.pcolormesh(vx,vz,np.ma.masked_where((f<=0)|~np.isfinite(v),v),
        cmap=plt.get_cmap('coolwarm').with_extremes(bad='#eeeeee'),
        norm=SymLogNorm(linthresh=.001,linscale=1.,vmin=-limit,vmax=limit,base=10),shading='nearest')
    if np.nanmin(v)<0<np.nanmax(v): ax.contour(vx,vz,v,levels=[0],colors='k',linewidths=.5)
    ax.set(aspect='equal',xlabel=r'MSO $v_x$ (km/s)',title=title)
    ax.text(0,1.06,letter,transform=ax.transAxes,fontweight='bold')
axs[0].set_ylabel(r'MSO $v_z$ (km/s)')
fig.suptitle(r'O$_2^+$: full-path energy gain, PSD-weighted over $v_y$'+'\nPreliminary velocity grid',fontsize=10)
fig.colorbar(im,ax=axs,location='bottom',shrink=.8,pad=.08,label=r'Mean forward-time energy gain (keV)',
    ticks=[-10,-1,-.1,-.01,-.001,0,.001,.01,.1,1,10])
image=Path(__file__).resolve().parent/'images'/'energy_gain_xz_preliminary.png'
fig.savefig(image,dpi=350)
scale=np.maximum.reduce([np.abs(b['total_eV']),np.abs(b['deltaK_eV']),np.ones(len(b))])
bad=abs(b['residual_eV'])/scale>.01
qa={'closure_over_1pct':int(bad.sum()),'closure_over_1pct_positive_psd':int(np.sum(bad & (b['psd_s3_m6']>0))),'trajectories':len(b),'statuses':{str(int(k)):int(np.sum(b['status']==k)) for k in np.unique(b['status'])},
    'max_abs_energy_residual_eV':float(np.max(abs(b['residual_eV']))),
    'relative_residual_scale':'max(abs(total work), abs(delta K), 1 eV)',
    'relative_residual_percentiles':np.percentile(abs(b['residual_eV'])/scale,[50,95,99,100]).tolist(),
    'max_field_sum_residual_eV':float(np.max(abs(b['total_eV']-b['conv_eV']-b['hall_eV']))),
    'color_scale':'SymLogNorm base10; linear within +/-1 eV','color_limit_keV':float(limit),'zero_psd_pixels':int((f<=0).sum()),
    'note':'Unconverged velocity grid. Full-path endpoint-to-detector energy, not birth-to-detector energy. Coolwarm; common symmetric-log limits; no smoothing.'}
(out/'energy_qa.json').write_text(json.dumps(qa,indent=2))
print(json.dumps(qa,indent=2))
