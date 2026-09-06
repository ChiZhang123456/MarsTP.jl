"""Plot the marginal VDF, not a Vy=0 slice or a trajectory histogram."""
from pathlib import Path
import sys,json
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont
out=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'data'
image_dir=Path(__file__).resolve().parent/'images'
image_dir.mkdir(exist_ok=True)
a=np.genfromtxt(out/'vdf.csv',delimiter=',',names=True)
vx=np.unique(a['vx_kms']);vz=np.unique(a['vz_kms'])
f=a['f_total_s2_m5'].reshape(len(vz),len(vx))
fvol=a['f_volume_s2_m5'].reshape(f.shape);fsheet=a['f_sheet_s2_m5'].reshape(f.shape)
assert np.isfinite(f).all() and (f>=0).all() and np.allclose(f,fvol+fsheet,rtol=1e-12,atol=0)
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':10,'axes.linewidth':.7})
# Quantitative single-panel figure: show the physical source-weighted VDF.
# Native 10 km/s grid, no smoothing. Finite Vy integration range is in the title.
positive=f[f>0]
assert len(positive)
upper=10.**np.ceil(np.log10(positive.max()));lower=max(10.**np.floor(np.log10(positive.min())),upper/1e8)
fig,ax=plt.subplots(figsize=(183/25.4,165/25.4),layout='constrained')
im=ax.pcolormesh(vx,vz,np.ma.masked_equal(f,0),cmap=plt.get_cmap('turbo').with_extremes(bad='#eeeeee'),norm=LogNorm(lower,upper),shading='nearest',rasterized=True)
ax.set(xlim=(vx.min()-5,vx.max()+5),ylim=(vz.min()-5,vz.max()+5),aspect='equal',
    xlabel=r'MSO $v_x$ (km/s)',ylabel=r'MSO $v_z$ (km/s)')
ax.set_title(r'O$_2^+$ at (0, 0, 2 $R_M$), preliminary grid'+'\n'+r'$f_{xz}=\int_{-200}^{200} f\,dv_y$  ($v_y$ in km/s)',fontsize=11,pad=10)
fig.colorbar(im,ax=ax,pad=.025,shrink=.83,extend='min' if positive.min()<lower else 'neither',label=r'$f_{xz}$ (s$^2$ m$^{-5}$)')
ax.axhline(0,color='#444444',alpha=.22,lw=.6);ax.axvline(0,color='#444444',alpha=.22,lw=.6)
fig.savefig(image_dir/'vdf_xz_preliminary.png',dpi=350)
peak=np.unravel_index(np.argmax(f),f.shape)
meta={'definition':'Vy-integrated marginal VDF, not Vy=0 slice','units':'s^2 m^-5','detector_Rm':[0,0,2],
 'vx_vz_range_kms':[-200,200],'dvx_dvz_kms':10,'vy_range_kms':[-200,200],'dvy_kms':20,
 'dt_s':-.1,'max_lookback_s':500,'source_altitude_km':400,'inner_altitude_km':200,
 'color_limits':[lower,upper],'zero_pixels':int((f==0).sum()),'positive_below_color_min':int(((f>0)&(f<lower)).sum()),
 'peak':{'vx_kms':float(vx[peak[1]]),'vz_kms':float(vz[peak[0]]),'value':float(f[peak])},
 'display':'turbo; no smoothing; exact zeros grey; values below color minimum use under-color',
 'finite_grid_number_density_m3':float(f.sum()*1e8),
 'volume_fraction_on_grid':float(fvol.sum()/f.sum()),'sheet_fraction_on_grid':float(fsheet.sum()/f.sum()),
 'convergence':'FAILED sampled Vy-spacing and timestep checks; see qa_refinement.csv', 'note':'finite velocity range and 500 s source history; numerical zeros may include Maxwellian underflow; not a converged density estimate'}
(out/'plot_metadata.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2))
