"""Compare O2+ and O+ flux at 500 km using the existing Julia/Python workflow."""
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.font_manager import findfont

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'outputs'/('outflow_500km_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
OUT.mkdir(parents=True, exist_ok=True)
if not (OUT/'source.csv').exists():
    with (OUT/'source.csv').open('w') as data, (OUT/'run.log').open('w') as log:
        subprocess.run(['julia','--startup-file=no',f'--project={ROOT}',str(Path(__file__).with_name('sample_outflow_500km.jl'))],stdout=data,stderr=log,check=True,cwd=ROOT)
a = np.genfromtxt(OUT/'source.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
assert len(a)==2*91*181
for key in a.dtype.names[1:]:
    assert np.isfinite(a[key]).all(), key
assert np.all(a['n_m3']>=0)
assert np.allclose(a['outward_flux_cm2_s'], a['n_m3']*np.maximum(a['ur_ms'],0)/1e4)
assert np.allclose(a['speed_flux_cm2_s'], a['n_m3']*np.sqrt(a['ux_ms']**2+a['uy_ms']**2+a['uz_ms']**2)/1e4)
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':12,'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
lon = np.unique(a['longitude_deg']); lat = np.unique(a['latitude_deg'])
meta = {'altitude_km':500,'Rm_km':3390,'coordinate':'MSO longitude=atan2(y,x), latitude=asin(z/r); not geographic',
        'input':'data/mars_fields_spherical_from_dat.vts','interpolation':'MarsTP spherical scalar interpolation of density and Cartesian velocity before flux calculation',
        'outward_flux':'n*max(U dot er,0), bulk outward number flux, cm^-2 s^-1',
        'speed_flux':'n*norm(U), cm^-2 s^-1','grid_step_deg':2,'colormap':'turbo','scale':'shared linear, zero to combined maximum',
        'mask':'south pole excluded because the existing model documentation flags its moments as suspect; grey',
        'zero_values':'shown at zero on color scale; inward flow has zero outward flux','statistics':{}}
# Contract: two equally sized quantitative maps compare species under identical
# color limits. No smoothing or percentile clipping. Preserve samples and units.
for field, stem, label in [('outward_flux_cm2_s','outflow_500km',r'Outward flux, $n\max(U_r,0)$ (cm$^{-2}$ s$^{-1}$)'),
                           ('speed_flux_cm2_s','speed_flux_500km',r'Flux magnitude, $n|\mathbf{U}|$ (cm$^{-2}$ s$^{-1}$)')]:
    vals = np.array([a[field][a['species']==s].reshape(len(lat),len(lon)) for s in ('O2+','O+')])
    norm = Normalize(0,float(vals[:,lat>-90,:].max()))
    fig, axes = plt.subplots(1,2,figsize=(10.4,3.8),layout='constrained',sharex=True,sharey=True)
    for i,(ax,title,v) in enumerate(zip(axes,[r'O$_2^+$',r'O$^+$'],vals)):
        im = ax.pcolormesh(lon,lat,np.ma.array(v,mask=np.broadcast_to((lat==-90)[:,None],v.shape)),
                           cmap=plt.get_cmap('turbo').with_extremes(bad='#d9d9d9'),norm=norm,shading='nearest',rasterized=True)
        ax.set(xlim=(-180,180),ylim=(-90,90),xticks=[-180,-90,0,90,180],yticks=[-90,-45,0,45,90],xlabel='MSO longitude (deg)')
        ax.set_title(f'({chr(97+i)}) {title}, 500 km',loc='left')
        ax.set_aspect('equal')
    axes[0].set_ylabel('MSO latitude (deg)')
    cb = fig.colorbar(im, ax=axes,location='bottom',fraction=.1,pad=.07,aspect=45,label=label)
    cb.formatter.set_powerlimits((0,0)); cb.update_ticks()
    for ext in ('png','pdf'):
        fig.savefig(OUT/f'{stem}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    meta['statistics'][field]={'vmin':0,'vmax':norm.vmax,'species_max':[float(v[lat>-90].max()) for v in vals]}
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(OUT,flush=True)
print(json.dumps(meta['statistics'],indent=2))
