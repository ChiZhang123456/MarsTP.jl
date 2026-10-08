"""MSO outward bulk flux at 200, 400 and 600 km, six shared-scale panels."""
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
HEIGHTS = (200, 400, 600)
OUT = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'outputs'/('outflow_200_400_600km_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
OUT.mkdir(parents=True, exist_ok=True)
if not (OUT/'source.csv').exists():
    with (OUT/'source.csv').open('w') as data, (OUT/'run.log').open('w') as log:
        subprocess.run(['julia','--startup-file=no',f'--project={ROOT}',str(Path(__file__).with_name('sample_outflow_500km.jl')),*map(str,HEIGHTS)],stdout=data,stderr=log,check=True,cwd=ROOT)
a = np.genfromtxt(OUT/'source.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
assert len(a)==6*91*181
assert all(np.isfinite(a[k]).all() for k in a.dtype.names[1:])
assert (a['n_m3']>=0).all()
assert np.allclose(a['outward_flux_cm2_s'],a['n_m3']*np.maximum(a['ur_ms'],0)/1e4)
lon=np.unique(a['longitude_deg']); lat=np.unique(a['latitude_deg'])
values=np.empty((3,2,len(lat),len(lon)))
for row,h in enumerate(HEIGHTS):
    for col,s in enumerate(('O2+','O+')):
        part=a[(a['altitude_km']==h)&(a['species']==s)]
        assert len(part)==91*181
        assert np.array_equal(part['longitude_deg'].reshape(91,181)[0],lon)
        assert np.array_equal(part['latitude_deg'].reshape(91,181)[:,0],lat)
        values[row,col]=part['outward_flux_cm2_s'].reshape(len(lat),len(lon))
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':12,'axes.linewidth':.7,'pdf.fonttype':42})
norm=Normalize(0,float(values[:,:,lat>-90,:].max()))
# Quantitative grid: compare altitude and species with identical linear limits,
# no smoothing or clipping, and the documented suspect south pole masked.
fig,axes=plt.subplots(3,2,figsize=(10.4,9.1),layout='constrained',sharex=True,sharey=True)
for row,h in enumerate(HEIGHTS):
    for col,title in enumerate((r'O$_2^+$',r'O$^+$')):
        ax=axes[row,col]; v=values[row,col]
        im=ax.pcolormesh(lon,lat,np.ma.array(v,mask=np.broadcast_to((lat==-90)[:,None],v.shape)),
            cmap=plt.get_cmap('turbo').with_extremes(bad='#d9d9d9'),norm=norm,shading='nearest',rasterized=True)
        ax.set(xlim=(-180,180),ylim=(-90,90),xticks=[-180,-90,0,90,180],yticks=[-90,-45,0,45,90],aspect='equal')
        ax.set_title(f'({chr(97+row*2+col)}) {title}, {h} km',loc='left')
        if col==0: ax.set_ylabel('MSO latitude (deg)')
        if row==2: ax.set_xlabel('MSO longitude (deg)')
cb=fig.colorbar(im,ax=axes,location='bottom',fraction=.05,pad=.035,aspect=50,
    label=r'Outward flux, $n\max(U_r,0)$ (cm$^{-2}$ s$^{-1}$)')
cb.formatter.set_powerlimits((0,0)); cb.update_ticks()
for ext in ('png','pdf'):
    fig.savefig(OUT/f'outflow_200_400_600km.{ext}',dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
meta={'altitude_km':HEIGHTS,'Rm_km':3390,'species':['O2+','O+'],
    'coordinate':'MSO longitude=atan2(y,x), latitude=asin(z/r), not geographic',
    'input':'data/mars_fields_spherical_from_dat.vts','flux':'n*max(U dot er,0), cm^-2 s^-1',
    'interpolation':'MarsTP spherical scalar interpolation of density and Cartesian velocity before radial projection; 200 km queries nudged 1e-8 m into domain by ionosphere_properties',
    'grid_step_deg':2,'colormap':'turbo','normalization':'shared linear for all six panels',
    'vmin':0,'vmax':norm.vmax,'panel_maxima':values[:,:,lat>-90,:].max(axis=(2,3)).tolist(),
    'mask':'documented suspect south-pole row shown grey, excluded from limits',
    'zero':'inward and zero radial bulk velocity have zero outward flux; displayed at colorbar minimum'}
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(OUT,flush=True)
print(json.dumps(meta['panel_maxima']))
