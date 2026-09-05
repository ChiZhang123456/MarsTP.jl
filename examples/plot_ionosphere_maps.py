"""Compare MHD O2+ flux and Ti on 200/400 km shells in native MSO angles."""
import subprocess,json,sys
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont
ROOT=Path(__file__).resolve().parents[1]
IMAGE=ROOT/'examples/images'
IMAGE.mkdir(exist_ok=True)
if len(sys.argv)>1:
    OUT=Path(sys.argv[1]).resolve()
else:
    OUT=ROOT/'outputs'/('ionosphere_maps_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    OUT.mkdir(parents=True,exist_ok=False)
    with (OUT/'source.csv').open('w') as data, (OUT/'run.log').open('w') as log:
        subprocess.run(['julia','--startup-file=no',f'--project={ROOT}',str(ROOT/'examples/sample_ionosphere_maps.jl')],stdout=data,stderr=log,check=True,cwd=ROOT)
a=np.genfromtxt(OUT/'source.csv',delimiter=',',names=True)
assert len(a)==2*91*181 and all(np.isfinite(a[k]).all() for k in a.dtype.names)
assert np.allclose(a['flux_m2_s'],a['n_m3']*np.sqrt(a['ux_ms']**2+a['uy_ms']**2+a['uz_ms']**2),rtol=1e-12)
lon=np.unique(a['longitude_deg']); lat=np.unique(a['latitude_deg'])
flux=a['flux_m2_s'].reshape(2,len(lat),len(lon))/1e4
Ti=a['Ti_K'].reshape(flux.shape)
assert (flux>=0).all() and (Ti>0).all()
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':10,
    'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
# Contract: quantitative grid; compare shell structure with shared logarithmic
# scales per quantity, preserving all positive values without clipping.
fig,axes=plt.subplots(2,2,figsize=(183/25.4,145/25.4),layout='constrained',sharex=True,sharey=True)
for col,(values,label,title) in enumerate([(flux,r'$n_i|\mathbf{U}_i|$ (cm$^{-2}$ s$^{-1}$)','O$_2^+$ flux'),(Ti,r'$T_i$ (K)','O$_2^+$ ion temperature')]):
    valid=values[:,lat>-90,:]
    norm=LogNorm(vmin=10**np.floor(np.log10(valid[valid>0].min())),vmax=10**np.ceil(np.log10(valid.max())))
    for row,h in enumerate([200,400]):
        ax=axes[row,col]
        im=ax.pcolormesh(lon,lat,np.ma.array(values[row],mask=(values[row]==0)|np.broadcast_to((lat==-90)[:,None],values[row].shape)),norm=norm,cmap=plt.get_cmap('turbo').with_extremes(bad='#dddddd'),shading='nearest',rasterized=True)
        ax.set(xlim=(-180,180),ylim=(-90,90),xticks=[-180,-90,0,90,180],yticks=[-90,-45,0,45,90])
        ax.set_title(f'({chr(97+2*row+col)}) {h} km'+(f' | {title}' if row==0 else ''),loc='left')
        if col==0: ax.set_ylabel('MSO latitude (deg)')
        if row==1: ax.set_xlabel('MSO longitude (deg)')
        ax.tick_params(direction='out',width=.6,length=3)
    fig.colorbar(im,ax=axes[:,col].tolist(),location='bottom',shrink=.9,pad=.06,label=label)
fig.savefig(IMAGE/'flux_Ti_200_400km.png',dpi=350)
fig.savefig(IMAGE/'flux_Ti_200_400km.pdf',dpi=600)
fig.savefig(IMAGE/'flux_Ti_200_400km.svg',dpi=600)
svg=IMAGE/'flux_Ti_200_400km.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
meta={'altitude_km':[200,400],'Rm_km':3390,'grid_deg':2,'coordinate':'MSO spherical angles, latitude=90-colatitude; longitude=atan2(y,x)',
    'flux_definition':'n*norm(Ui), no radial projection','flux_plot_unit':'cm^-2 s^-1','temperature_plot_unit':'K',
    'zero_flux_count':[int((x==0).sum()) for x in flux], 'zero_flux_color':'light grey', 'masked_latitude_deg':-90, 'mask_reason':'suspect south-pole moments, temperature near 1e-10 K and zero bulk speed; no interpolation repair', 'color_scale_excludes_masked_pole':True, 'scale':'shared logarithmic per column; no clipping','flux_range_cm2_s':[[float(x.min()),float(x.max())] for x in flux],
    'Ti_range_K':[[float(x.min()),float(x.max())] for x in Ti],
    'source':'data/mars_fields_spherical_from_dat.vts; MarsTP ionosphere_properties',
    'interpolation':'linear spherical-grid scalar interpolation of n, Ti and each Cartesian Ui component; flux computed afterwards',
    'note':'native angular coordinates, not planetographic geography; 200 km query nudged 1e-8 m into table for roundoff'}
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
print(json.dumps(meta,indent=2));print('OUTPUT:',OUT,flush=True)
