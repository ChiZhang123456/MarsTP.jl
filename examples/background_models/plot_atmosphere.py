"""Plot prepared GITM/AMPS inputs and the explicitly marked MarsTP extension.

Figure contract: compare vertical and horizontal structure, without treating the
hot-O population as GITM thermal O or inventing an AMPS temperature. Three-row, two-column
quantitative grid with GITM at 200 km and AMPS at 500 km; PNG output only.
"""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import subprocess
import numpy as np
import scipy
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.font_manager import findfont

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)
(OUT/'images').mkdir(exist_ok=True)
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 10,
    'axes.spines.top': False,
    'axes.spines.right': False, 'legend.frameon': False})
findfont('Arial', fallback_to_default=False)
RM = 3390e3
KB, MP, GRAV = 1.380649e-23, 1.67262192369e-27, 3.71
data = {name: loadmat(ROOT / 'data' / f'{name}_sph.mat') for name in ('gitm', 'amps')}
interps = {}
qa = {}
for model, d in data.items():
    h, lat, lon = [d[k].ravel().astype(float) for k in ('altitude', 'latitude', 'longitude')]
    assert np.allclose(d['r'].ravel() - h, RM)
    assert np.allclose(d['theta'].ravel(), np.deg2rad(90-lat))
    assert np.allclose(d['phi'].ravel(), np.deg2rad(lon))
    assert all(np.all(np.diff(a) > 0) for a in (h, lat, lon))
    for field in (('nCO2', 'nO', 'Tn') if model == 'gitm' else ('nO_hot',)):
        a = d[field]
        assert a.shape == (len(h),len(lat),len(lon))
        assert np.isfinite(a).all() and (a > 0).all()
        interps[field] = RegularGridInterpolator((h, lat, lon), a, bounds_error=True)
        ix = (len(h)//2, len(lat)//2, len(lon)//2)
        assert np.isclose(interps[field]([[h[ix[0]],lat[ix[1]],lon[ix[2]]]])[0],a[ix])
        qa[field] = {'shape': list(a.shape), 'range': [float(a.min()),float(a.max())],
            'longitude_endpoint_max_relative_difference': float(np.max(np.abs(a[...,0]-a[...,-1])/a[...,0]))}

def evaluate(field, h, lat, lon):
    h, lat, lon = np.broadcast_arrays(h,lat,lon)
    pts = np.column_stack((h.ravel(),lat.ravel(),lon.ravel()))
    if field == 'nO_hot':
        return interps[field](pts).reshape(h.shape)
    extended = pts[:,0] > 220e3
    base = pts.copy()
    base[extended,0] = 220e3-5.0  # Exact reference height in neutral_properties.
    val = interps[field](base)
    if field != 'Tn' and extended.any():
        temp = interps['Tn'](base[extended])
        mass = (44.01 if field == 'nCO2' else 15.999)*MP
        exponent = (pts[extended,0]-220e3)/(KB*temp/(mass*GRAV))
        val[extended] *= np.exp(-exponent)
        # Retain the full exponential tail, with no artificial zero cutoff.
    return val.reshape(h.shape)


# Comparison figure, not a statistical inference. Density maps use independent
# logarithmic scales; temperature uses a linear scale. All maps use turbo.
heights = np.unique(np.r_[np.arange(100,221,5),np.arange(225,601,5)])
profiles = {k:evaluate(k,heights*1e3,0,0) for k in interps}
colors = {'nCO2':'#247b98','nO':'#bb6b25','nO_hot':'#79549e'}
labels = {'nCO2':r'GITM CO$_2$','nO':'GITM O','nO_hot':'AMPS hot O','Tn':'GITM neutral temperature'}
fig, axes = plt.subplots(3,2,figsize=(9.2,10.5),layout='constrained')
axn,axt=axes[0]
for field in colors:
    for mask,style in ([(heights>=100,'-')] if field=='nO_hot' else [(heights<=220,'-'),(heights>=220,'--')]):
        axn.plot(profiles[field][mask]/1e6,heights[mask],style,color=colors[field],lw=1.6,
                 label=labels[field] if style=='-' else None)
for mask,style in [(heights<=220,'-'),(heights>=220,'--')]:
    axt.plot(profiles['Tn'][mask],heights[mask],style,color=colors['nCO2'],lw=1.6)
axn.set_xscale('log');axn.set_xlabel(r'Number density (cm$^{-3}$)')
axt.set_xlabel('Temperature (K)')
axn.legend(loc='upper right',fontsize=8)
for ax,title in [(axn,'a  Density'),(axt,'b  Neutral temperature')]:
    ax.set_title(title,loc='left',fontweight='bold')
    ax.set_ylabel('Altitude (km)');ax.set_ylim(100,600)
    ax.axhline(220,color='0.65',ls=':',lw=.8)
    ax.grid(alpha=.15)
    ax.text(.03,.97,'Longitude = latitude = 0°',transform=ax.transAxes,va='top',fontsize=8)
axt.text(.04,.57,'GITM: dashed above 220 km\nConstant temperature extension\nAMPS temperature unavailable',
         transform=axt.transAxes,va='top',fontsize=8)
# Keep density legend away from the longitude annotation.
axn.legend(loc='lower left',fontsize=8)
map_data={}
for ax,field,letter in zip(axes[1:].flat,['nCO2','nO','Tn','nO_hot'],'cdef'):
    d=data['amps' if field=='nO_hot' else 'gitm']
    lon,lat=np.meshgrid(d['longitude'].ravel(),d['latitude'].ravel())
    map_height = 500 if field=='nO_hot' else 200
    z=evaluate(field,map_height*1e3,lat,lon)/(1 if field=='Tn' else 1e6)
    assert np.isfinite(z).all() and (z>0).all(), field
    map_data[field]=z
    qa[field]['map_altitude_km']=map_height
    if field!='nO_hot':
        native=d[field][np.flatnonzero(d['altitude'].ravel()==200e3)[0]]/(1 if field=='Tn' else 1e6)
        np.testing.assert_allclose(z,native,rtol=1e-13)
        qa[field]['native_200km_slice_verified']=True
    density=field!='Tn'
    norm=LogNorm(z.min(),z.max()) if density else Normalize(z.min(),z.max())
    mesh=ax.pcolormesh(lon,lat,z,shading='nearest',cmap='turbo',norm=norm,rasterized=True)
    if density:
        levels=10.0**np.arange(np.ceil(np.log10(z.min())),np.floor(np.log10(z.max()))+1)
        if len(levels)<3:levels=np.geomspace(z.min(),z.max(),6)[1:-1]
        if len(levels)>6:levels=levels[::int(np.ceil(len(levels)/6))]
        fmt=lambda v: f'{v:.0e}' if v<.01 or v>=1e4 else f'{v:.2g}'
    else:
        levels=np.arange(np.ceil(z.min()/50)*50,z.max(),50)
        fmt='%d'
    cs=ax.contour(lon,lat,z,levels=levels,colors='black',linewidths=.55,alpha=.8)
    ax.clabel(cs,inline=True,fontsize=7,fmt=fmt)
    ax.set_title(f'{letter}  {labels[field]} ({map_height} km)',loc='left',fontweight='bold')
    ax.set(xlim=(0,360),ylim=(-90,90),xlabel='Longitude (°)',ylabel='Latitude (°)',
           xticks=[0,90,180,270,360],yticks=[-90,-45,0,45,90])
    cb=fig.colorbar(mesh,ax=ax,pad=.025,fraction=.045,aspect=24)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(r'Density (cm$^{-3}$)' if density else 'Temperature (K)',fontsize=9)
    qa[field]['map_min_max']=[float(z.min()),float(z.max())]
    qa[field]['map_zero_count']=int(np.sum(z==0))
fig.suptitle('GITM and AMPS atmosphere',fontsize=15,fontweight='bold')
fig.supxlabel('Maps: GITM 200 km (native), AMPS 500 km; Mars radius: 3390 km. Profiles above 220 km:\n'
              'GITM exponential density and constant temperature extension, no zero cutoff.',fontsize=8)
fig.savefig(OUT/'images'/'atmosphere_gitm200km_amps500km.png',dpi=300)
plt.close(fig)
np.savez_compressed(OUT/'atmosphere_source_data.npz',altitude_km=heights,
                    **{'profile_'+k:v for k,v in profiles.items()},
                    **{'map_'+k:v for k,v in map_data.items()},
                    gitm_latitude=data['gitm']['latitude'],gitm_longitude=data['gitm']['longitude'],
                    amps_latitude=data['amps']['latitude'],amps_longitude=data['amps']['longitude'])
metadata={'radius_m':RM,'map_altitude_km':{'gitm':200,'amps':500},'zero_cutoff':False,
          'extension':'n=n(219995 m)*exp(-(h-220000 m)/H); H=kB*T/(m*g); T=T(219995 m).',
          'mass_convention':'15.999 and 44.01 times proton mass, as in MarsTP',
          'gravity_m_s2':GRAV,'input_density_unit':'m^-3','plot_density_unit':'cm^-3',
          'AMPS':'No temperature. Retain original longitude endpoints without enforcing periodicity.',
          'qa':qa,'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':mpl.__version__},
          'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').glob('*sph.mat')}}
(OUT/'metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
print(OUT)
print(json.dumps(qa,indent=2))
