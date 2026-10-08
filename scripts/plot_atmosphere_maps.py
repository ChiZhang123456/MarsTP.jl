"""Plot prepared GITM/AMPS inputs and the explicitly marked MarsTP extension.

Figure contract: compare vertical and horizontal structure, without treating the
hot-O population as GITM thermal O or inventing an AMPS temperature. Two-row
quantitative grids, one per altitude; PNG previews and editable PDF/SVG exports.
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

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs' / ('atmosphere_maps_' + datetime.now().strftime('%Y%m%d_%H%M%S'))
OUT.mkdir(parents=True, exist_ok=False)
mpl.rcParams.update({'font.family': 'Arial', 'font.size': 10,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False,
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
        val[np.flatnonzero(extended)[exponent > 20]] = 0.0
    return val.reshape(h.shape)

heights = np.unique(np.r_[np.arange(100,221,5),np.arange(225,601,5)])
profiles = {k:evaluate(k,heights*1e3,0,0) for k in interps}
np.savetxt(OUT/'profiles_lon0_lat0.csv',np.column_stack([heights]+[profiles[k] for k in interps]),
    delimiter=',',header='altitude_km,nCO2_m-3,nO_m-3,Tn_K,nO_hot_m-3',comments='')
colors = {'nCO2':'#247b98','nO':'#bb6b25','nO_hot':'#79549e'}
labels = {'nCO2':r'GITM CO$_2$','nO':'GITM O','nO_hot':'AMPS hot O'}
maps = {}
for height in (150,500):
    maps[height] = {}
    for field in interps:
        d = data['amps' if field=='nO_hot' else 'gitm']
        lon,lat = np.meshgrid(d['longitude'].ravel(),d['latitude'].ravel())
        val = evaluate(field,height*1e3,lat,lon)
        maps[height][field] = (lon,lat,val)
        np.savetxt(OUT/f'map_{height}km_{field}.csv',np.column_stack((lon.ravel(),lat.ravel(),val.ravel())),
            delimiter=',',header='longitude_deg,latitude_deg,'+field+('_K' if field=='Tn' else '_m-3'),comments='')

for height in (150,500):
    fig = plt.figure(figsize=(16,8.4),layout='constrained')
    gs = fig.add_gridspec(2,4,height_ratios=(1.15,1))
    axn,axt = fig.add_subplot(gs[0,:2]),fig.add_subplot(gs[0,2:])
    for field in colors:
        y=profiles[field]/1e6
        if field=='nO_hot':
            axn.plot(y,heights,color=colors[field],label=labels[field],lw=2)
        else:
            for mask,style,label in [(heights<=220,'-',labels[field]),(heights>=220,'--',None)]:
                axn.plot(np.where(y[mask]>0,y[mask],np.nan),heights[mask],style,color=colors[field],label=label,lw=2)
    for mask,style in [(heights<=220,'-'),(heights>=220,'--')]:
        axt.plot(profiles['Tn'][mask],heights[mask],style,color='#247b98',lw=2,
            label=r'GITM neutral temperature (CO$_2$ and O)' if style=='-' else None)
    axn.set_xscale('log')
    axn.set_xlabel('Number density (cm$^{-3}$)')
    axt.set_xlabel('Temperature (K)')
    for ax,title in [(axn,'a  Density profile'),(axt,'b  Temperature profile')]:
        ax.set_title(title+'  |  longitude = latitude = 0°',loc='left',fontweight='bold')
        ax.set_ylabel('Altitude (km)');ax.set_ylim(100,600)
        ax.axhspan(220,600,color='#ededed',zorder=-5)
        ax.axhline(height,color='#444444',lw=.8,ls=':')
        ax.grid(alpha=.18);ax.legend(loc='lower right',fontsize=10)
    axt.text(.03,.95,'Above 220 km: MarsTP extension (dashed)\nAMPS temperature: not available',
        transform=axt.transAxes,va='top',fontsize=10)
    fields = ['nCO2','nO','nO_hot','Tn']
    for i,field in enumerate(fields):
        ax = fig.add_subplot(gs[1,i]);lon,lat,val=maps[height][field]
        density=field!='Tn'
        z=val/1e6 if density else val
        allz=np.concatenate([maps[h][field][2].ravel()/(1e6 if density else 1) for h in (150,500)])
        pos=allz[allz>0]
        norm=LogNorm(pos.min(),pos.max()) if density else Normalize(pos.min(),pos.max())
        cmap=mpl.colormaps['viridis' if density else 'inferno'].copy();cmap.set_bad('#d9d9d9')
        mesh=ax.pcolormesh(lon,lat,np.ma.masked_less_equal(z,0),shading='nearest',cmap=cmap,norm=norm,rasterized=True)
        title=labels.get(field,'GITM neutral T')
        if height>220 and field!='nO_hot':title+=' (extended)'
        ax.set_title(chr(99+i)+'  '+title,loc='left',fontsize=10,fontweight='bold')
        ax.set_xlim(0,360);ax.set_ylim(-90,90);ax.set_xticks([0,90,180,270,360]);ax.set_yticks([-90,-45,0,45,90])
        ax.set_xlabel('Longitude (°)');ax.set_ylabel('Latitude (°)' if i==0 else '')
        cb=fig.colorbar(mesh,ax=ax,orientation='horizontal',pad=.08,fraction=.07,aspect=28)
        cb.set_label('Number density (cm$^{-3}$)' if density else 'Temperature (K)')
    fig.suptitle(f'GITM and AMPS atmosphere  |  maps at {height} km',fontsize=17,fontweight='bold')
    fig.supxlabel('Altitude = r − 3390 km. File longitude/latitude coordinates. Shared color scales between heights.\n'
        'GITM above 220 km: constant-T exponential extension; grey map cells: model cutoff (Δh/H > 20). AMPS 150 km: linear interpolation.',fontsize=10)
    for ext in ('png','pdf','svg'):
        fig.savefig(OUT/f'atmosphere_{height}km.{ext}',dpi=200)
    plt.close(fig)

metadata={'purpose':'Vertical and horizontal GITM/AMPS comparison, not a fit or statistical inference',
 'radius_m':RM,'profile_altitude_km':[100,600], 'map_altitudes_km':[150,500],
 'units':'Input SI as consumed by MarsTP chemistry; densities plotted in cm^-3; T in K.',
 'coordinates':'File longitude/latitude. No unsupported geographic/MSO interpretation imposed.',
 'interpolation':'Linear in altitude, latitude and longitude; no out-of-grid interpolation.',
 'GITM_extension':'Matches src/data/atmosphere.jl prescription: reference at 219995 m, dh=h-220000 m, constant g=3.71, species mass = 15.999 or 44.01 times proton mass; exp(-dh/H), zero for dh/H>20; T held constant.',
 'AMPS':'Only nO_hot provided; no temperature inferred. File values at 100 km retained (not tracing boundary zero). Longitude endpoints retained separately, not forced periodic.',
 'QA':qa,'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':mpl.__version__},
 'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
 'input_sha256':{str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').glob('*sph.mat')}}
(OUT/'metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
print(OUT)
print(json.dumps(qa,indent=2))
