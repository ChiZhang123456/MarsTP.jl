"""Three altitudes by four columns, separate shared scales for speed/outward flux."""
import json
import sys
from datetime import datetime
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(sys.argv[1]).resolve()
OUT = ROOT/'outputs'/('speed_outward_log_200_400_600km_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
OUT.mkdir(parents=True,exist_ok=False)
a=np.genfromtxt(SOURCE/'source.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
assert len(a)==6*91*181
assert all(np.isfinite(a[k]).all() for k in a.dtype.names[1:])
assert (a['n_m3']>=0).all()
assert np.allclose(a['speed_flux_cm2_s'],a['n_m3']*np.sqrt(a['ux_ms']**2+a['uy_ms']**2+a['uz_ms']**2)/1e4)
lon=np.unique(a['longitude_deg']); lat=np.unique(a['latitude_deg'])
heights=(200,400,600)
speed=np.empty((3,2,len(lat),len(lon)))
radial=np.empty_like(speed); vr=np.empty_like(speed)
for row,h in enumerate(heights):
    for col,s in enumerate(('O2+','O+')):
        p=a[(a['altitude_km']==h)&(a['species']==s)]
        assert len(p)==91*181
        assert np.array_equal(p['longitude_deg'].reshape(91,181),np.broadcast_to(lon,(91,181)))
        assert np.array_equal(p['latitude_deg'].reshape(91,181),np.broadcast_to(lat[:,None],(91,181)))
        speed[row,col]=p['speed_flux_cm2_s'].reshape(91,181)
        vr[row,col]=p['ur_ms'].reshape(91,181)
        radial[row,col]=(p['n_m3']*p['ur_ms']/1e4).reshape(91,181)
pole=np.broadcast_to((lat==-90)[None,None,:,None],speed.shape)
speed=np.ma.array(speed,mask=pole|(speed<=0))
outward=np.ma.array(radial,mask=pole|(vr<=0)|(radial<=0))
assert outward.count()>0 and (outward.compressed()>=0).all()
norms=[LogNorm(float(v.min()),float(v.max())) for v in (speed,outward)]
findfont('Arial',fallback_to_default=False)
mpl.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.linewidth':.7,'pdf.fonttype':42})
# Quantitative grid: left and right groups each use one shared logarithmic scale
# across species and altitudes; only positive radial velocities contribute to the right group.
fig=plt.figure(figsize=(18,9),layout='constrained')
subfigs=fig.subfigures(1,2,wspace=.04)
for group,(sub,vals,norm,label) in enumerate(zip(subfigs,(speed,outward),norms,
        (r'$n|\mathbf{V}|$ (cm$^{-2}$ s$^{-1}$)',r'$nV_r$, $V_r>0$ (cm$^{-2}$ s$^{-1}$)'))):
    sub.suptitle(r'Flux magnitude, $n|\mathbf{V}|$' if group==0 else r'Outward radial flux, $nV_r$ ($V_r>0$)',fontsize=14)
    axes=sub.subplots(3,2,sharex=True,sharey=True)
    for row,h in enumerate(heights):
        for col,title in enumerate((r'O$_2^+$',r'O$^+$')):
            ax=axes[row,col]
            im=ax.pcolormesh(lon,lat,vals[row,col],norm=norm,
                cmap=plt.get_cmap('turbo').with_extremes(bad='white'),shading='nearest',rasterized=True)
            ax.axhspan(-90,-89,color='#d9d9d9',lw=0)
            ax.set(xlim=(-180,180),ylim=(-90,90),xticks=[-180,-90,0,90,180],yticks=[-90,-45,0,45,90],aspect='equal')
            ax.set_title(f'({chr(97+row*4+group*2+col)}) {title}, {h} km',loc='left')
            if col==0: ax.set_ylabel('MSO latitude (deg)')
            if row==2: ax.set_xlabel('MSO longitude (deg)')
    cb=sub.colorbar(im,ax=axes,location='bottom',fraction=.05,pad=.035,aspect=45,label=label)
    cb.update_ticks()
for ext in ('png','pdf'):
    fig.savefig(OUT/f'speed_outward_200_400_600km.{ext}',dpi=300,bbox_inches='tight',pad_inches=.12)
plt.close(fig)
meta=json.loads((SOURCE/'metadata.json').read_text(encoding='utf-8'))
meta.update({'source_samples':str(SOURCE/'source.csv'),'layout':'3 rows (200,400,600 km), 4 columns (speed O2+, speed O+, outward O2+, outward O+)',
    'flux':'left: n*norm(V); right: n*Vr only where Vr>0, positive outward flux; cm^-2 s^-1',
    'normalization':'two independent logarithmic scales, each shared across its six panels; limits include all unmasked positive values',
    'color_limits':[[n.vmin,n.vmax] for n in norms],
    'mask':'south-pole row grey; right-hand Vr<=0 cells white',
    'zero':'nonpositive flux masked for log scale; right Vr<=0 excluded'})
for key in ('vmin','vmax','panel_maxima'): meta.pop(key,None)
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
np.savez_compressed(OUT/'plot_data.npz',longitude_deg=lon,latitude_deg=lat,altitude_km=heights,
    speed_flux_cm2_s=speed.data,speed_mask=np.ma.getmaskarray(speed),radial_flux_cm2_s=radial,outward_mask=np.ma.getmaskarray(outward))
print(OUT)
print(meta['color_limits'])
