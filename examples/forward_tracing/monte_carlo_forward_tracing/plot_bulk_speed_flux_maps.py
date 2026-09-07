"""Reproduce the six bulk-speed flux maps from bundled data or fresh sampler CSV."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont

HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv',type=Path,help='CSV from sample_bulk_speed_flux_maps.jl')
    parser.add_argument('--output',type=Path,default=HERE/'bulk_speed_flux_200_400_600km.png')
    args=parser.parse_args()
    if args.csv:
        a=np.genfromtxt(args.csv,delimiter=',',names=True,dtype=None,encoding='utf-8')
        lon=np.unique(a['longitude_deg']);lat=np.unique(a['latitude_deg']);h=np.unique(a['altitude_km'])
        assert list(h)==[200,400,600]
        values=np.empty((3,2,len(lat),len(lon)))
        for row,height in enumerate(h):
            for col,species in enumerate(('O2+','O+')):
                p=a[(a['altitude_km']==height)&(a['species']==species)]
                assert len(p)==len(lat)*len(lon)
                p=np.sort(p,order=['latitude_deg','longitude_deg'])
                v=p['n_m3']*np.sqrt(p['ux_ms']**2+p['uy_ms']**2+p['uz_ms']**2)/1e4
                assert np.isfinite(v).all() and (v>=0).all()
                values[row,col]=v.reshape(len(lat),len(lon))
        mask=np.broadcast_to((lat==-90)[None,None,:,None],values.shape)|(values<=0)
    else:
        with np.load(HERE/'bulk_speed_flux_maps.npz') as a:
            lon=a['longitude_deg'];lat=a['latitude_deg'];h=a['altitude_km']
            values=a['speed_flux_cm2_s'];mask=a['speed_mask']|(values<=0)
    v=np.ma.array(values,mask=mask)
    assert np.isfinite(v.compressed()).all() and v.min()>0
    norm=LogNorm(float(v.min()),float(v.max()))
    findfont('Arial',fallback_to_default=False)
    mpl.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.linewidth':.7,'pdf.fonttype':42})
    fig,axes=plt.subplots(3,2,figsize=(9,9),layout='constrained',sharex=True,sharey=True)
    fig.suptitle(r'Flux magnitude, $n|\mathbf{V}|$',fontsize=14)
    for row,height in enumerate(h):
        for col,species in enumerate((r'O$_2^+$',r'O$^+$')):
            ax=axes[row,col]
            im=ax.pcolormesh(lon,lat,v[row,col],norm=norm,cmap=plt.get_cmap('turbo').with_extremes(bad='white'),shading='nearest',rasterized=True)
            ax.axhspan(-90,-89,color='#d9d9d9',lw=0)
            ax.set(xlim=(-180,180),ylim=(-90,90),xticks=[-180,-90,0,90,180],yticks=[-90,-45,0,45,90],aspect='equal')
            ax.set_title(f'({chr(97+row*2+col)}) {species}, {int(height)} km',loc='left')
            if col==0: ax.set_ylabel('MSO latitude (deg)')
            if row==2: ax.set_xlabel('MSO longitude (deg)')
    fig.colorbar(im,ax=axes,location='bottom',fraction=.05,pad=.035,aspect=45,label=r'$n|\mathbf{V}|$ (cm$^{-2}$ s$^{-1}$)')
    fig.savefig(args.output,dpi=300,bbox_inches='tight',pad_inches=.12)
    plt.close(fig)
    print(f'{args.output}: shared log limits {norm.vmin:g}, {norm.vmax:g} cm^-2 s^-1')

if __name__=='__main__':
    main()
