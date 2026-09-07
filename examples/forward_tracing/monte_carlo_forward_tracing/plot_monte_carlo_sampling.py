"""Three-panel 2D O2+ importance-sampling illustration, PNG only.

This is a normalized two-dimensional teaching model with vz=0, not a slice
of a 3D VDF. The production MarsTP source sampler remains three-dimensional.
Colors in panels b/c are per-particle weights, not histogram sums or PSD.
"""
from pathlib import Path
import argparse
import math
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.font_manager import findfont


def sample_demo(count=100_000, seed=20260907):
    if count < 2:
        raise ValueError('At least two proposal draws are required')
    n, mass, qe, temperature_ev = 5e6, 5.352390155808e-26, 1.602176634e-19, 10.
    bulk = np.array([-10e3, 0.])
    factor, area = 4., 1.
    sigma = math.sqrt(qe*temperature_ev/mass)
    rng = np.random.default_rng(seed)
    velocity = bulk+sigma*math.sqrt(factor)*rng.standard_normal((count,2))
    d2 = np.sum((velocity-bulk)**2,axis=1)
    # d=2: g2/gs2 = factor * exp[-d2/(2 sigma^2)*(1-1/factor)].
    importance = factor*np.exp(-d2/(2*sigma**2)*(1-1/factor))
    density_weight = n*importance/importance.sum()
    # Normal is (-1,0,0); inward draws retained with zero rate, no rejection.
    rate_weight = n*area*np.maximum(-velocity[:,0],0)*importance/count
    un = -bulk[0]
    phi = math.exp(-.5*(un/sigma)**2)/math.sqrt(2*math.pi)
    cdf = .5*math.erfc(-un/sigma/math.sqrt(2))
    analytic_rate = n*area*(sigma*phi+un*cdf)
    rate_se = math.sqrt(count)*np.std(rate_weight,ddof=1)
    neff = importance.sum()**2/np.sum(importance**2)
    weighted_mean = np.sum(density_weight[:,None]*velocity,axis=0)/n
    assert np.isclose(density_weight.sum(),n,rtol=1e-14)
    assert np.all(rate_weight[velocity[:,0]>=0]==0)
    assert np.all(abs(weighted_mean-bulk)<6*sigma/math.sqrt(neff))
    assert abs(rate_weight.sum()-analytic_rate)<6*rate_se
    return dict(n=n,mass=mass,qe=qe,temperature_ev=temperature_ev,bulk=bulk,
        factor=factor,area=area,sigma=sigma,velocity=velocity,importance=importance,
        density_weight=density_weight,rate_weight=rate_weight,analytic_rate=analytic_rate,
        rate_se=rate_se,neff=neff,weighted_mean=weighted_mean)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('monte_carlo_sampling.png'))
    parser.add_argument('--count',type=int,default=100_000)
    parser.add_argument('--seed',type=int,default=20260907)
    args=parser.parse_args()
    d=sample_demo(args.count,args.seed)
    mpl.rcParams.update({'font.family':'Arial','font.size':10,'mathtext.fontset':'dejavusans'})
    findfont('Arial',fallback_to_default=False)
    fig,axes=plt.subplots(1,3,figsize=(15,5.1),layout='constrained')
    # Include all samples; use the same square velocity window in all panels.
    radius=math.ceil(np.max(abs(d['velocity']-d['bulk']))/5000)*5000
    xlim=(d['bulk'][0]-radius,d['bulk'][0]+radius);ylim=(-radius,radius)
    x=np.linspace(*xlim,361);y=np.linspace(*ylim,361)
    xx,yy=np.meshgrid(x,y,indexing='xy')
    f2=d['n']/(2*math.pi*d['sigma']**2)*np.exp(-((xx-d['bulk'][0])**2+yy**2)/(2*d['sigma']**2))
    im=axes[0].pcolormesh(x/1000,y/1000,f2,cmap='turbo',shading='auto',norm=Normalize(0,f2.max()))
    axes[0].set_title('(a) Analytic Maxwellian',loc='left')
    fig.colorbar(im,ax=axes[0],label=r'$f_{xy}$ (s$^2$ m$^{-5}$)',shrink=.79)
    for ax,weights,title,label in zip(axes[1:],(d['density_weight'],d['rate_weight']),
            ('(b) MC density weight','(c) MC flux (rate) weight'),
            (r'$W_{n,i}$ (m$^{-3}$)',r'$Q_i$ (s$^{-1}$)')):
        positive=weights>0
        if np.any(~positive):
            ax.scatter(*d['velocity'][~positive].T/1000,s=1,c='#cccccc',linewidths=0,rasterized=True)
        order=np.flatnonzero(positive)[np.argsort(weights[positive])]
        im=ax.scatter(*d['velocity'][order].T/1000,c=weights[order],s=1,
                      cmap='turbo',norm=Normalize(0,weights.max()),linewidths=0,rasterized=True)
        ax.set_title(title,loc='left')
        fig.colorbar(im,ax=ax,label=label,shrink=.79)
    axes[2].axvline(0,color='#555555',ls=':',lw=.8)
    for ax in axes:
        ax.set(xlim=np.array(xlim)/1000,ylim=np.array(ylim)/1000,
               xlabel=r'$v_x$ (km/s)',ylabel=r'$v_y$ (km/s)',aspect='equal')
    fig.suptitle('O$_2^+$: $n$ = 5 cm$^{-3}$, $U_x$ = -10 km/s, $T$ = 10 eV, $v_z$ = 0\n'
                 f'{args.count:,} proposal draws; '+r'$T_s=4T$; source area = 1 m$^2$, normal = $-\hat{x}$')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=220,bbox_inches='tight')
    plt.close(fig)
    print(f'density sum = {d["density_weight"].sum():.12g} m^-3')
    print(f'rate sum = {d["rate_weight"].sum():.12g} s^-1; analytic = {d["analytic_rate"]:.12g} s^-1')
    print(f'rate discrepancy = {(d["rate_weight"].sum()-d["analytic_rate"])/d["rate_se"]:.4f} MC standard errors')
    print(f'weighted mean velocity = {d["weighted_mean"]/1000} km/s; Neff = {d["neff"]:.1f}')
    print(f'sigma = {d["sigma"]/1000:.6f} km/s; saved {args.output}')


if __name__=='__main__':
    main()
