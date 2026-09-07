"""Re-bin finite-volume O2+ Monte Carlo residence records in SI units.

Usage: python analyze_monte_carlo.py RUN_DIR --output-dir NEW_DIR
Defaults: 5 km/s bins on each axis from -500 to 500 km/s.
The panels integrate over the entire omitted velocity axis; no slices.
"""
from pathlib import Path
import argparse
import json
import math
import tomllib
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.font_manager import findfont


def load_csv(path):
    a=np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf8", ndmin=1)
    # Internal compatibility alias only. New files explicitly store rate_weight_s1.
    a.dtype.names=tuple("weight_s1" if n=="rate_weight_s1" else n for n in a.dtype.names)
    return a


def velocity_segments(v0, v1, edges):
    """Exact bin residence fractions for linearly interpolated velocity."""
    cuts = [0., 1.]
    for k in range(3):
        if v1[k] != v0[k]:
            f = (edges[k] - v0[k]) / (v1[k] - v0[k])
            cuts.extend(f[(f > 0) & (f < 1)])
    cuts = np.unique(cuts)
    for a, b in zip(cuts[:-1], cuts[1:]):
        vmid = v0 + .5 * (a+b) * (v1-v0)
        idx = tuple(int(np.searchsorted(edges[k], vmid[k], side="right")-1) for k in range(3))
        if any(i < 0 or i >= len(edges[k])-1 for k, i in enumerate(idx)):
            raise ValueError("Velocity outside requested bins; enlarge --vmax-kms")
        yield idx, b-a


def bin_residence(records, edges, volume):
    shape = tuple(len(e)-1 for e in edges)
    numbers = np.zeros(shape)
    # Each particle is one independent draw; multiple time samples are correlated.
    per_particle_bin = {}
    total_number = 0.
    current = np.zeros(3)
    for r in records:
        dt = r["t1_s"]-r["t0_s"]
        if dt <= 0 or r["weight_s1"] < 0:
            raise ValueError("Invalid dwell time/weight")
        v0 = np.array([r[f"v{k}0_ms"] for k in "xyz"])
        v1 = np.array([r[f"v{k}1_ms"] for k in "xyz"])
        number = r["weight_s1"] * dt
        total_number += number
        current += number*(v0+v1)/2/volume
        for idx, fraction in velocity_segments(v0, v1, edges):
            amount = number*fraction
            numbers[idx] += amount
            key = (int(r["particle_id"]), *idx)
            per_particle_bin[key] = per_particle_bin.get(key, 0.)+amount
    sumsq = np.zeros(shape)
    unique = np.zeros(shape, dtype=np.int32)
    for key, amount in per_particle_bin.items():
        idx = key[1:]
        sumsq[idx] += amount**2
        unique[idx] += 1
    neff = np.divide(numbers**2, sumsq, out=np.zeros(shape), where=sumsq > 0)
    dv3 = np.diff(edges[0])[:,None,None]*np.diff(edges[1])[None,:,None]*np.diff(edges[2])[None,None,:]
    f3d = numbers/volume/dv3
    assert np.isclose(numbers.sum(), total_number, rtol=1e-12, atol=1e-100)
    return f3d, numbers, unique, neff, current


def draw_panels(path, edges, a, b, label, subtitle, hits, neff, meta, zoom=False):
    mpl.rcParams.update({"font.family":"Arial", "font.size":10, "pdf.fonttype":42, "svg.fonttype":"none"})
    findfont("Arial", fallback_to_default=False)
    positive = np.r_[a[a>0], b[b>0]]
    norm = None
    if positive.size:
        vmin, vmax = positive.min(), positive.max()
        if vmax <= vmin:
            vmin, vmax = vmin/2, vmax*2
        norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, axes = plt.subplots(1,2,figsize=(10.4,4.8),layout="constrained")
    cmap=plt.get_cmap("turbo").copy()
    cmap.set_bad("#f0f0f0")
    supports=[]
    for matrix,ydim in zip((a,b),(1,2)):
        ix,iy=np.nonzero(matrix>0)
        supports.append(None if len(ix)==0 else (edges[0][ix.min()],edges[0][ix.max()+1],edges[ydim][iy.min()],edges[ydim][iy.max()+1]))
    available=[s for s in supports if s is not None]
    if zoom and available:
        xlo=min(s[0] for s in available);xhi=max(s[1] for s in available)
        span=max(xhi-xlo,max(s[3]-s[2] for s in available))+4*np.diff(edges[0])[0]
        xmid=(xlo+xhi)/2
    for ax, matrix, ydim, panel in zip(axes,(a,b),(1,2),("a","b")):
        im=ax.pcolormesh(edges[0]/1e3,edges[ydim]/1e3,np.ma.masked_less_equal(matrix.T,0),cmap=cmap,norm=norm,rasterized=True)
        ax.set(xlabel=r"$v_x$ (km/s)",ylabel=rf"$v_{'y' if ydim==1 else 'z'}$ (km/s)",aspect="equal")
        ax.set_title(f"({panel}) "+("Vx-Vy" if ydim==1 else "Vx-Vz"),loc="left")
        ax.axhline(0,lw=.5,c="gray",alpha=.4);ax.axvline(0,lw=.5,c="gray",alpha=.4)
        ax.set_facecolor("#f0f0f0")
        if not np.any(matrix>0):
            ax.text(.5,.5,"No sampled contribution in this view",ha="center",transform=ax.transAxes,fontsize=9)
        elif zoom:
            s=supports[ydim-1];ymid=(s[2]+s[3])/2
            ax.set_xlim((xmid-span/2)/1e3,(xmid+span/2)/1e3)
            ax.set_ylim((ymid-span/2)/1e3,(ymid+span/2)/1e3)
        if meta.get("plot_limit_kms") is not None:
            lim=meta["plot_limit_kms"]
            ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
            mx=(edges[0][:-1]<lim*1e3)&(edges[0][1:]>-lim*1e3)
            my=(edges[ydim][:-1]<lim*1e3)&(edges[ydim][1:]>-lim*1e3)
            if np.any(matrix>0) and not np.any(matrix[np.ix_(mx,my)]>0):
                ax.text(.5,.5,f"Sampled signal lies outside\nthe displayed ±{lim:g} km/s limits",ha="center",va="center",transform=ax.transAxes,fontsize=9)
    if positive.size:
        fig.colorbar(im,ax=axes,label=label,shrink=.78)
    probe=", ".join(f"{v:g}" for v in meta["detector_Rm"])
    fig.suptitle(f"O$_2^+$, probe ({probe}) $R_M$, cube side {meta['cube_side_Rm']:g} $R_M$\n"+subtitle)
    fig.savefig(path.with_suffix(".png"),dpi=190,bbox_inches="tight",pad_inches=.12)
    plt.close(fig)


def main():
    # The production estimator stores the full 200^3 grid sparsely.
    # bin_residence remains available as a small dense reference for tests.
    from analyze_probe import main as sparse_main
    sparse_main()


if __name__ == "__main__":
    main()
