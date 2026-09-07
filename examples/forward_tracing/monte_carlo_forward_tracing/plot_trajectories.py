"""Plot a reproducible uniform sample of 5000 propagated trajectories (PNG only)."""
from pathlib import Path
import argparse
import csv
import json
import tomllib
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, Rectangle
from matplotlib.lines import Line2D
import sys
_VIS_ROOT = next(p for p in Path(__file__).resolve().parents if (p / "src" / "MarsTP.jl").exists())
sys.path.insert(0, str(_VIS_ROOT / "src"))
from visualization.mars import bs_mpb, plot_mars


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--output', type=Path, default=Path('trajectories_5000.png'))
    parser.add_argument('--count', type=int, default=5000)
    parser.add_argument('--seed', type=int, default=20260906)
    parser.add_argument('--max-points', type=int, default=300)
    args = parser.parse_args()
    if args.count < 1 or args.max_points < 2:
        raise ValueError('count must be positive and max-points >= 2')
    meta = tomllib.loads((args.run/'metadata.toml').read_text())
    # Reservoir sampling: equal selection probability, no preference for probe hits.
    rng = np.random.default_rng(args.seed)
    selected = []
    eligible = 0
    with (args.run/'particles.csv').open(newline='') as stream:
        for row in csv.DictReader(stream):
            if float(row['rate_weight_s1']) <= 0 or float(row['end_time_s']) <= 0:
                continue
            eligible += 1
            item = (int(row['particle_id']), row['status'])
            if len(selected) < args.count:
                selected.append(item)
            else:
                j = int(rng.integers(eligible))
                if j < args.count:
                    selected[j] = item
    if len(selected) != args.count:
        raise ValueError(f'Only {eligible} eligible trajectories')
    print(f'Selected {len(selected)} of {eligible} eligible particles; reading saved paths', flush=True)
    # Discover group locations without assuming the run batch size.
    wanted = dict(selected)
    paths, statuses, found = [], [], []
    for batch_index, file in enumerate(sorted(args.run.glob('trajectories_*.jld2')), 1):
        with h5py.File(file, 'r') as data:
            for key in data.keys():
                if not key.startswith('p') or not key[1:].isdigit():
                    continue
                pid = int(key[1:])
                if pid not in wanted:
                    continue
                state = data[f'{key}/state']
                assert state.shape[1] == 7 and state.shape[0] > 1
                # Display decimation only; retain both endpoints and raw files.
                ix = np.unique(np.linspace(0, state.shape[0]-1,
                                           min(args.max_points, state.shape[0]), dtype=int))
                xyz = state[ix, 1:4]/meta['Rm_m']
                assert np.isfinite(xyz).all()
                paths.append(xyz)
                statuses.append(wanted[pid])
                found.append(pid)
        if batch_index % 100 == 0:
            print(f'Read {batch_index} batches, {len(found)} selected paths', flush=True)
    assert len(found) == args.count and len(set(found)) == args.count
    mpl.rcParams.update({'font.family': 'Arial', 'font.size': 11})
    colors = dict(zip(('inner', 'outer', 'time_limit'), plt.get_cmap('turbo')([.13, .53, .88])))
    assert set(statuses) <= colors.keys(), set(statuses)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.3), layout='constrained')
    center = meta['detector_Rm']
    side = meta['cube_side_Rm']
    for ax, (i,j), name in zip(axes, ((0,2),(0,1),(1,2)), ('XZ','XY','YZ')):
        ax.add_collection(LineCollection([p[:,[i,j]] for p in paths],
                          colors=[colors[s] for s in statuses], linewidths=.22, alpha=.16,
                          rasterized=True, zorder=1))
        if i == 0:
            bs_mpb(ax=ax, draw_bs=True, draw_mpb=False, sphere=False,
                   boundary_color='black', boundary_ls='--', boundary_lw=1, mars_lw=0, mars_ls='-')
            bs_mpb(ax=ax, draw_bs=False, draw_mpb=True, sphere=False,
                   boundary_color='black', boundary_ls=':', boundary_lw=1, mars_lw=0, mars_ls='-')
        plot_mars(ax=ax, texture=False, facecolor='#b9a296', edgecolor='black', lw=.7, zorder=5)
        ax.add_patch(Circle((0,0), 1+meta['source_altitude_km']*1000/meta['Rm_m'],
                            fill=False, ec='gray', lw=.8, zorder=6))
        ax.add_patch(Rectangle((center[i]-side/2,center[j]-side/2),side,side,
                               fill=False, ec='magenta', lw=1.4, zorder=7))
        ax.set(xlim=(-4.2,4.2), ylim=(-4.2,4.2), aspect='equal',
               xlabel=f'{"XYZ"[i]} ($R_M$)', ylabel=f'{"XYZ"[j]} ($R_M$)')
        ax.set_title(name, loc='left')
    handles = [Line2D([],[],color=c,lw=2,label=s) for s,c in colors.items()]
    handles += [Line2D([],[],color='k',ls='--',label='BS'),
                Line2D([],[],color='k',ls=':',label='MPB'),
                Line2D([],[],color='magenta',label='Probe')]
    axes[1].legend(handles=handles, loc='lower left', fontsize=8, ncol=2,
                   frameon=True, facecolor='white', framealpha=1, edgecolor='none')
    fig.suptitle(f'O$_2^+$: {args.count:,} sampled trajectories, $R_M$ = {meta["Rm_m"]/1000:g} km')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches='tight')
    plt.close(fig)
    args.output.with_suffix('.selection.json').write_text(json.dumps({
        'seed':args.seed, 'eligible_particles':eligible, 'count':len(found),
        'max_display_points':args.max_points, 'particle_ids':sorted(found),
        'coordinate_system':meta['coordinate_system'],
        'boundary_note':'BS/MPB X-rho curves shown only on XY/XZ, assuming +X sunward; contextual, not fitted MHD boundaries.'}, indent=2))
    print(f'Saved {args.output}: {len(found)} trajectories from {eligible} eligible particles')


if __name__ == '__main__':
    main()
