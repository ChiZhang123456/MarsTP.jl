"""Intersect every saved trajectory with three cubes, without reintegration.

Writes clipped positions/times and the Boris drift velocity. Run
synchronize_reprobe.jl next to evaluate synchronized endpoint velocities.
"""
from pathlib import Path
import argparse
import csv
import json
import tomllib
import shutil
import h5py
import numpy as np
from numba import njit

PROBES = {'probe_1_0_2': [1., 0., 2.],
          'probe_0_0_2': [0., 0., 2.],
          'probe_m1p5_0_1': [-1.5, 0., 1.]}


@njit
def intersections(state, centers, side):
    hits = []
    for j in range(len(state)-1):
        for probe in range(len(centers)):
            near, far, entry, leave = 0., 1., 0, 0
            for k in range(3):
                a = state[j, k+1]
                d = state[j+1, k+1]-a
                lo, hi = centers[probe,k]-side/2, centers[probe,k]+side/2
                if d == 0:
                    if a < lo or a > hi:
                        far = -1.
                        break
                else:
                    q0, q1 = (lo-a)/d, (hi-a)/d
                    f0, f1 = -(k+1), k+1
                    if q0 > q1:
                        q0, q1, f0, f1 = q1, q0, f1, f0
                    if q0 >= near:
                        near, entry = q0, f0
                    if q1 <= far:
                        far, leave = q1, f1
                    if near >= far:
                        break
            if near < far:
                hits.append((probe, j, near, far, entry, leave))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run', type=Path)
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    meta = tomllib.loads((args.run/'metadata.toml').read_text())
    shutil.copyfile(args.run/'metadata.toml',args.output/'source_metadata.toml')
    assert tomllib.loads((args.run/'completion.toml').read_text())['complete']
    weights = np.zeros(meta['n_particles']+1)
    with (args.run/'particles.csv').open(newline='') as f:
        for row in csv.DictReader(f):
            weights[int(row['particle_id'])] = float(row['rate_weight_s1'])
    centers = np.array(list(PROBES.values()))*meta['Rm_m']
    handles, writers = [], []
    header = ['particle_id','rate_weight_s1','t0_s','t1_s',
              'x0_m','y0_m','z0_m','x1_m','y1_m','z1_m',
              'step_t0_s','half_vx_ms','half_vy_ms','half_vz_ms','entry_face','exit_face']
    for name, center in PROBES.items():
        p = args.output/name
        p.mkdir()
        info = dict(meta, detector_Rm=center, source_run_directory=str(args.run.resolve()),
                    extraction='piecewise-linear saved positions; exact slab intersections')
        (p/'metadata.json').write_text(json.dumps(info, indent=2))
        handle = (p/'clipped_segments.csv').open('w', newline='')
        handles.append(handle)
        writer = csv.writer(handle); writer.writerow(header); writers.append(writer)
    groups = 0
    counts = np.zeros(len(PROBES), dtype=int)
    try:
        for batch, file in enumerate(sorted(args.run.glob('trajectories_*.jld2')), 1):
            with h5py.File(file, 'r') as data:
                for key in data.keys():
                    if not key.startswith('p') or not key[1:].isdigit():
                        continue
                    groups += 1
                    pid = int(key[1:])
                    if weights[pid] == 0:
                        continue
                    state = data[f'{key}/state'][()]
                    if len(state) < 2:
                        continue
                    assert state.shape[1] == 7 and np.isfinite(state).all()
                    for probe, j, near, far, entry, leave in intersections(state, centers, meta['cube_side_Rm']*meta['Rm_m']):
                        a, b = state[j], state[j+1]
                        assert b[0] > a[0]
                        pa, pb = a+near*(b-a), a+far*(b-a)
                        half = (b[1:4]-a[1:4])/(b[0]-a[0])
                        writers[probe].writerow([pid, weights[pid], pa[0], pb[0],
                            *pa[1:4], *pb[1:4], a[0], *half, entry, leave])
                        counts[probe] += 1
            if batch % 25 == 0:
                for h in handles: h.flush()
                print(f'{batch} batches, {groups} particle groups, segments={counts.tolist()}', flush=True)
    finally:
        for h in handles: h.close()
    assert groups == meta['n_particles']
    (args.output/'extraction_complete.json').write_text(json.dumps({
        'groups_checked':groups, 'segments':dict(zip(PROBES,counts.tolist()))}, indent=2))
    print('Complete:', counts.tolist())


if __name__ == '__main__':
    main()
