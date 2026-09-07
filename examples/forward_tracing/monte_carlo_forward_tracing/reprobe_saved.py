"""Compatibility CLI delegating to analyze_saved_probes.jl.

The intersections helper is retained only as an independent test reference.
Production calculations use MarsTP.ForwardPSDAccumulator.
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
                    if a < lo or a >= hi:
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
    # Compatibility entry point: delegate to the single Julia production path.
    import subprocess
    parser=argparse.ArgumentParser(description="Analyze saved trajectories using MarsTP.forward_psd's shared core")
    parser.add_argument('run',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    repo=Path(__file__).resolve().parents[3]
    subprocess.run(['julia','--startup-file=no',f'--project={repo}',
                    str(Path(__file__).with_name('analyze_saved_probes.jl')),
                    str(args.run.resolve()),str(args.output.resolve())],check=True)


if __name__=='__main__':
    main()
