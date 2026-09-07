# Preliminary Vy-integrated detector VDF

The detector is at (0,0,2 Rm), with a traversable 400 km source surface and a 200 km absorbing boundary. The initial calculation uses 41×41×21=35,301 trajectories and rectangular integration over Vy with dVy=20 km/s. The output f_xz has units s² m⁻⁵.

Vx and Vz span −200 to 200 km/s at 10 km/s spacing; Vy spans the same range. The Boris step is −0.1 s and the maximum lookback time is 500 s. Of these trajectories, 426 reached the inner boundary and 34,875 reached the outer boundary, without nonfinite states or time limits.

## Unconverged numerical checks

Five high-value points were compared using Vy spacings of 20, 10, 5 and 2.5 km/s and time steps of 0.1 and 0.05 s. Results vary strongly with Vy spacing, and some are sensitive to the time step. Peaks, integrated density and shape must not be treated as converged scientific results.

At (Vx,Vz)=(0,0) km/s and dt=0.1 s, decreasing Vy spacing from 20 to 2.5 km/s changes f_xz from 5.877e-10 to 7.346e-11 s² m⁻⁵. At (10,150) km/s, the 2.5 km/s grid captures a contribution missed by coarser grids. Narrow Vy peaks need to be resolved before refining Vx, Vz and time. Sampled contributions at |Vy| from 200 to 300 km/s were numerically zero at these five points only; this is not a full tail-convergence check.

The figure uses turbo with logarithmic normalization, without smoothing or gap filling. Light gray denotes numerical zero, potentially including Maxwellian-tail underflow; positive values below the color range are purple. Numerical data are in vdf.csv/JLD2 and refinement results in qa_refinement.csv.
