# Run from any directory with julia --project=<MarsTP.jl> <this script>.
using MarsTP
using StaticArrays
using JLD2

# An isolated output directory preserves earlier runs.
mkpath(project_path("output"))
output_dir = mktempdir(project_path("output"); prefix = "smoke_", cleanup = false)
cfg = BacktraceConfig(
    detector_Rm = SA[-1.5, 0.0, 1.0],
    vx_min_kms = -1.0, vx_max_kms = 1.0,
    vy_min_kms = -1.0, vy_max_kms = 1.0,
    vz_min_kms = -1.0, vz_max_kms = 1.0,
    dv_kms = 1.0, dvy_kms = 1.0,
    tspan = (0.0, -2.0), dt = -0.2,
    solver = :boris, stream_vy = true, progress_interval = 9,
    checkpoint_file = joinpath(output_dir, "checkpoint.jld2"),
    output_file = joinpath(output_dir, "vdf.jld2"),
)
result = run_backtrace_vdf(cfg)
@assert size(result.f2d_xz) == (3, 3)
@assert all(isfinite, result.f2d_xz)
@assert all(>=(0), result.f2d_xz)
@assert sum(result.status_counts) == 27
@assert result.status_counts[4] == 0 "Nonfinite trajectory encountered"
@assert load(cfg.output_file, "f2d_xz") == result.f2d_xz
println("SMOKE_OK: 27 trajectories; status_counts = ", result.status_counts)
println("Output: ", cfg.output_file)
println("This checks execution and output integrity, not physical convergence.")
