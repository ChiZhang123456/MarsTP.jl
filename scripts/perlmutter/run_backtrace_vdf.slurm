#!/bin/bash
#SBATCH --job-name=marstp_backtrace
#SBATCH --account=mXXXX
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=logs/backtrace_%j.out
#SBATCH --error=logs/backtrace_%j.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs output

export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-128}"

julia --project=. -e '
using MarsTP, StaticArrays

cfg = BacktraceConfig(
    detector_Rm = SA[
        parse(Float64, get(ENV, "DETECTOR_X_RM", "-1.5")),
        parse(Float64, get(ENV, "DETECTOR_Y_RM", "0.0")),
        parse(Float64, get(ENV, "DETECTOR_Z_RM", "1.0")),
    ],
    vx_min_kms = parse(Float64, get(ENV, "VX_MIN_KMS", "-500")),
    vx_max_kms = parse(Float64, get(ENV, "VX_MAX_KMS", "500")),
    vy_min_kms = parse(Float64, get(ENV, "VY_MIN_KMS", "-500")),
    vy_max_kms = parse(Float64, get(ENV, "VY_MAX_KMS", "500")),
    vz_min_kms = parse(Float64, get(ENV, "VZ_MIN_KMS", "-500")),
    vz_max_kms = parse(Float64, get(ENV, "VZ_MAX_KMS", "500")),
    dv_kms = parse(Float64, get(ENV, "DV_KMS", "1")),
    dvy_kms = parse(Float64, get(ENV, "DVY_KMS", "1")),
    dt = -abs(parse(Float64, get(ENV, "DT_SECONDS", "0.2"))),
    stream_vy = true,
    output_file = get(ENV, "OUTPUT_FILE", MarsTP.project_path("output", "backtrace_vdf.jld2")),
    checkpoint_file = get(ENV, "CHECKPOINT_FILE", MarsTP.project_path("output", "backtrace_checkpoint.jld2")),
)

run_backtrace_vdf(cfg)
'
