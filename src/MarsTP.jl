module MarsTP

using LinearAlgebra
using StaticArrays
using MAT
using JLD2
using ReadVTK
using FastInterpolations
using CairoMakie
using TestParticle
using VelocityDistributionFunctions

import TestParticle as TP
import VelocityDistributionFunctions as VDF

export Rm, Rinner, Router
export project_path, data_path, resolve_project_path
export MHDFields, SpeciesMoments, load_mhd_fields, load_mhd_moments
export load_gitm, load_amps, neutral_properties, hot_oxygen_density
export O2plusSourceRates, load_o2plus_source_rates, assert_same_grid
export BacktraceConfig, run_backtrace_vdf, velocity_axes
export ForwardTraceConfig, trace_forward
export build_field_work_interpolators, field_work
export Electric_field_work_profile, field_work_profile, particle_field_work
export plot_vx_vz, load_detector_vdf
export o2plus_production_density

include("constants.jl")
include("paths.jl")
include("data/mhd.jl")
include("data/atmosphere.jl")
include("data/source_rates.jl")
include("chemistry/reactions.jl")
include("tracing/back_tracing.jl")
include("tracing/forward_tracing.jl")
include("analysis/electric_field_work.jl")
include("analysis/detector_vdf.jl")

end
