using MarsTP, StaticArrays, JLD2
out=ARGS[1]
cfg=BacktraceConfig(detector_Rm=SA[0.,0.,2.],ionosphere_altitude_km=400.,
    vx_min_kms=-200.,vx_max_kms=200.,vz_min_kms=-200.,vz_max_kms=200.,
    vy_min_kms=-200.,vy_max_kms=200.,dv_kms=10.,dvy_kms=20.,
    tspan=(0.,-500.),dt=-.1,progress_interval=1000,
    checkpoint_file=joinpath(out,"checkpoint.jld2"),output_file=joinpath(out,"vdf.jld2"))
println("Julia $VERSION; TestParticle $(pkgversion(MarsTP.TP)); threads=$(Threads.nthreads())")
r=run_backtrace_vdf(cfg)
@assert sum(r.status_counts)==41*41*21
@assert r.status_counts[4]==0
@assert all(isfinite,r.f2d_xz) && all(>=(0),r.f2d_xz)
@assert r.f2d_xz ≈ r.f2d_volume+r.f2d_ionosphere
open(joinpath(out,"vdf.csv"),"w") do io
    println(io,"vx_kms,vz_kms,f_total_s2_m5,f_volume_s2_m5,f_sheet_s2_m5,crossings")
    for j in eachindex(r.vz_km),i in eachindex(r.vx_km)
        println(io,join((r.vx_km[i],r.vz_km[j],r.f2d_xz[i,j],r.f2d_volume[i,j],r.f2d_ionosphere[i,j],r.ionosphere_crossings[i,j]),','))
    end
end
println("STATUS ",r.status_counts)
println("RANGE ",extrema(r.f2d_xz))
println("OUTPUT ",out)
