using MarsTP, StaticArrays, LinearAlgebra, Test
source = load_ionosphere_source(;altitude_km=400.)
for h in (200.,400.,800.)
    s=IonosphereSource(source.n,source.temperature,source.velocity,Rm+h*1e3)
    for d in (SA[1.,0.,0.],SA[0.,0.,1.],SA[1.,1.,1.])
        println("h=$h dir=$d ",ionosphere_properties(s,d))
    end
end
cfg=BacktraceConfig(detector_Rm=SA[1.12,0.,0.],ionosphere_altitude_km=400.,
    vx_min_kms=10.,vx_max_kms=10.,vy_min_kms=0.,vy_max_kms=0.,vz_min_kms=0.,vz_max_kms=0.,
    tspan=(0.,-2.),dt=-.05,stream_vy=false,progress_interval=1,
    output_file=joinpath(mktempdir(project_path("outputs");cleanup=false),"ionosphere_smoke.jld2"))
r=run_backtrace_vdf(cfg)
@test sum(r.status_counts)==1
@test r.status_counts[3]==1
@test r.f2d_ionosphere[1]>0
@test r.f2d_xz ≈ r.f2d_volume+r.f2d_ionosphere
@test all(isfinite,r.f2d_xz)
println("MHD_SMOKE: ",r.status_counts," f_iono=",r.f2d_ionosphere," output=",cfg.output_file)
