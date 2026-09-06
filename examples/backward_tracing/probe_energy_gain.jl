using MarsTP, StaticArrays, JLD2
out=ARGS[1]
smoke=length(ARGS)>1 && ARGS[2]=="smoke"
cfg=BacktraceConfig(detector_Rm=SA[0.,0.,2.],include_work=true,
    vx_min_kms=-200.,vx_max_kms=200.,vz_min_kms=-200.,vz_max_kms=200.,
    vy_min_kms=-200.,vy_max_kms=200.,dv_kms=smoke ? 200. : 10.,dvy_kms=smoke ? 200. : 20.,
    dt=-.1,tspan=(0.,-500.),progress_interval=1000,
    checkpoint_file=joinpath(out,"checkpoint.jld2"),output_file=joinpath(out,"energy.jld2"))
println("Julia $VERSION; threads=$(Threads.nthreads()); smoke=$smoke")
r=run_backtrace_vdf(cfg)
@assert r.status_counts[4]==0
# Independent slice sum detects inconsistencies in threaded per-trajectory output.
@assert all(isapprox.(dropdims(sum(r.psd3d;dims=2);dims=2)*(cfg.dvy_kms*1e3),r.f2d_xz;rtol=1e-12,atol=0))
open(joinpath(out,"energy_xz.csv"),"w") do io
    println(io,"vx_kms,vz_kms,f_s2_m5,total_eV,conv_eV,hall_eV")
    for iz in eachindex(r.vz_km),ix in eachindex(r.vx_km)
        println(io,join((r.vx_km[ix],r.vz_km[iz],r.f2d_xz[ix,iz],r.mean_work2d_eV[ix,iz,:]...),','))
    end
end
open(joinpath(out,"energy_3d.csv"),"w") do io
    println(io,"vx_kms,vy_kms,vz_kms,psd_s3_m6,status,total_eV,conv_eV,hall_eV,deltaK_eV,residual_eV")
    for iz in eachindex(r.vz_km),iy in eachindex(r.vy_km),ix in eachindex(r.vx_km)
        w=r.work3d_eV[ix,iy,iz,:];dk=r.delta_kinetic3d_eV[ix,iy,iz]
        println(io,join((r.vx_km[ix],r.vy_km[iy],r.vz_km[iz],r.psd3d[ix,iy,iz],r.status3d[ix,iy,iz],w...,dk,w[1]-dk),','))
    end
end
println("STATUS ",r.status_counts)
println("MAX ENERGY RESIDUAL eV ",maximum(abs,r.work3d_eV[:,:,:,1]-r.delta_kinetic3d_eV))

open(joinpath(out,"refine_points.csv"),"w") do io
    println(io,"vx_kms,vy_kms,vz_kms")
    for iz in eachindex(r.vz_km),iy in eachindex(r.vy_km),ix in eachindex(r.vx_km)
        w=r.work3d_eV[ix,iy,iz,1]; dk=r.delta_kinetic3d_eV[ix,iy,iz]
        abs(w-dk)>0.01max(abs(w),abs(dk),1.) || continue
        println(io,join((r.vx_km[ix],r.vy_km[iy],r.vz_km[iz]),','))
    end
end
