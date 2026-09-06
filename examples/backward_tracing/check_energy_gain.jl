using MarsTP, StaticArrays, LinearAlgebra
out=ARGS[1];mkpath(out)
fields=load_mhd_fields();param=MarsTP.mhd_param(fields;species="O2+")
itp=build_field_work_interpolators()
open(joinpath(out,"energy_step_check.csv"),"w") do io
    println(io,"vx_kms,vy_kms,vz_kms,dt_s,total_eV,conv_eV,hall_eV,deltaK_eV,residual_eV,status")
    for v in (SA[0.,0.,0.],SA[0.,0.,20.],SA[-100.,0.,100.],SA[200.,200.,200.],SA[10.,-20.,-20.]), dt in (-.1,-.05,-.025)
        cfg=BacktraceConfig(include_ionosphere=false,dt=dt)
        r=MarsTP._trace_sources(SA[0.,0.,2Rm],v*1000,param,cfg,(p,v)->0.,nothing;work_itp=itp,outer=last(fields.r))
        println(io,join((v...,dt,r.work_eV...,r.delta_kinetic_eV,r.work_eV[1]-r.delta_kinetic_eV,r.status),','))
        flush(io)
    end
end
