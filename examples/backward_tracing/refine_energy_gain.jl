using MarsTP, StaticArrays, DelimitedFiles
out=ARGS[1]
detector=length(ARGS)>1 && ARGS[2]=="tail" ? SA[-2Rm,0.,0.] : SA[0.,0.,2Rm]
points=readdlm(joinpath(out,"refine_points.csv"),',',Float64;skipstart=1)
fields=load_mhd_fields();param=MarsTP.mhd_param(fields;species="O2+")
itp=build_field_work_interpolators()
open(joinpath(out,"energy_refinement.csv"),"w") do io
    println(io,"vx_kms,vy_kms,vz_kms,dt_s,total_eV,conv_eV,hall_eV,deltaK_eV,residual_eV,status")
    for row in eachrow(points), dt in (-.05,-.025,-.0125,-.00625)
        v=SVector{3,Float64}(row)*1000
        cfg=BacktraceConfig(include_ionosphere=false,dt=dt)
        r=MarsTP._trace_sources(detector,v,param,cfg,(p,v)->0.,nothing;work_itp=itp,outer=last(fields.r))
        println(io,join((row...,dt,r.work_eV...,r.delta_kinetic_eV,r.work_eV[1]-r.delta_kinetic_eV,r.status),','))
        flush(io)
    end
end
