using MarsTP, StaticArrays, LinearAlgebra, Random
out=ARGS[1];mkpath(out)
n=length(ARGS)>1 ? parse(Int,ARGS[2]) : 5000
fields=load_mhd_fields();param=MarsTP.mhd_param(fields;species="O2+")
itp=build_field_work_interpolators()
cfg=BacktraceConfig(include_ionosphere=false,dt=-.05,tspan=(0.,-500.))
rng=Xoshiro(20260905)
open(joinpath(out,"segments.csv"),"w") do io
 open(joinpath(out,"particles.csv"),"w") do summary
    println(io,"id,x0_Rm,z0_Rm,x1_Rm,z1_Rm,conv_eV_s,hall_eV_s,total_eV_s")
    println(summary,"id,vx_kms,vy_kms,vz_kms,status,end_time_s,total_eV,conv_eV,hall_eV,deltaK_eV,residual_eV")
    for i in 1:n
        speed=(10+190rand(rng))*1000;angle=2pi*rand(rng)
        v0=SA[speed*cos(angle),0.,speed*sin(angle)]
        saved=Ref(SA[0.,0.,2Rm]);saved_t=Ref(0.);last=Ref(saved[]);last_t=Ref(0.)
        bucket=zeros(3)
        function save_segment()
            elapsed=saved_t[]-last_t[]
            elapsed>0 || return
            power=bucket/elapsed
            println(io,join((i,saved[][1]/Rm,saved[][3]/Rm,last[][1]/Rm,last[][3]/Rm,power[2],power[3],power[1]),','))
            saved[]=last[];saved_t[]=last_t[];fill!(bucket,0.)
        end
        function observe(a,b,ta,tb,dw)
            bucket .+= dw;last[]=b;last_t[]=tb
            saved_t[]-tb >= .5-1e-8 && save_segment()
        end
        result=MarsTP._trace_sources(SA[0.,0.,2Rm],v0,param,cfg,(p,v)->0.,nothing;
            outer=fields.r[end],work_itp=itp,segment_observer=observe)
        save_segment()
        result.status != 3 || error("Nonfinite trajectory $i")
        println(summary,join((i,(v0/1000)...,result.status,result.time,result.work_eV...,result.delta_kinetic_eV,result.work_eV[1]-result.delta_kinetic_eV),','))
        i%500==0 && println("Completed $i/$n")
    end
 end
end
