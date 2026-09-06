using MarsTP, StaticArrays, LinearAlgebra, Random
out=ARGS[1];mkpath(out)
n=length(ARGS)>1 ? parse(Int,ARGS[2]) : 5000
mode=length(ARGS)>2 ? ARGS[3] : "original"
mode in ("original","tail_isotropic") || error("Unknown velocity sampling mode")
detector=mode=="original" ? SA[0.,0.,2Rm] : SA[-2Rm,0.,0.]
dt=length(ARGS)>3 ? parse(Float64,ARGS[4]) : -.05
dt<0 || error("dt must be negative")
open(joinpath(out,"case.json"),"w") do io
    print(io,"{\"mode\":\"",mode,"\",\"detector_Rm\":[",join(detector/Rm,','),"],\"dt_s\":",dt,",\"seed\":20260905,\"particles\":",n,"}")
end
fields=load_mhd_fields();param=MarsTP.mhd_param(fields;species="O2+")
itp=build_field_work_interpolators()
cfg=BacktraceConfig(include_ionosphere=false,dt=dt,tspan=(0.,-500.))
rng=Xoshiro(20260905)
open(joinpath(out,"segments.csv"),"w") do io
 open(joinpath(out,"particles.csv"),"w") do summary
    println(io,"id,x0_Rm,z0_Rm,x1_Rm,z1_Rm,conv_eV_s,hall_eV_s,total_eV_s,elapsed_s,lookback_conv_eV,lookback_hall_eV,lookback_total_eV")
    println(summary,"id,vx_kms,vy_kms,vz_kms,status,end_time_s,total_eV,conv_eV,hall_eV,deltaK_eV,residual_eV")
    for i in 1:n
        v0 = if mode=="original"
            speed=(10+190rand(rng))*1000;angle=2pi*rand(rng)
            SA[speed*cos(angle),0.,speed*sin(angle)]
        else
            speed=rand(rng)*1e5;mu=2rand(rng)-1;phi=2pi*rand(rng)
            speed*SA[sqrt(1-mu^2)*cos(phi),sqrt(1-mu^2)*sin(phi),mu]
        end
        saved=Ref(detector);saved_t=Ref(0.);last=Ref(saved[]);last_t=Ref(0.)
        bucket=zeros(3)
        cumulative=zeros(3)
        function save_segment()
            elapsed=saved_t[]-last_t[]
            elapsed>0 || return
            power=bucket/elapsed
            cumulative .+= bucket
            println(io,join((i,saved[][1]/Rm,saved[][3]/Rm,last[][1]/Rm,last[][3]/Rm,power[2],power[3],power[1],elapsed,cumulative[2],cumulative[3],cumulative[1]),','))
            saved[]=last[];saved_t[]=last_t[];fill!(bucket,0.)
        end
        function observe(a,b,ta,tb,dw)
            bucket .+= dw;last[]=b;last_t[]=tb
            saved_t[]-tb >= .5-1e-8 && save_segment()
        end
        result=MarsTP._trace_sources(detector,v0,param,cfg,(p,v)->0.,nothing;
            outer=fields.r[end],work_itp=itp,segment_observer=observe)
        save_segment()
        result.status != 3 || error("Nonfinite trajectory $i")
        println(summary,join((i,(v0/1000)...,result.status,result.time,result.work_eV...,result.delta_kinetic_eV,result.work_eV[1]-result.delta_kinetic_eV),','))
        i%500==0 && println("Completed $i/$n")
    end
 end
end
