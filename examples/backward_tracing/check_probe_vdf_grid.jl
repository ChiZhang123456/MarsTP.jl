using MarsTP, StaticArrays, LinearAlgebra, DelimitedFiles
out=ARGS[1]
points=readdlm(joinpath(out,"qa_points.csv"),',',Float64)
fields=load_mhd_fields();sr=load_o2plus_source_rates(;fill_invalid=false);gitm=load_gitm()
shell=load_ionosphere_source(;altitude_km=400.)
param=MarsTP.mhd_param(fields;species="O2+")
mass=MarsTP.TP.SpeciesDict["O2+"].m
q=MarsTP.TP.build_interpolator(MarsTP.TP.StructuredGrid,sr.production_density,fields.r,fields.theta,fields.phi)
volume(p,v)=MarsTP._volume_source(q,gitm,mass,p,v)
ion(p,v)=ionosphere_distribution(shell,p,v;mass)
open(joinpath(out,"qa_refinement.csv"),"w") do io
    println(io,"vx_kms,vz_kms,dt_s,vy_kms,f3d")
    for row in eachrow(points),dt in (-.1,-.05)
        cfg=BacktraceConfig(dt=dt,tspan=(0.,-500.))
        vyy=collect(-300.:2.5:300.)
        f=zeros(length(vyy))
        Threads.@threads for i in eachindex(vyy)
            r=MarsTP._trace_sources(SA[0.,0.,2Rm],SA[row[1],vyy[i],row[2]]*1e3,param,cfg,volume,ion;outer=last(fields.r))
            @assert r.status!=3
            f[i]=r.volume+r.ionosphere
        end
        for i in eachindex(vyy)
            println(io,join((row[1],row[2],dt,vyy[i],f[i]),','))
        end
        flush(io)
    end
end
println("REFINEMENT_OK")
