using MarsTP, StaticArrays, CUDA, Test, LinearAlgebra, Dates, TOML, JLD2, KernelAbstractions

CUDA.functional() || error("CUDA device required")
root = project_path("outputs","actual_mhd_domain_"*Dates.format(now(),"yyyymmdd_HHMMSS")*"_"*string(time_ns()))
mkpath(root)
# Default: the legacy field. Pass the legacy and 100 km field paths to validate both.
paths = isempty(ARGS) ? [data_path("mars_fields_spherical_from_dat.vts")] :
    [MarsTP.resolve_project_path(path) for path in ARGS]
directions = [SA[1.,0.,0.],SA[-1.,0.,0.],SA[0.,1.,0.],SA[0.,-1.,0.],SA[0.,0.,1.],SA[0.,0.,-1.]]
cfg = ForwardTraceConfig(species="O2+",solver=:boris,dt=.005,tspan=(0.,1.))
reports = Dict{String,Any}[]
println("output=",root)
@testset "Actual MHD serial/CUDA radial boundaries" begin
    for (file_id,path) in enumerate(paths)
        fields=load_mhd_fields(path)
        mesh_lo,hi=first(fields.r),last(fields.r)
        lo=max(Rinner,mesh_lo)
        @test all(isfinite,fields.E) && all(isfinite,fields.B)
        file_id==2 && @test mesh_lo == Rm+100e3
        states=SVector{6,Float64}[]
        expected=Symbol[]
        for d in directions
            push!(states,vcat((lo+1500.)*d,-1e5*d));push!(expected,:inner)
            push!(states,vcat((hi-1500.)*d,1e5*d));push!(expected,:outer)
            push!(states,vcat((lo+500e3)*d,1e4*SA[d[2],d[3],d[1]]));push!(expected,:time_limit)
        end
        # On the lower boundary, an inward first step must terminate at t=0.
        push!(states,SA[lo,0.,0.,-1e5,0.,0.]);push!(expected,:inner)
        push!(states,SA[Rm+99e3,0.,0.,-1e5,0.,0.]);push!(expected,:initial_outside)
        for (cadence_id,saveat) in enumerate(((),.02))
            serial=only.(getproperty.(trace_forward(states;config=cfg,fields,saveat),:u))
            gpu=only.(getproperty.(trace_forward(states;config=cfg,fields,saveat,backend=CUDA.CUDABackend()),:u))
            rows=String["particle,expected,serial_retcode,gpu_retcode,serial_end_s,gpu_end_s,time_error_s,max_position_error_m,max_velocity_error_m_s,finite"]
            maxx=maxv=maxt=0.
            for (i,(a,b)) in enumerate(zip(serial,gpu))
                finite=all(u->all(isfinite,u),a.u)&&all(u->all(isfinite,u),b.u)
                @test finite
                @test a.retcode==b.retcode
                @test a.retcode==(expected[i]==:time_limit ? MarsTP.TP.ReturnCode.Success : MarsTP.TP.ReturnCode.Terminated)
                @test length(a.t)==length(b.t)
                @test a.t≈b.t atol=1e-10 rtol=0
                @test all(>(0),diff(a.t)) && all(>(0),diff(b.t))
                dx=maximum(norm(x[1:3]-y[1:3]) for (x,y) in zip(a.u,b.u))
                dv=maximum(norm(x[4:6]-y[4:6]) for (x,y) in zip(a.u,b.u))
                dt=abs(last(a.t)-last(b.t))
                @test dx<.001
                @test dv<.001
                @test dt<1e-10
                @test expected[i]==:time_limit ? last(a.t)≈1. : last(a.t)<1.
                maxx=max(maxx,dx);maxv=max(maxv,dv);maxt=max(maxt,dt)
                push!(rows,join((i,expected[i],a.retcode,b.retcode,last(a.t),last(b.t),dt,dx,dv,finite),','))
            end
            @test last(serial[end].t)==last(gpu[end].t)==0.
            stem="field$(file_id)_cadence$(cadence_id)"
            write(joinpath(root,stem*".csv"),join(rows,"\n")*"\n")
            # ODESolution also references the full field through prob.p.
            # Save only trajectories and statuses, not another copy of the mesh.
            histories(sols)=[(;t=s.t,u=s.u,retcode=string(s.retcode)) for s in sols]
            JLD2.jldsave(joinpath(root,stem*".jld2");states,expected,
                serial=histories(serial),gpu=histories(gpu))
            push!(reports,Dict("file"=>path,"bytes"=>filesize(path),"size"=>collect(size(fields.E)),
                "minimum_grid_altitude_km"=>(mesh_lo-Rm)/1e3,"absorption_altitude_km"=>200.,"maximum_radius_m"=>hi,"particles"=>length(states),
                "saveat"=>string(saveat),"max_position_error_m"=>maxx,"max_velocity_error_m_s"=>maxv,"max_time_error_s"=>maxt))
            println(last(reports));flush(stdout)
        end
        cpu_bounded=[only(trace_forward_bounded([state];config=cfg,fields,backend=CPU(),save_every=1)) for state in states]
        gpu_bounded=trace_forward_bounded(states;config=cfg,fields,backend=CUDA.CUDABackend(),save_every=1)
        bounded_rows=String["particle,status,cpu_end_s,gpu_end_s,end_altitude_km,max_position_error_m,max_velocity_error_m_s,finite"]
        bx=bv=bt=0.
        for (i,(a,b)) in enumerate(zip(cpu_bounded,gpu_bounded))
            status=expected[i]==:initial_outside ? :inner : expected[i]
            @test a.status==b.status==status
            @test a.t≈b.t atol=1e-10 rtol=0
            @test all(u->all(isfinite,u),a.u) && all(u->all(isfinite,u),b.u)
            dx=maximum(norm(x[1:3]-y[1:3]) for (x,y) in zip(a.u,b.u))
            dv=maximum(norm(x[4:6]-y[4:6]) for (x,y) in zip(a.u,b.u))
            dt=abs(last(a.t)-last(b.t))
            @test dx<.001 && dv<.001 && dt<1e-10
            altitude=(norm(last(a.u)[1:3])-Rm)/1e3
            expected[i]==:inner && @test altitude≈200. atol=1e-8
            expected[i]==:initial_outside && @test a.t==b.t==[0.]
            bx=max(bx,dx);bv=max(bv,dv);bt=max(bt,dt)
            push!(bounded_rows,join((i,a.status,last(a.t),last(b.t),altitude,dx,dv,true),','))
        end
        write(joinpath(root,"field$(file_id)_bounded.csv"),join(bounded_rows,"\n")*"\n")
        JLD2.jldsave(joinpath(root,"field$(file_id)_bounded.jld2");states,cpu_bounded,gpu_bounded)
        push!(reports,Dict("file"=>path,"mode"=>"exact spherical CPU/CUDA stopping","particles"=>length(states),
            "absorption_altitude_km"=>200.,"max_position_error_m"=>bx,"max_velocity_error_m_s"=>bv,"max_time_error_s"=>bt))
        # Outside-grid initial conditions stop before querying fields.
        for altitude in (99.,50.,0.,-1.),backend in (nothing,CUDA.CUDABackend())
            start=SA[Rm+altitude*1e3,0.,0.,-1e5,0.,0.]
            sol=only(only(trace_forward([start];config=cfg,fields,backend)).u)
            @test sol.retcode == MarsTP.TP.ReturnCode.Terminated
            @test sol.t == [0.]
            @test sol.u == [start]
            @test all(isfinite,only(sol.u))
        end
    end
end
report=Dict("julia"=>string(VERSION),"TestParticle"=>string(pkgversion(MarsTP.TP)),
    "gpu"=>CUDA.name(CUDA.device()),"mass_species"=>cfg.species,"dt_s"=>cfg.dt,"tspan_s"=>collect(cfg.tspan),
    "radius_m"=>Rm,"field"=>"static total electric and magnetic fields; Cartesian simulation xyz, no new MSO/MSE transform",
    "runs"=>reports,"below_100km_initial"=>"Terminated at t0 before field query, unchanged initial state, no extrapolation", "status"=>"passed")
open(joinpath(root,"summary.toml"),"w") do io;TOML.print(io,report);end
println("PASSED ",root)
