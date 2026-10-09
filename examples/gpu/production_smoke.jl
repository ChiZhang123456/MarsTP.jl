# Full production path, 60 particles/backend, with source weights, CSV, JLD2,
# termination counts and three work summaries. Every run gets a new directory.
using MarsTP, Test, Dates, TOML, LinearAlgebra
include("../forward_tracing/monte_carlo_forward_tracing/monte_carlo_shell.jl")
root=project_path("outputs","gpu_production_validation_"*Dates.format(now(),"yyyymmdd_HHMMSS")*"_"*string(time_ns()))
tmax=isempty(ARGS) ? .2 : parse(Float64,ARGS[1])
options=(;cell_stride=2000,per_cell=10,tmax,dt=.1,batch_size=32,work_mode=:summary)
for backend in (:cpu,:cuda)
    ShellMonteCarlo.run_monte_carlo(joinpath(root,String(backend)),ShellMonteCarlo.Config(;options...,tracing_backend=backend))
end
records=Dict{Symbol,Dict{Int,Any}}()
for backend in (:cpu,:cuda)
    records[backend]=Dict{Int,Any}()
    foreach_saved_trajectory(r->(records[backend][r.particle_id]=r),joinpath(root,String(backend)))
end
@test length(records[:cpu])==length(records[:cuda])==60
@test keys(records[:cpu])==keys(records[:cuda])
maxwork=0.;maxposition=0.;maxvelocity=0.
@testset "Actual-MHD production CPU/GPU output" begin
    for id in keys(records[:cpu])
        a,b=records[:cpu][id],records[:cuda][id]
        @test a.termination_code==b.termination_code
        @test a.rate_weight_s==b.rate_weight_s
        @test a.source_density_weight_m3==b.source_density_weight_m3
        @test a.trajectory.t ≈ b.trajectory.t atol=1e-9
        for (x,y) in zip(a.trajectory.u,b.trajectory.u)
            global maxposition=max(maxposition,norm(x[1:3]-y[1:3]))
            global maxvelocity=max(maxvelocity,norm(x[4:6]-y[4:6]))
        end
        @test maximum(norm(x[1:3]-y[1:3]) for (x,y) in zip(a.trajectory.u,b.trajectory.u))<.01
        @test maximum(norm(x[4:6]-y[4:6]) for (x,y) in zip(a.trajectory.u,b.trajectory.u))<.01
        for name in keys(a.work.summary)
            @test a.work.summary[name] ≈ b.work.summary[name] atol=1e-7 rtol=1e-8
        end
        for name in ("total_eV","conv_eV","hall_eV")
            global maxwork=max(maxwork,abs(a.work.summary[name]-b.work.summary[name]))
        end
    end
    cpu=TOML.parsefile(joinpath(root,"cpu","completion.toml"))
    gpu=TOML.parsefile(joinpath(root,"cuda","completion.toml"))
    @test cpu["complete"] && gpu["complete"]
    @test cpu["status_counts"]==gpu["status_counts"]
    # The CSV preserves existing columns and adds convection/Hall work.
    function csvrows(path)
        lines=readlines(path);headers=split(first(lines),',')
        return headers,[split(line,',') for line in lines[2:end]]
    end
    h,a=csvrows(joinpath(root,"cpu","particles.csv"));g,b=csvrows(joinpath(root,"cuda","particles.csv"))
    @test h==g
    @test length(a)==length(b)==60
    for name in ("work_eV","work_conv_eV","work_hall_eV","field_sum_residual_eV","probe_residence_s")
        column=findfirst(==(name),h)
        @test column!==nothing
        @test parse.(Float64,getindex.(a,column)) ≈ parse.(Float64,getindex.(b,column)) atol=1e-7 rtol=1e-8
    end
end
println((output=root,max_work_difference_eV=maxwork,max_position_difference_m=maxposition,max_velocity_difference_ms=maxvelocity))
