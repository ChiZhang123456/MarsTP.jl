using JLD2, TOML
@testset "Shared streaming PSD and batch IO" begin
    # Keep diagnostic files; this test never recursively deletes directories.
    folder=mktempdir(;cleanup=false)
    a=(;t=[0.,1.,2.],u=[SA[-2.,0.,0.,-4.,0.,0.],SA[0.,0.,0.,0.,0.,0.],SA[2.,0.,0.,4.,0.,0.]])
    b=(;t=[0.,2.],u=[SA[-2.,1.,0.,1.,0.,0.],SA[2.,1.,0.,1.,0.,0.]])
    c=(;t=[0.],u=[a.u[1]])
    kw=(;detector_m=zeros(3),side_m=2.,vlim=2.,vgrid=4)
    mem=forward_psd([a,b,c];kw...,rate_weights_s=[16.,8.,0.])
    for (i,traj) in enumerate((a,b,c))
        path=joinpath(folder,"trajectories_$(lpad(i,5,'0')).jld2")
        write_trajectory_batch(path,[traj];particle_ids=[i],rate_weights_s=[[16.,8.,0.][i]],
            source_density_weights_m3=[1.],cell_ids=[i],termination_codes=["time_limit"])
    end
    open(joinpath(folder,"metadata.toml"),"w") do io
        TOML.print(io,Dict("n_particles"=>3,"species"=>"O2+"))
    end
    open(joinpath(folder,"completion.toml"),"w") do io;TOML.print(io,Dict("complete"=>true));end
    disk=forward_psd_saved(folder;kw...)
    @test disk.psd == mem.psd
    @test disk.residence_s == mem.residence_s
    @test disk.outside_vlim_residence_s == mem.outside_vlim_residence_s
    @test disk.particle_ids == [1,2,3]
    @test disk.retcodes == fill("time_limit",3)
    @test disk.density_total_m3 == mem.density_total_m3
    for option in (:xy,:xz,:yz)
        @test forward_psd_saved(folder;kw...,option).psd == forward_psd([a,b,c];kw...,rate_weights_s=[16.,8.,0.],option).psd
    end
    sparse=forward_psd_saved(folder;kw...,storage=:sparse)
    @test all(mem.psd[k...]==v for (k,v) in sparse.psd)
    @test length(sparse.psd)==count(>(0),mem.psd)
    acc=ForwardPSDAccumulator(;kw...)
    observed=[]
    accumulate_forward_psd!(acc,a;rate_weight_s=16.,segment_observer=x->push!(observed,x))
    @test length(observed)==2
    @test observed[1].v0_ms == SA[-2.,0.,0.]
    @test observed[2].v1_ms == SA[2.,0.,0.]
    @test sum(x.t1_s-x.t0_s for x in observed)==acc.residence[1]
    @test all(n==1 for n in values(acc.unique)) # time samples not independent particles
    @test finish_forward_psd(acc).psd == finish_forward_psd(acc).psd # non-mutating
    firstfile=joinpath(folder,"trajectories_00001.jld2")
    @test_throws ArgumentError write_trajectory_batch(firstfile,[a];particle_ids=[1],rate_weights_s=[1.])
    @test_throws ArgumentError foreach_saved_trajectory(_->nothing,[firstfile,firstfile])
    @test_throws ArgumentError forward_psd_saved(folder;kw...,species="H+")
    @test_throws ArgumentError write_trajectory_batch(joinpath(folder,"bad.jld2"),[a];particle_ids=[1],rate_weights_s=[-1.])
    @test !isfile(joinpath(folder,"bad.jld2"))
    incomplete=joinpath(folder,"incomplete.jld2")
    jldopen(incomplete,"w") do f;f["format_version"]=1;end
    @test_throws ArgumentError foreach_saved_trajectory(_->nothing,[incomplete])
    legacy=joinpath(folder,"legacy.jld2")
    jldopen(legacy,"w") do f
        f["p9/state"]=vcat(permutedims(a.t),reduce(hcat,a.u));f["p9/rate_weight_s1"]=16.
    end
    @test forward_psd_saved([legacy];kw...).psd==forward_psd(a;kw...,rate_weights_s=[16.]).psd
    compressed=joinpath(folder,"compressed.jld2")
    write_trajectory_batch(compressed,[a];particle_ids=[10],rate_weights_s=[16.],compress=true)
    @test forward_psd_saved([compressed];kw...).psd==forward_psd(a;kw...,rate_weights_s=[16.]).psd
    println("IO test artifacts: $folder")
end
