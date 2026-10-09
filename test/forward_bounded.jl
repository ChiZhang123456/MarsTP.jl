using KernelAbstractions, LinearAlgebra

function bounded_checks(backend)
    inner,outer = Rm+200e3,Router
    r = collect(range(inner,outer;length=8))
    th = collect(range(0.,pi;length=9))
    ph = collect(range(0.,2pi;length=17))
    E = zeros(3,length(r),length(th),length(ph))
    fields = MHDFields(r,th,ph,E,copy(E),:total,"zero_field")
    cfg = ForwardTraceConfig(species="H+",solver=:boris,dt=1.,tspan=(5.,15.))
    states = [SA[inner+250.,0.,0.,-100.,0.,0.],
              SA[outer-750.,0.,0.,100.,0.,0.],
              SA[inner+1000.,0.,0.,0.,0.,0.],
              SA[inner,0.,0.,-100.,0.,0.],
              SA[outer,0.,0.,100.,0.,0.],
              SA[inner,0.,0.,100.,0.,0.]]
    sols = trace_forward_bounded(states;config=cfg,fields,backend,save_every=3)
    @test getproperty.(sols,:status) == [:inner,:outer,:time_limit,:inner,:outer,:time_limit]
    @test last.(getproperty.(sols,:t)) ≈ [7.5,12.5,15.,5.,5.,15.]
    for (i,sol) in enumerate(sols)
        @test all(isfinite,s for u in sol.u for s in u)
        @test all(>(0),diff(sol.t))
        @test last(sol.u)[4:6] ≈ states[i][4:6]
        @test last(sol.u)[1:3] ≈ states[i][1:3]+(last(sol.t)-5)*states[i][4:6]
    end
    @test norm(last(sols[1].u)[1:3]) ≈ inner
    @test norm(last(sols[2].u)[1:3]) ≈ outer
    rounded_r = copy(r)
    rounded_r[end] = outer-8eps(outer)
    rounded = MHDFields(rounded_r,th,ph,E,copy(E),:total,"rounded_outer_axis")
    on_outer = trace_forward_bounded([states[5]];config=cfg,fields=rounded,backend)[1]
    @test on_outer.status == :outer
    @test on_outer.t == [5.]
    # Cross the entire inner sphere in one step; endpoint-only checks miss it.
    big = ForwardTraceConfig(species="H+",solver=:boris,dt=1.,tspan=(0.,1.))
    through = trace_forward_bounded([SA[inner+100.,0.,0.,-2inner-200.,0.,0.]];
        config=big,fields,backend,save_every=0)[1]
    @test through.status == :inner
    @test last(through.t) ≈ 100/(2inner+200)
    @test last(through.u)[1] ≈ inner
    # Nonzero fields: unbounded upstream Boris and bounded kernel agree before events.
    B = copy(E)
    for k in eachindex(ph),j in eachindex(th),i in eachindex(r)
        B[:,i,j,k] .= [0.,0.,1e-8]
    end
    MarsTP._rotate_vectors_to_spherical!(B,r,th,ph)
    magnetic = MHDFields(r,th,ph,E,B,:total,"uniform_Bz")
    cfgB = ForwardTraceConfig(species="H+",solver=:boris,dt=.01,tspan=(0.,1.))
    start = [SA[1.5Rm,.1Rm,.4Rm,1e4,0.,0.]]
    bounded = trace_forward_bounded(start;config=cfgB,fields=magnetic,backend)[1]
    reference = only(trace_forward(start;config=cfgB,fields=magnetic)[1].u)
    @test bounded.t ≈ reference.t
    @test bounded.u ≈ reference.u
    @test maximum(abs(sum(abs2,s[4:6])/1e8-1) for s in bounded.u) < 1e-10
    @test last(bounded.u)[5] < 0
    # Check a fractional magnetic-field event against the analytic proton orbit.
    event_start = SA[inner+100.,0.,0.,-1000.,0.,0.]
    omega = MarsTP.TP.SpeciesDict["H+"].q/MarsTP.TP.SpeciesDict["H+"].m*1e-8
    analytic(t) = SA[event_start[1]-1000sin(omega*t)/omega,
                     1000(1-cos(omega*t))/omega,0.]
    lo,hi = 0.,.2
    for _ in 1:60
        mid=(lo+hi)/2
        if norm(analytic(mid)) > inner; lo=mid; else; hi=mid; end
    end
    exact = (lo+hi)/2
    errors = Float64[]
    for dt in (.02,.01)
        c=ForwardTraceConfig(species="H+",solver=:boris,dt=dt,tspan=(0.,.2))
        event=trace_forward_bounded([event_start];config=c,fields=magnetic,backend)[1]
        @test event.status == :inner
        @test norm(last(event.u)[1:3]) ≈ inner
        push!(errors,abs(last(event.t)-exact))
    end
    @test errors[2] < errors[1]
    # Invalid fields remain a numerical failure, never an escape or absorption.
    invalid = MHDFields(r,th,ph,fill(NaN,size(E)),copy(E),:total,"invalid")
    @test_throws ErrorException trace_forward_bounded(start;config=cfgB,fields=invalid,backend)
    failed = trace_forward_bounded(start;config=cfgB,fields=invalid,backend,throw_on_failure=false)[1]
    @test failed.status == :numerical_failure
    @test failed.retcode == MarsTP.TP.ReturnCode.Failure
    @test_throws ArgumentError trace_forward_bounded(states;fields,backend,config=ForwardTraceConfig())
    @test_throws ArgumentError trace_forward_bounded(states;fields,backend,config=cfg,inner_radius_m=inner-1)
end

function bounded_work_checks(backend)
    r=collect(range(Rinner,Router;length=8));th=collect(range(0.,pi;length=9));ph=collect(range(0.,2pi;length=17))
    E=zeros(3,length(r),length(th),length(ph))
    E[1,:,:,:].=1e-5
    MarsTP._rotate_vectors_to_spherical!(E,r,th,ph)
    fields=MHDFields(r,th,ph,E,zero(E),:total,"uniform_Ex")
    interp(A)=MarsTP.TP.build_interpolator(MarsTP.TP.StructuredGrid,A,r,th,ph)
    itp=MarsTP.FieldWorkInterpolators(interp(E),interp(2E),interp(-E))
    cfg=ForwardTraceConfig(species="O2+",solver=:boris,dt=.1,tspan=(0.,1.))
    states=[SA[Rinner+150.,0.,0.,-1000.,0.,0.],SA[Router-150.,0.,0.,1000.,0.,0.],
        SA[Rm+500e3,0.,0.,1000.,0.,0.]]
    full=trace_forward_bounded(states;config=cfg,fields,backend,work_itp=itp)
    sparse=trace_forward_bounded(states;config=cfg,fields,backend,work_itp=itp,save_every=0)
    @test getproperty.(full,:status)==[:inner,:outer,:time_limit]
    for (a,b) in zip(full,sparse)
        reference=trajectory_work(a,itp;species="O2+",mode=:summary).summary
        for key in keys(reference)
            @test getproperty(a.work,key) ≈ getproperty(reference,key) atol=1e-8
            @test getproperty(a.work,key) ≈ getproperty(b.work,key) atol=1e-8
        end
        @test a.work.conv_eV ≈ 2a.work.total_eV
        @test a.work.hall_eV ≈ -a.work.total_eV
        @test abs(a.work.energy_residual_eV)<1e-8
        @test length(b.u)==2
    end
    @test full[1].work.total_eV<0
    @test full[2].work.total_eV>0
end

@testset "Per-particle spherical stopping (CPU kernel)" begin
    bounded_checks(KernelAbstractions.CPU())
    bounded_work_checks(KernelAbstractions.CPU())
end

# Run explicitly when validating a CUDA implementation; ordinary tests remain
# usable on machines without an NVIDIA device.
if get(ENV,"MARSTP_TEST_CUDA","false") == "true"
    using CUDA
    CUDA.functional() || error("CUDA validation requested but GPU is unavailable")
    @testset "Per-particle spherical stopping (CUDA)" begin
        bounded_checks(CUDA.CUDABackend())
        bounded_work_checks(CUDA.CUDABackend())
    end
end
