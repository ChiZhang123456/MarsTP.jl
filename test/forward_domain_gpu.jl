using MarsTP, StaticArrays, Test, CUDA, KernelAbstractions, LinearAlgebra

@testset "Mixed radial exits and interior NaNs on CPU/GPU" begin
    r = collect(range(Rinner,Router;length=8))
    th = collect(range(0.,pi;length=9))
    ph = collect(range(0.,2pi;length=17))
    E = zeros(3,length(r),length(th),length(ph))
    fields = MHDFields(r,th,ph,E,copy(E),:total,"zero_field_domain_test")
    invalid = MHDFields(r,th,ph,fill(NaN,size(E)),copy(E),:total,"invalid_interior")
    cfg = ForwardTraceConfig(species="H+",solver=:boris,dt=.01,tspan=(0.,.1))
    states = [SA[Router-1500.,0.,0.,1e5,0.,0.],
              SA[Rinner+1500.,0.,0.,-1e5,0.,0.],
              SA[2Rm,0.,0.,0.,0.,0.]]
    backends = CUDA.functional() ? (nothing,CPU(),CUDA.CUDABackend()) : (nothing,CPU())
    for backend in backends, saveat in ((),.03)
        sols = only.(getproperty.(trace_forward(states;config=cfg,fields,backend,saveat),:u))
        @test getproperty.(sols,:retcode) == [MarsTP.TP.ReturnCode.Terminated,
            MarsTP.TP.ReturnCode.Terminated,MarsTP.TP.ReturnCode.Success]
        @test last.(getproperty.(sols,:t)) ≈ [.01,.01,.1]
        @test all(s -> all(>(0),diff(s.t)),sols)
        @test all(s -> all(u -> all(isfinite,u),s.u),sols)
        @test all(s -> all(u -> Rinner <= norm(u[1:3]) <= Router,s.u),sols)
        for (start,sol) in zip(states,sols)
            @test last(sol.u)[1:3] ≈ start[1:3]+last(sol.t)*start[4:6]
            @test last(sol.u)[4:6] ≈ start[4:6]
        end
        immediate = only(only(trace_forward([SA[Router-1.,0.,0.,1e5,0.,0.]];
            config=cfg,fields,backend,saveat)).u)
        @test immediate.retcode == MarsTP.TP.ReturnCode.Terminated
        @test immediate.t == [0.]
        @test all(u -> all(isfinite,u),immediate.u)
        @test_throws ErrorException trace_forward([states[3]];config=cfg,fields=invalid,backend,saveat)
        low = SA[Rm+99e3,0.,0.,-1e5,0.,0.]
        mixed = only.(getproperty.(trace_forward([low,states[3]];config=cfg,fields,backend,saveat),:u))
        @test mixed[1].retcode == MarsTP.TP.ReturnCode.Terminated
        @test mixed[1].t == [0.]
        @test mixed[1].u == [low]
        @test mixed[2].retcode == MarsTP.TP.ReturnCode.Success
        @test last(mixed[2].t) ≈ .1
    end
    # Adaptive Boris needs a nonzero magnetic field to set its gyro timestep.
    B = copy(E)
    B[1,:,:,:] .= 1e-8
    MarsTP._rotate_vectors_to_spherical!(B,r,th,ph)
    magnetic = MHDFields(r,th,ph,E,B,:total,"uniform_Bx_domain_test")
    adaptive = ForwardTraceConfig(species="H+",solver=:adaptive,safety=.001,tspan=(0.,.1))
    for saveat in ((),.03)
        sols = only.(getproperty.(trace_forward(states;config=adaptive,fields=magnetic,saveat),:u))
        @test getproperty.(sols,:retcode) == [MarsTP.TP.ReturnCode.Terminated,
            MarsTP.TP.ReturnCode.Terminated,MarsTP.TP.ReturnCode.Success]
        @test all(s -> all(u -> all(isfinite,u),s.u),sols)
        @test all(s -> all(>(0),diff(s.t)),sols)
        @test all(s -> all(u -> Rinner <= norm(u[1:3]) <= Router,s.u),sols)
        @test all(s -> last(s.t) < .015,sols[1:2])
    end
end
