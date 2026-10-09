# Small actual-MHD comparison. No data or previous results are overwritten.
using MarsTP, StaticArrays, KernelAbstractions, CUDA, Test, LinearAlgebra
CUDA.functional() || error("A working CUDA device is required")
fields = load_mhd_fields()
cfg = ForwardTraceConfig(species="O2+",solver=:boris,dt=.01,tspan=(0.,10.))
inner,outer = Rm+200e3,Router
# Include near-boundary particles and a particle farther inside the domain.
states = [SA[inner+100.,0.,0.,-1e4,0.,0.],
          SA[outer-100.,0.,0.,1e4,0.,0.],
          SA[1.5Rm,.1Rm,.4Rm,1e4,0.,0.]]
cpu = trace_forward_bounded(states;config=cfg,fields,backend=CPU(),save_every=10)
gpu = trace_forward_bounded(states;config=cfg,fields,backend=CUDA.CUDABackend(),save_every=10)
@test getproperty.(cpu,:status) == getproperty.(gpu,:status)
@test getproperty.(gpu,:status) == [:inner,:outer,:time_limit]
for (a,b) in zip(cpu,gpu)
    @test a.t ≈ b.t atol=1e-8
    @test maximum(norm(x[1:3]-y[1:3]) for (x,y) in zip(a.u,b.u)) < .01
    @test maximum(norm(x[4:6]-y[4:6]) for (x,y) in zip(a.u,b.u)) < .01
    println((status=b.status,time_s=last(b.t),altitude_km=(norm(last(b.u)[1:3])-Rm)/1e3))
end
