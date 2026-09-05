using MarsTP, StaticArrays, Test, LinearAlgebra
const TP = MarsTP.TP
@testset "MHD ionosphere definition" begin
    source = IonosphereSource(_->2., _->1000., (_->-3., _->4., _->0.), Rm+400e3)
    p = ionosphere_properties(source, SA[1.,0.,0.])
    @test p.flux == 10 # no radial selection, despite negative Ux
    @test p.Ui == SA[-3.,4.,0.]
    sp = TP.SpeciesDict["O2+"]
    vth = sqrt(2TP.kB*1000/sp.m)
    @test ionosphere_distribution(source,SA[1.,0.,0.],p.Ui).f ≈ 2/(pi^1.5*vth^3)
    @test_throws ArgumentError MarsTP._ionosphere_radius(199.)
    @test_throws ArgumentError MarsTP._ionosphere_radius(801.)
    @test MarsTP._ionosphere_radius(200.) == Rinner
    @test MarsTP._ionosphere_radius(800.) == Rm+800e3
end
@testset "Backtrace crossings and source quadrature" begin
    param = TP.prepare(_->SA[0.,0.,0.], _->SA[0.,0.,1e-9]; species=TP.SpeciesDict["O2+"])
    # Motion parallel to B: exact straight line. Constant source tests full last segment.
    for alg in (:boris,:adaptive), h in (200.,400.,800.)
        cfg=BacktraceConfig(ionosphere_altitude_km=h,solver=alg,dt=-.2,tspan=(0.,-2.))
        r=Rm+h*1e3
        result=MarsTP._trace_sources(SA[0.,0.,r+1000],SA[0.,0.,1000.],param,cfg,
            (p,v)->3.,(p,v)->(;f=7.,flux=8.))
        @test result.status == 2
        @test result.time ≈ -1 atol=1e-7
        @test norm(result.position) ≈ r
        @test result.volume ≈ 3 atol=1e-7
        @test result.ionosphere == 7
        @test result.velocity ≈ SA[0.,0.,1000.]
    end
    cfg=BacktraceConfig(tspan=(0.,-.93),dt=-.2)
    result=MarsTP._trace_sources(SA[0.,0.,2Rm],SA[0.,0.,1000.],param,cfg,(p,v)->3.,(p,v)->(;f=7.,flux=8.))
    @test result.status == 1
    @test result.time ≈ -.93
    @test result.volume ≈ 2.79
    @test result.ionosphere == 0
    result=MarsTP._trace_sources(SA[0.,0.,Router-100],SA[0.,0.,-1000.],param,cfg,(p,v)->3.,(p,v)->(;f=7.,flux=8.))
    @test result.status == 4
    @test result.time ≈ -.1 atol=1e-8
    @test result.volume ≈ .3 atol=1e-8
    @test MarsTP._backtrace_crossing(SA[2.,0.,0.],SA[-2.,0.,0.],1.,4.) == (.25,2)
end

@testset "Synchronized boundary velocity with electric acceleration" begin
    sp=TP.SpeciesDict["O2+"]
    acceleration=100.
    param=TP.prepare(_->SA[0.,0.,acceleration*sp.m/sp.q], _->SA[0.,0.,1e-9]; species=sp)
    radius=Rm+400e3
    exact_time=(-1000+sqrt(1000^2-2acceleration*1000))/acceleration
    errors=Float64[]
    for dt in (-.1,-.05)
        cfg=BacktraceConfig(ionosphere_altitude_km=400.,dt=dt,tspan=(0.,-2.))
        r=MarsTP._trace_sources(SA[0.,0.,radius+1000],SA[0.,0.,1000.],param,cfg,
            (p,v)->0.,(p,v)->(;f=v[3],flux=0.))
        @test r.status==2
        @test r.velocity[3] ≈ 1000+acceleration*r.time atol=1e-7
        @test r.ionosphere==r.velocity[3]
        push!(errors,abs(r.time-exact_time))
    end
    @test errors[2] < errors[1]
    @test errors[2] < 1e-4
end
