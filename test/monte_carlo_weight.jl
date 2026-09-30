using Random, LinearAlgebra
@testset "Maxwellian importance sampling and outward flux" begin
    mass = MarsTP.TP.SpeciesDict["O2+"].m
    @test thermal_speed_from_temperature_ev(1.,mass)^2 ≈ 2*1.602176634e-19/mass
    @test maxwellian_importance_weight_3d((0.,0.,0.),(100.,200.,300.),1.,1.) == 1
    @test_throws ArgumentError thermal_speed_from_temperature_ev(NaN)
    @test_throws ArgumentError particle_density_weight(1.,1.,0.)
    args = (;position_m=(Rm+400e3,0.,0.),bulk_velocity_m_s=(0.,0.,0.),temperature_ev=1.)
    a = sample_maxwellian_source(100; args...,rng=MersenneTwister(4))
    b = sample_maxwellian_source(100; args...,rng=MersenneTwister(4))
    @test a.initial_states == b.initial_states
    @test a.density_weights_m3 === nothing
    @test all(==(1),a.macro_weights)
    for factor in (1.,4.)
        n,A,N = 1e6,2e6,40000
        s = sample_maxwellian_source(N;args...,rng=MersenneTwister(42),
            weights=MonteCarloWeight(sampling_temperature_factor=factor,source_number_density_m3=n),
            normal=(1.,0.,0.),area_m2=A,flux_model=:reservoir)
        @test sum(s.density_weights_m3) ≈ n
        @test abs(sum(s.importance_weights)/N-1) < 0.025
        sigma = thermal_speed_from_temperature_ev(1.,mass)/sqrt(2)
        @test isapprox(sum(s.rate_weights_s),n*A*sigma/sqrt(2pi);rtol=0.03)
        v2 = sum(s.importance_weights[i]*sum(abs2,s.initial_states[i][4:6]) for i in 1:N)/N
        @test isapprox(v2,3sigma^2;rtol=0.03)
        @test all(i -> s.initial_states[i][4]>0 || s.rate_weights_s[i]==0,1:N)
    end
end

@testset "Prescribed bulk-speed source rate" begin
    n,A,N=2e6,3e4,20000
    U=(-1000.,200.,300.)
    args=(;position_m=(Rm+500e3,0.,0.),bulk_velocity_m_s=U,temperature_ev=1.,
        weights=MonteCarloWeight(sampling_temperature_factor=4.,source_number_density_m3=n),area_m2=A,flux_model=:bulk_speed)
    s=sample_maxwellian_source(N;args...,rng=MersenneTwister(18))
    @test sum(s.rate_weights_s) ≈ n*A*norm(U)
    @test s.rate_weights_s ≈ A*norm(U).*s.density_weights_m3
    @test any(v->v[4]<0,s.initial_states) && any(v->v[4]>0,s.initial_states)
    @test all(>(0),s.rate_weights_s)
    # Changing the normal must not change states or bulk-speed weights.
    b=sample_maxwellian_source(N;args...,normal=(-1.,0.,0.),rng=MersenneTwister(18))
    @test s.initial_states==b.initial_states
    @test s.rate_weights_s==b.rate_weights_s
    z=sample_maxwellian_source(100;args...,bulk_velocity_m_s=(0.,0.,0.),rng=MersenneTwister(18))
    @test all(iszero,z.rate_weights_s)
    @test_throws ArgumentError sample_maxwellian_source(10;args...,flux_model=:invalid)
end
@testset "Reservoir drift and zero-field slab density recovery" begin
    mass=MarsTP.TP.SpeciesDict["O2+"].m
    sigma=thermal_speed_from_temperature_ev(1.,mass)/sqrt(2)
    n,A,N=1e6,2e6,100000
    phi(z)=exp(-z*z/2)/sqrt(2pi)
    Phi(z)=0.5*ccall((:erfc,Base.Math.libm),Cdouble,(Cdouble,),-z/sqrt(2))
    for factor in (1.,4.), ur in (-sigma,0.,sigma), ut in (0.,3sigma)
        s=sample_maxwellian_source(N;position_m=(0.,0.,0.),
            bulk_velocity_m_s=(ur,ut,0.),temperature_ev=1.,normal=(1.,0.,0.),area_m2=A,
            weights=MonteCarloWeight(source_number_density_m3=n,sampling_temperature_factor=factor),
            rng=MersenneTwister(20260929))
        expected=n*A*(sigma*phi(ur/sigma)+ur*Phi(ur/sigma))
        @test isapprox(sum(s.rate_weights_s),expected;rtol=.035)
        # In a zero-field slab of thickness L, tau=L/vr for outward draws.
        # Q*tau/(A*L) must recover the outgoing Maxwellian half-space density.
        outgoing=findall(i->s.initial_states[i][4]>0,1:N)
        contributions=[s.rate_weights_s[i]/(A*s.initial_states[i][4]) for i in outgoing]
        @test isapprox(sum(contributions),n*Phi(ur/sigma);rtol=.035)
        recovered_tangent=sum(contributions[k]*s.initial_states[i][5] for (k,i) in enumerate(outgoing))/sum(contributions)
        @test abs(recovered_tangent-ut)<.035sigma
        recovered_second=sum(contributions[k]*(s.initial_states[i][5]-ut)^2 for (k,i) in enumerate(outgoing))/sum(contributions)
        @test isapprox(recovered_second,sigma^2;rtol=.04)
        @test any(iszero,s.rate_weights_s)
        @test all(i->s.initial_states[i][4]>0 || iszero(s.rate_weights_s[i]),1:N)
    end
    @test_throws ArgumentError sample_maxwellian_source(10;position_m=(0.,0.,0.),
        bulk_velocity_m_s=(0.,0.,0.),temperature_ev=1.,area_m2=A)
end
