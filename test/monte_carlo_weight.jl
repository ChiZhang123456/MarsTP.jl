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
            normal=(1.,0.,0.),area_m2=A)
        @test sum(s.density_weights_m3) ≈ n
        @test abs(sum(s.importance_weights)/N-1) < 0.025
        sigma = thermal_speed_from_temperature_ev(1.,mass)/sqrt(2)
        @test isapprox(sum(s.rate_weights_s),n*A*sigma/sqrt(2pi);rtol=0.03)
        v2 = sum(s.importance_weights[i]*sum(abs2,s.initial_states[i][4:6]) for i in 1:N)/N
        @test isapprox(v2,3sigma^2;rtol=0.03)
        @test all(i -> s.initial_states[i][4]>0 || s.rate_weights_s[i]==0,1:N)
    end
end
