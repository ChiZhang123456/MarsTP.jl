using MarsTP, StaticArrays, Test
@testset "Cartesian electric work" begin
    p = SA[0.0, 2.0, 0.0]
    v = SA[3.0, 0.0, 0.0]
    @test MarsTP._work_one(_ -> SA[1.0, 0.0, 0.0], p, v, 2.0, 1.0) == 6.0
    @test MarsTP._work_one(_ -> SA[-1.0, 0.0, 0.0], p, v, 2.0, 1.0) == -6.0
    @test_throws ErrorException MarsTP._work_one(_ -> SA[NaN, 0.0, 0.0], p, v, 2.0, 1.0)
    fields = MarsTP.FieldWorkInterpolators(_ -> SA[1.0, 0.0, 0.0],
        _ -> SA[2.0, 0.0, 0.0], _ -> SA[-1.0, 0.0, 0.0])
    trajectory = (; t = [0.0, 1.0],
        u = [SA[0., 2., 0., 3., 0., 0.], SA[3., 2., 0., 3., 0., 0.]])
    result = field_work(trajectory, fields)
    @test result.total_eV ≈ 3.0
    @test result.conv_eV ≈ 6.0
    @test result.hall_eV ≈ -3.0
    @test field_work((; u = [trajectory]), fields) == result
end
