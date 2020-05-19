using Jute
using Curve25519
using BenchmarkTools
using DarkCurves


function ref_mul(p, y)
    res = p
    for i in 2:y
        res += p
    end
    res
end


@testcase "Basic operations" begin

    stp = curve_scalar_type(RistrettoCurveVT)
    ptp = curve_point_type(RistrettoCurveVT)

    b = one(ptp)

    @test ref_mul(b + b, 23) == (b + b) * convert(stp, 23)
    @test b * zero(stp) == zero(ptp)
    @test iszero(b * zero(stp))
    @test b * one(stp) == b
    @test iszero(b * (-one(stp)) + b)
end


@testcase tags=[:performance] "Addition performance" begin

    ptp = curve_point_type(RistrettoCurveVT)

    b1 = one(ptp)
    b2 = b1 + b1
    b4 = b2 + b2

    trial = @benchmark $b2 + $b4
    @test_result benchmark_result(trial)
end


@testcase "Multiplication" begin

    stp = curve_scalar_type(RistrettoCurveVT)
    ptp = curve_point_type(RistrettoCurveVT)

    b1 = one(ptp)
    p = b1 + b1

    x_i = 123
    x = convert(stp, x_i)

    @test p * x == ref_mul(p, x_i)
end


@testcase tags=[:performance] "Multiplication performance" begin

    stp = curve_scalar_type(RistrettoCurveVT)
    ptp = curve_point_type(RistrettoCurveVT)

    b1 = one(ptp)
    b2 = b1 + b1

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(stp, x_bi)

    trial = @benchmark $b2 * $x
    @test_result benchmark_result(trial)
end


@testcase "Random scalar" for rng in fixed_rng(123)

    stp = curve_scalar_type(RistrettoCurveVT)
    x = rand(rng, stp)
    @test typeof(x) == stp

    a = rand(rng, stp, 10)
    @test eltype(a) == stp
end


@testcase "Random point" for rng in fixed_rng(123)

    ptp = curve_point_type(RistrettoCurveVT)
    x = rand(rng, ptp)
    @test typeof(x) == ptp

    a = rand(rng, ptp, 10)
    @test eltype(a) == ptp
end


@testcase tags=[:performance] "Random scalar, performance" for rng in fixed_rng(123)
    stp = curve_scalar_type(RistrettoCurveVT)
    trial = @benchmark rand($rng, $stp)
    @test_result benchmark_result(trial)
end


@testcase tags=[:performance] "Random point, performance" for rng in fixed_rng(123)
    ptp = curve_point_type(RistrettoCurveVT)
    trial = @benchmark rand($rng, $ptp)
    @test_result benchmark_result(trial)
end


exit(runtests())
