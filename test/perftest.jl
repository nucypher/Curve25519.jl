using Random
using DarkIntegers
using Curve25519: base_point, curve_order, RistrettoPointCT, RistrettoScalarCT, RistrettoScalarVT, RistrettoPointVT
using Curve25519
using BenchmarkTools
using ConstantTime
using DarkCurves
using Profile

const CT = ConstantTime


function test_single_mul_25519()
    b = base_point(RistrettoPointCT)
    b2 = b + b

    order = curve_order(RistrettoScalarVT)

    z = zero(RistrettoPointCT)

    x1 = b * (order - 1)
    x2 = b2 * (order ÷ 2)
    x3 = b * (order)
    x4 = b * (order + 1)

    @assert x1 == x2
    @assert x1 + b == x3
    @assert x3 == z
    @assert x4 == b


    display(@benchmark $b * $order)
    println()
end


function test_single_mul_secp()

    curve_type = Curve_secp256k1
    point_type = JacobianPoint
    stp = curve_scalar_type(curve_type, MgModUInt, MLUInt{2, UInt128})
    ptp = point_type{curve_type, stp}

    b1 = one(ptp)

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325686
    x = convert(MLUInt{2, UInt128}, x_bi)

    display(@benchmark $b1 * $x)
    println()
end


function test_rand_25519()
    base = base_point(RistrettoPointVT)
    table = Curve25519.RistrettoBasepointTable(base)

    b1 = base * 123
    b2 = table * 123
    @assert b1 == b2

    x_bi = 115047236638587805833081834189719086745649315857841928574581145752217906325
    x = convert(RistrettoScalarVT, x_bi)

    @assert table * 0 == zero(RistrettoPointVT)

    display(@benchmark $base * $x)
    println()

    display(@benchmark $table * $x)
    println()

    rng = MersenneTwister(123)
    rand(rng, RistrettoPointVT)

    display(@benchmark rand($rng, RistrettoPointVT))
    println()
end


#test_single_mul_25519()
#test_single_mul_secp()
test_rand_mul_25519()
