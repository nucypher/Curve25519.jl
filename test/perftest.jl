using DarkIntegers
using Curve25519: base_point, curve_order, RistrettoPointCT, RistrettoScalarCT, RistrettoScalarVT
using Curve25519
using BenchmarkTools
using ConstantTime

const CT = ConstantTime


b = base_point(RistrettoPointCT)
b2 = b + b

order = curve_order(RistrettoScalarVT)

z = zero(RistrettoPointCT)

x1 = b * CT.wrap(order - 1)
x2 = b2 * CT.wrap(order ÷ 2)
x3 = b * CT.wrap(order)
x4 = b * CT.wrap(order + 1)

@assert x1 == x2
@assert x1 + b == x3
@assert x3 == z
@assert x4 == b


ct_order = CT.wrap(order)

display(@benchmark b * ct_order)
println()
