using DarkIntegers
using Curve25519: RistrettoScalar, base_point, RistrettoPoint, curve_order, EdwardsPoint
using Curve25519
using BenchmarkTools


b = base_point()
b2 = b + b

order = curve_order()

z = zero(RistrettoPoint{RistrettoScalar})

x1 = b * (order - 1)
x2 = b2 * (order ÷ 2)
x3 = b * order
x4 = b * (order + 1)

@assert x1 == x2
@assert x1 + b == z
@assert x3 == z
@assert x4 == b


display(@benchmark b * order)
println()
