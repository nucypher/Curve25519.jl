using DarkIntegers
using Curve25519: base_point, curve_order, RistrettoPointT
using Curve25519
using BenchmarkTools

b = base_point()
b2 = b + b

order = curve_order()

z = zero(RistrettoPointT)

x1 = b * (order - 1)
x2 = b2 * (order ÷ 2)
x3 = b * order
x4 = b * (order + 1)

@assert x1 == x2
@assert x1 + b == x3
@assert x3 == z
@assert x4 == b


display(@benchmark b * order)
println()
