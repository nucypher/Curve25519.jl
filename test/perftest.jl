using DarkIntegers
using Curve25519: RistrettoScalar, base_point, RistrettoPoint, curve_order, EdwardsPoint, FieldElement51
using Curve25519
using BenchmarkTools


function from_fe51(x::FieldElement51)
    xx = big.(x.limbs)
    xx[1] + (xx[2] << 51) + (xx[3] << 102) + (xx[4] << 153) + (xx[5] << 204)
end


b = base_point()
b2 = b + b

order = curve_order()
order_big = from_fe51(order)

z = zero(RistrettoPoint{RistrettoScalar})

x1 = b * (order_big - 1)
x2 = b2 * (order_big ÷ 2)
x3 = b * order_big
x4 = b * (order_big + 1)

@assert x1 == x2
@assert x1 + b == x3
@assert x3 == z
@assert x4 == b


display(@benchmark b * order_big)
println()
