#=
Ported from backend/serial/u64/field.rs

Integers modulo 2^255-19, used as coordinates in the curve.
Represented with 5 UInt64 integers, with 51 bits in each.

The limbs can grow up to 54 bits, afterwards arithmetic operations will stop working,
and reduction is needed. Currently reduction only happens during negation (explicitly)
or multiplication (implicitly), and it is assumed that it is enough for all curve operations.
The limb size is not tracked, use with care.
=#

const LOW_51_BIT_MASK = (one(UInt64) << 51) - one(UInt64)


struct FieldElement51
    limbs :: NTuple{5, UInt64}
end


function CT.select(choice::CT.Choice, p::CT.Value{FieldElement51}, q::CT.Value{FieldElement51})
    pp = CT.unwrap(p)
    qq = CT.unwrap(q)
    CT.wrap(FieldElement51((
        CT.unwrap(CT.select(choice, CT.wrap(pp.limbs[1]), CT.wrap(qq.limbs[1]))),
        CT.unwrap(CT.select(choice, CT.wrap(pp.limbs[2]), CT.wrap(qq.limbs[2]))),
        CT.unwrap(CT.select(choice, CT.wrap(pp.limbs[3]), CT.wrap(qq.limbs[3]))),
        CT.unwrap(CT.select(choice, CT.wrap(pp.limbs[4]), CT.wrap(qq.limbs[4]))),
        CT.unwrap(CT.select(choice, CT.wrap(pp.limbs[5]), CT.wrap(qq.limbs[5]))),
        )))
end


Base.getindex(x::FieldElement51, i) = x.limbs[i]


Base.zero(::Type{FieldElement51}) = FieldElement51((0, 0, 0, 0, 0))


Base.one(::Type{FieldElement51}) = FieldElement51((1, 0, 0, 0, 0))


Base.:+(x::FieldElement51, y::FieldElement51) =
    FieldElement51((x[1] + y[1], x[2] + y[2], x[3] + y[3], x[4] + y[4], x[5] + y[5]))


function Base.:-(x::FieldElement51, y::FieldElement51)
    # Assuming that each limb is at most 54 bits in length,
    # and adding a corresponding multiple of the modulus before subtraction.
    weak_reduce(FieldElement51((
        (x[1] + 0x007ffffffffffed0) - y[1],
        (x[2] + 0x007ffffffffffff0) - y[2],
        (x[3] + 0x007ffffffffffff0) - y[3],
        (x[4] + 0x007ffffffffffff0) - y[4],
        (x[5] + 0x007ffffffffffff0) - y[5],
        )))
end


function Base.:-(x::FieldElement51)
    # See the comment in the binary `-` method
    weak_reduce(FieldElement51((
        0x007ffffffffffed0 - x[1],
        0x007ffffffffffff0 - x[2],
        0x007ffffffffffff0 - x[3],
        0x007ffffffffffff0 - x[4],
        0x007ffffffffffff0 - x[5],
        )))
end


# Given 64-bit input limbs, reduce to enforce the bound 2^(51 + epsilon).
function weak_reduce(x::FieldElement51)
    #=
    Since the input limbs are bounded by 2^64, the biggest
    carry-out is bounded by 2^13.

    The biggest carry-in is c5 * 19, resulting in

        2^51 + 19*2^13 < 2^51.0000000001

    Because we don't need to canonicalize, only to reduce the
    limb sizes, it's OK to do a "weak reduction", where we
    compute the carry-outs in parallel.
    =#

    c1 = x[1] >> 51
    c2 = x[2] >> 51
    c3 = x[3] >> 51
    c4 = x[4] >> 51
    c5 = x[5] >> 51

    x1 = x[1] & LOW_51_BIT_MASK
    x2 = x[2] & LOW_51_BIT_MASK
    x3 = x[3] & LOW_51_BIT_MASK
    x4 = x[4] & LOW_51_BIT_MASK
    x5 = x[5] & LOW_51_BIT_MASK

    x1 += c5 * 19
    x2 += c1
    x3 += c2
    x4 += c3
    x5 += c4

    FieldElement51((x1, x2, x3, x4, x5))
end


function Base.:*(x::FieldElement51, y::FieldElement51)

    # Alias self, _rhs for more readable formulas
    a = x
    b = y

    #=
    Precondition: assume input limbs x[i], y[i] are bounded as

        x[i], y[i] < 2^(51 + e)

    where `e` is a real parameter measuring the "bit excess" of the limbs.
    =#

    #=
    This fits into a u64 whenever 51 + e + lg(19) < 64.
    Since 51 + e + lg(19) < 51 + 4.25 + e = 55.25 + e,
    this fits if e < 8.75.
    =#
    b1_19 = y[2] * 19
    b2_19 = y[3] * 19
    b3_19 = y[4] * 19
    b4_19 = y[5] * 19

    # Multiply to get 128-bit coefficients of output
    c0 = (
        widemul(x[1],y[1]) + widemul(x[5],b1_19) + widemul(x[4],b2_19) +
        widemul(x[3],b3_19) + widemul(x[2],b4_19))
    c1 = (
        widemul(x[2],y[1]) + widemul(x[1],y[2])  + widemul(x[5],b2_19) +
        widemul(x[4],b3_19) + widemul(x[3],b4_19))
    c2 = (
        widemul(x[3],y[1]) + widemul(x[2],y[2]) + widemul(x[1],y[3]) +
        widemul(x[5],b3_19) + widemul(x[4],b4_19))
    c3 = (
        widemul(x[4],y[1]) + widemul(x[3],y[2]) + widemul(x[2],y[3]) +
        widemul(x[1],y[4]) + widemul(x[5],b4_19))
    c4 = (
        widemul(x[5],y[1]) + widemul(x[4],y[2]) + widemul(x[3],y[3]) +
        widemul(x[2],y[4]) + widemul(x[1],y[5]))

    #=
    How big are the c[i]? We have

        c[i] < 2^(102 + 2*e) * (1+i + (4-i)*19)
             < 2^(102 + lg(1 + 4*19) + 2*e)
             < 2^(108.27 + 2*e)

    The carry (c[i] >> 51) fits into a u64 when
        108.27 + 2*e - 51 < 64
        2*e < 6.73
        e < 3.365.

    So we require e <= 3 to ensure this fits.
    =#

    c1 += (c0 >> 51) % UInt64
    out0 = (c0 % UInt64) & LOW_51_BIT_MASK

    c2 += (c1 >> 51) % UInt64
    out1 = (c1 % UInt64) & LOW_51_BIT_MASK

    c3 += (c2 >> 51) % UInt64
    out2 = (c2 % UInt64) & LOW_51_BIT_MASK

    c4 += (c3 >> 51) % UInt64
    out3 = (c3 % UInt64) & LOW_51_BIT_MASK

    carry = (c4 >> 51) % UInt64
    out4 = (c4 % UInt64) & LOW_51_BIT_MASK

    #=
    To see that this does not overflow, we need out[0] + carry * 19 < 2^64.

        c4 < x1*y5 + x2*y4 + x3*y3 + x4*y2 + x5*y1 + (carry from c3)
           < 5*(2^(51 + e) * 2^(51 + e)) + (carry from c3)
           < 2^(102 + 2*e + lg(5)) + 2^64.

    When e <= 3 we get

        c4 < 2^110.33  so that carry < 2^59.33

    so that

        out[0] + carry * 19 < 2^51 + 19 * 2^59.33 < 2^63.58

    and there is no overflow.
    =#

    out0 += carry * 19

    # Now out[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
    out1 += out0 >> 51
    out0 &= LOW_51_BIT_MASK

    # Now out[i] < 2^(51 + epsilon) for all i.
    FieldElement51((out0, out1, out2, out3, out4))
end


Base.convert(::Type{FieldElement51}, x::Integer) = FieldElement51((
    (x & 0x0007ffffffffffff) % UInt64,
    ((x >> 51) & 0x0007ffffffffffff) % UInt64,
    ((x >> 102) & 0x0007ffffffffffff) % UInt64,
    ((x >> 153) & 0x0007ffffffffffff) % UInt64,
    ((x >> 204) & 0x0007ffffffffffff) % UInt64,
    ))


function Base.convert(::Type{T}, x::FieldElement51) where T <: Integer
    xx = convert.(T, x.limbs)
    xx[1] + (xx[2] << 51) + (xx[3] << 102) + (xx[4] << 153) + (xx[5] << 204)
end


# Arithmetic operations on FieldElement51 are constant-time


Base.:+(x::CT.Value{FieldElement51}, y::CT.Value{FieldElement51}) = CT.Value(CT.unwrap(x) + CT.unwrap(y))

Base.:-(x::CT.Value{FieldElement51}) = CT.Value(-CT.unwrap(x))
Base.:-(x::CT.Value{FieldElement51}, y::CT.Value{FieldElement51}) = CT.Value(CT.unwrap(x) - CT.unwrap(y))

Base.:*(x::CT.Value{FieldElement51}, y::CT.Value{FieldElement51}) = CT.Value(CT.unwrap(x) * CT.unwrap(y))
Base.:*(x::CT.Value{FieldElement51}, y::FieldElement51) = CT.Value(CT.unwrap(x) * y)


