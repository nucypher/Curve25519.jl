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


@inline function CT.select(choice::CT.Choice, p::CT.Value{FieldElement51}, q::CT.Value{FieldElement51})
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


@inline Base.getindex(x::FieldElement51, i) = x.limbs[i]


@inline Base.zero(::Type{FieldElement51}) = FieldElement51((0, 0, 0, 0, 0))


@inline Base.one(::Type{FieldElement51}) = FieldElement51((1, 0, 0, 0, 0))


@inline Base.:+(x::FieldElement51, y::FieldElement51) =
    FieldElement51((x[1] + y[1], x[2] + y[2], x[3] + y[3], x[4] + y[4], x[5] + y[5]))


@inline function Base.:-(x::FieldElement51, y::FieldElement51)
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


@inline function Base.:-(x::FieldElement51)
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
@inline function weak_reduce(x::FieldElement51)
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


@inline function Base.:*(x::FieldElement51, y::FieldElement51)

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


#/// Given `k > 0`, return `self^(2^k)`.
function pow2k(x::FieldElement51, k::Int)

    #/// Multiply two 64-bit integers with 128 bits of output.
    a = x

    a0 = a[1]
    a1 = a[2]
    a2 = a[3]
    a3 = a[4]
    a4 = a[5]

    for i in 1:k
        #=
        // Precondition: assume input limbs a[i] are bounded as
        //
        // a[i] < 2^(51 + b)
        //
        // where b is a real parameter measuring the "bit excess" of the limbs.

        // Precomputation: 64-bit multiply by 19.
        //
        // This fits into a u64 whenever 51 + b + lg(19) < 64.
        //
        // Since 51 + b + lg(19) < 51 + 4.25 + b
        //                       = 55.25 + b,
        // this fits if b < 8.75.
        =#
        a3_19 = 19 * a[3+1]
        a4_19 = 19 * a[4+1]

        #=
        // Multiply to get 128-bit coefficients of output.
        //
        // The 128-bit multiplications by 2 turn into 1 slr + 1 slrd each,
        // which doesn't seem any better or worse than doing them as precomputations
        // on the 64-bit inputs.
        =#
        c0 = widemul(a[0+1], a[0+1]) + 2 * (widemul(a[1+1], a4_19) + widemul(a[2+1], a3_19))
        c1 = widemul(a[3+1], a3_19) + 2 * (widemul(a[0+1], a[1+1]) + widemul(a[2+1], a4_19))
        c2 = widemul(a[1+1], a[1+1]) + 2 * (widemul(a[0+1], a[2+1]) + widemul(a[4+1], a3_19))
        c3 = widemul(a[4+1], a4_19) + 2 * (widemul(a[0+1], a[3+1]) + widemul(a[1+1], a[2+1]))
        c4 = widemul(a[2+1], a[2+1]) + 2 * (widemul(a[0+1], a[4+1]) + widemul(a[1+1], a[3+1]))

        #=
        // Same bound as in multiply:
        //    c[i] < 2^(102 + 2*b) * (1+i + (4-i)*19)
        //         < 2^(102 + lg(1 + 4*19) + 2*b)
        //         < 2^(108.27 + 2*b)
        //
        // The carry (c[i] >> 51) fits into a u64 when
        //    108.27 + 2*b - 51 < 64
        //    2*b < 6.73
        //    b < 3.365.
        //
        // So we require b < 3 to ensure this fits.
        debug_assert!(a[0] < (1 << 54));
        debug_assert!(a[1] < (1 << 54));
        debug_assert!(a[2] < (1 << 54));
        debug_assert!(a[3] < (1 << 54));
        debug_assert!(a[4] < (1 << 54));
        =#

        #// Casting to u64 and back tells the compiler that the carry is bounded by 2^64, so
        #// that the addition is a u128 + u64 rather than u128 + u128.
        c1 += (c0 >> 51) % UInt64
        a0 = (c0 % UInt64) & LOW_51_BIT_MASK

        c2 += (c1 >> 51) % UInt64
        a1 = (c1 % UInt64) & LOW_51_BIT_MASK

        c3 += (c2 >> 51) % UInt64
        a2 = (c2 % UInt64) & LOW_51_BIT_MASK

        c4 += (c3 >> 51) % UInt64
        a3 = (c3 % UInt64)

        carry = (c4 >> 51) % UInt64
        a4 = (c4 % UInt64) & LOW_51_BIT_MASK

        #=
        // To see that this does not overflow, we need a[0] + carry * 19 < 2^64.
        //
        // c4 < a2^2 + 2*a0*a4 + 2*a1*a3 + (carry from c3)
        //    < 2^(102 + 2*b + lg(5)) + 2^64.
        //
        // When b < 3 we get
        //
        // c4 < 2^110.33  so that carry < 2^59.33
        //
        // so that
        //
        // a[0] + carry * 19 < 2^51 + 19 * 2^59.33 < 2^63.58
        //
        // and there is no overflow.
        =#
        a0 = a0 + carry * 19

        #// Now a[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
        a1 += a0 >> 51
        a0 &= LOW_51_BIT_MASK

        #// Now all a[i] < 2^(51 + epsilon) and a = self^(2^k).
    end

    FieldElement51((a0, a1, a2, a3, a4))
end


# Compute (self^(2^250-1), self^11), used as a helper function
# within invert() and pow22523().
function pow22501(x::T) where T
    #=
    // Instead of managing which temporary variables are used
    // for what, we define as many as we need and leave stack
    // allocation to the compiler
    //
    // Each temporary variable t_i is of the form (self)^e_i.
    // Squaring t_i corresponds to multiplying e_i by 2,
    // so the pow2k function shifts e_i left by k places.
    // Multiplying t_i and t_j corresponds to adding e_i + e_j.
    //
    // Temporary t_i                      Nonzero bits of e_i
    //
    =#
    t0  = square(x)           # 1         e_0 = 2^1
    t1  = square(square(t0))    # 3         e_1 = 2^3
    t2  = x * t1              # 3,0       e_2 = 2^3 + 2^0
    t3  = t0 * t2               # 3,1,0
    t4  = square(t3)             # 4,2,1
    t5  = t2 * t4               # 4,3,2,1,0
    t6  = pow2k(t5, 5)             # 9,8,7,6,5
    t7  = t6 * t5               # 9,8,7,6,5,4,3,2,1,0
    t8  = pow2k(t7, 10)            # 19..10
    t9  = t8 * t7               # 19..0
    t10 = pow2k(t9, 20)            # 39..20
    t11 = t10 * t9              # 39..0
    t12 = pow2k(t11, 10)           # 49..10
    t13 = t12 * t7              # 49..0
    t14 = pow2k(t13, 50)           # 99..50
    t15 = t14 * t13             # 99..0
    t16 = pow2k(t15, 100)          # 199..100
    t17 = t16 * t15             # 199..0
    t18 = pow2k(t17, 50)           # 249..50
    t19 = t18 * t13             # 249..0

    (t19, t3)
end


#=
/// Given a nonzero field element, compute its inverse.
///
/// The inverse is computed as self^(p-2), since
/// x^(p-2)x = x^(p-1) = 1 (mod p).
///
/// This function returns zero on input zero.
=#
function Base.inv(x::FieldElement51)
    # The bits of p-2 = 2^255 -19 -2 are 11010111111...11.
    (t19, t3) = pow22501(x)   # t19: 249..0 ; t3: 3,1,0
    t20 = pow2k(t19, 5)         # 254..5
    t21 = t20 * t3              # 254..5,3,1,0

    t21
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


