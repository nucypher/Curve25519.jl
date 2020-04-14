#=
/// Determine if this `FieldElement` is negative, in the sense
/// used in the ed25519 paper: `x` is negative if the low bit is
/// set.
///
/// # Return
///
/// If negative, return `Choice(1)`.  Otherwise, return `Choice(0)`.
=#
function is_negative(x)
    isodd(x)
end


#/// Returns 2 times the square of this field element.
square2(x) = let s = square(x)
    s + s
end


# Given `k > 0`, return `self^(2^k)`.
function pow2k(x::T, k::Int) where T

    x ^ (one(BigInt) << k)

    # TODO: the Dalek implementation is optimized for the specific representation they're using

    #=
    @assert k > 0

    a = x

    while true
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
        let a3_19 = 19 * a[3];
        let a4_19 = 19 * a[4];

        #=
        // Multiply to get 128-bit coefficients of output.
        //
        // The 128-bit multiplications by 2 turn into 1 slr + 1 slrd each,
        // which doesn't seem any better or worse than doing them as precomputations
        // on the 64-bit inputs.
        =#
        let     c0: u128 = m(a[0],  a[0]) + 2*( m(a[1], a4_19) + m(a[2], a3_19) );
        let mut c1: u128 = m(a[3], a3_19) + 2*( m(a[0],  a[1]) + m(a[2], a4_19) );
        let mut c2: u128 = m(a[1],  a[1]) + 2*( m(a[0],  a[2]) + m(a[4], a3_19) );
        let mut c3: u128 = m(a[4], a4_19) + 2*( m(a[0],  a[3]) + m(a[1],  a[2]) );
        let mut c4: u128 = m(a[2],  a[2]) + 2*( m(a[0],  a[4]) + m(a[1],  a[3]) );

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
        =#
        debug_assert!(a[0] < (1 << 54));
        debug_assert!(a[1] < (1 << 54));
        debug_assert!(a[2] < (1 << 54));
        debug_assert!(a[3] < (1 << 54));
        debug_assert!(a[4] < (1 << 54));

        const LOW_51_BIT_MASK: u64 = (1u64 << 51) - 1;

        # Casting to u64 and back tells the compiler that the carry is bounded by 2^64, so
        # that the addition is a u128 + u64 rather than u128 + u128.
        c1 += ((c0 >> 51) as u64) as u128;
        a[0] = (c0 as u64) & LOW_51_BIT_MASK;

        c2 += ((c1 >> 51) as u64) as u128;
        a[1] = (c1 as u64) & LOW_51_BIT_MASK;

        c3 += ((c2 >> 51) as u64) as u128;
        a[2] = (c2 as u64) & LOW_51_BIT_MASK;

        c4 += ((c3 >> 51) as u64) as u128;
        a[3] = (c3 as u64) & LOW_51_BIT_MASK;

        let carry: u64 = (c4 >> 51) as u64;
        a[4] = (c4 as u64) & LOW_51_BIT_MASK;

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
        a[0] = a[0] + carry * 19;

        // Now a[1] < 2^51 + 2^(64 -51) = 2^51 + 2^13 < 2^(51 + epsilon).
        a[1] += a[0] >> 51;
        a[0] &= LOW_51_BIT_MASK;

        // Now all a[i] < 2^(51 + epsilon) and a = self^(2^k).

        k = k - 1;
        if k == 0 {
            break;
        }
    }

    FieldElement51(a)
    =#
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
function invert(x::T) where T
    # The bits of p-2 = 2^255 -19 -2 are 11010111111...11.
    (t19, t3) = pow22501(x)   # t19: 249..0 ; t3: 3,1,0
    t20 = pow2k(t19, 5)         # 254..5
    t21 = t20 * t3              # 254..5,3,1,0

    t21
end


# Raise this field element to the power (p-5)/8 = 2^252 -3.
function pow_p58(x::T) where T
    # The bits of (p-5)/8 are 101111.....11.
    (t19, _) = pow22501(x)    # 249..0
    t20 = pow2k(t19, 2)       # 251..2
    t21 = x * t20             # 251..2,0

    t21
end

#=
/// Given `FieldElements` `u` and `v`, compute either `sqrt(u/v)`
/// or `sqrt(i*u/v)` in constant time.
///
/// This function always returns the nonnegative square root.
///
/// # Return
///
/// - `(Choice(1), +sqrt(u/v))  ` if `v` is nonzero and `u/v` is square;
/// - `(Choice(1), zero)        ` if `u` is zero;
/// - `(Choice(0), zero)        ` if `v` is zero and `u` is nonzero;
/// - `(Choice(0), +sqrt(i*u/v))` if `u/v` is nonsquare (so `i*u/v` is square).
///
=#
function sqrt_ratio_i(u::T, v::T) where T
    #=
    // Using the same trick as in ed25519 decoding, we merge the
    // inversion, the square root, and the square test as follows.
    //
    // To compute sqrt(α), we can compute β = α^((p+3)/8).
    // Then β^2 = ±α, so multiplying β by sqrt(-1) if necessary
    // gives sqrt(α).
    //
    // To compute 1/sqrt(α), we observe that
    //    1/β = α^(p-1 - (p+3)/8) = α^((7p-11)/8)
    //                            = α^3 * (α^7)^((p-5)/8).
    //
    // We can therefore compute sqrt(u/v) = sqrt(u)/sqrt(v)
    // by first computing
    //    r = u^((p+3)/8) v^(p-1-(p+3)/8)
    //      = u u^((p-5)/8) v^3 (v^7)^((p-5)/8)
    //      = (uv^3) (uv^7)^((p-5)/8).
    //
    // If v is nonzero and u/v is square, then r^2 = ±u/v,
    //                                     so vr^2 = ±u.
    // If vr^2 =  u, then sqrt(u/v) = r.
    // If vr^2 = -u, then sqrt(u/v) = r*sqrt(-1).
    //
    // If v is zero, r is also zero.
    =#

    v3 = square(v) * v
    v7 = square(v3) * v
    r = (u * v3) * pow_p58(u * v7)
    check = v * square(r)

    i = SQRT_M1

    correct_sign_sqrt = ct_eq(check, u)
    flipped_sign_sqrt   = ct_eq(check, -u)
    flipped_sign_sqrt_i = ct_eq(check, (-u)*i)

    r_prime = SQRT_M1 * r
    r = conditional_assign(r, r_prime, flipped_sign_sqrt | flipped_sign_sqrt_i)

    # Choose the nonnegative square root.
    r_is_negative = is_negative(r)
    r = conditional_negate(r, r_is_negative)

    was_nonzero_square = correct_sign_sqrt | flipped_sign_sqrt

    (was_nonzero_square, r)
end
