#=
/// A lookup table of precomputed multiples of a point \\(P\\), used to
/// compute \\( xP \\) for \\( -8 \leq x \leq 8 \\).
///
/// The computation of \\( xP \\) is done in constant time by the `select` function.
///
/// Since `LookupTable` does not implement `Index`, it's more difficult
/// to accidentally use the table directly.  Unfortunately the table is
/// only `pub(crate)` so that we can write hardcoded constants, so it's
/// still technically possible.  It would be nice to prevent direct
/// access to the table.
///
/// XXX make this generic with respect to table size
#[derive(Copy, Clone)]
=#
struct LookupTable{P}
    vals :: Array{P, 1}

    function LookupTable{ProjectiveNielsPoint{T}}(p::EdwardsPoint{T}) where T
        points = Array{ProjectiveNielsPoint{T}}(undef, 8)
        points[1] = to_projective_niels(p)
        for j in 0:6
            points[j+1+1] = to_projective_niels(to_extended(p + points[j+1]))
        end
        new{ProjectiveNielsPoint{T}}(points)
    end

    function LookupTable{AffineNielsPoint{T}}(p::EdwardsPoint{T}) where T
        points = Array{AffineNielsPoint{T}}(undef, 8)
        points[1] = to_affine_niels(p)
        for j in 0:6
            points[j+1+1] = to_affine_niels(to_extended(p + points[j+1]))
        end
        new{AffineNielsPoint{T}}(points)
    end

end


# Given \\(-8 \leq x \leq 8\\), return \\(xP\\) in constant time.
@inline function Base.getindex(table::LookupTable{P}, x::Union{CT.Value, Integer}) where P
    @assert CT.unwrap(x) >= -8 && CT.unwrap(x) <= 8
    # Compute xabs = |x|
    xmask = x >> 7
    xabs = xor(x + xmask, xmask)

    t = get(table.vals, xabs, zero(P))
    CT.select(isodd(xmask), -t, t)
end


# From edwards.rs

struct EdwardsBasepointTable{T}

    tables :: Array{LookupTable{AffineNielsPoint{T}}, 1}

    function EdwardsBasepointTable(basepoint::EdwardsPoint{T}) where T
        tables = Array{LookupTable{AffineNielsPoint{T}}}(undef, 32)
        P = basepoint
        for i in 0:31
            # P = (16^2)^i * B
            tables[i+1] = LookupTable{AffineNielsPoint{T}}(P)
            P = mul_by_pow_2(P, 8)
        end
        new{T}(tables)
    end
end


@Base.propagate_inbounds @inline function Base.:*(table::EdwardsBasepointTable{T}, scalar::Z) where {T, Z}
    a = to_radix_16(scalar)
    tables = table.tables

    P = zero(EdwardsPoint{T})
    for i in 1:2:63
        P = to_extended(P + tables[(i>>1)+1][a[i+1]+1])
    end

    P = mul_by_pow_2(P, 4)

    for i in 0:2:63
        P = to_extended(P + tables[(i>>1)+1][a[i+1]+1])
    end

    P
end


struct RistrettoBasepointTable{T}

    table :: EdwardsBasepointTable{T}

    function RistrettoBasepointTable(basepoint::RistrettoPoint{T}) where T
        new{T}(EdwardsBasepointTable(basepoint.ep))
    end
end


@inline function Base.:*(table::RistrettoBasepointTable{T}, scalar::Z) where {T, Z}
    RistrettoPoint{T}(table.table * scalar)
end
