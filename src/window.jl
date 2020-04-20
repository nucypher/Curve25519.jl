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
struct LookupTable{T}
    vals :: Array{T, 1}

    function LookupTable(p::EdwardsPoint{T}) where T
        points = Array{ProjectiveNielsPoint{T}}(undef, 8)
        points[1] = to_projective_niels(p)
        for j in 0:6
            points[j+1+1] = to_projective_niels(to_extended(p + points[j+1]))
        end
        new{ProjectiveNielsPoint{T}}(points)
    end
end


# Given \\(-8 \leq x \leq 8\\), return \\(xP\\) in constant time.
function Base.getindex(table::LookupTable, x::Union{CT.Value, Integer})
    @assert CT.unwrap(x) >= -8 && CT.unwrap(x) <= 8
    # Compute xabs = |x|
    xmask = x >> 7
    xabs = xor(x + xmask, xmask)

    t = table.vals[xabs]
    CT.select(isodd(xmask), -t, t)
end
