struct RistrettoPoint{T}
    ep :: EdwardsPoint{T}
end


Base.zero(::Type{RistrettoPoint{T}}) where T = RistrettoPoint{T}(zero(EdwardsPoint{T}))


Base.:+(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T = RistrettoPoint{T}(p.ep + q.ep)


Base.:*(p::RistrettoPoint{T}, s::Z) where {T, Z} = RistrettoPoint{T}(p.ep * s)


function Base.:(==)(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T
    X1Y2 = p.ep.X * p.ep.Y
    Y1X2 = p.ep.Y * p.ep.X
    X1X2 = p.ep.X * p.ep.X
    Y1Y2 = p.ep.Y * p.ep.Y

    ct_eq(X1Y2, Y1X2) | ct_eq(X1X2, Y1Y2)
end


const RistrettoPointT = RistrettoPoint{InternalScalar}
