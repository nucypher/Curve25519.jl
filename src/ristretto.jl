struct RistrettoCurve <: DarkCurves.EllipticCurve end


struct RistrettoPoint{T} <: DarkCurves.EllipticCurvePoint{RistrettoCurve, T}
    ep :: EdwardsPoint{T}

    RistrettoPoint(ep::EdwardsPoint{T}) where T = new{T}(ep)
    RistrettoPoint{T}(ep::EdwardsPoint{T}) where T = new{T}(ep)
end


CT.wrap(p::RistrettoPoint{T}) where T = RistrettoPoint(CT.wrap(p.ep))


Base.zero(::Type{RistrettoPoint{T}}) where T = RistrettoPoint{T}(zero(EdwardsPoint{T}))


Base.:+(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T = RistrettoPoint{T}(p.ep + q.ep)


Base.:-(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T = RistrettoPoint{T}(p.ep - q.ep)


Base.:-(p::RistrettoPoint{T}) where T = RistrettoPoint{T}(-p.ep)


Base.:*(p::RistrettoPoint{T}, s::Z) where {T, Z<:Union{MgModUInt, ModUInt}} = RistrettoPoint{T}(p.ep * s)


function Base.:(==)(p::RistrettoPoint{T}, q::RistrettoPoint{T}) where T
    X1Y2 = p.ep.X * p.ep.Y
    Y1X2 = p.ep.Y * p.ep.X
    X1X2 = p.ep.X * p.ep.X
    Y1Y2 = p.ep.Y * p.ep.Y

    (X1Y2 == Y1X2) | (X1X2 == Y1Y2)
end
