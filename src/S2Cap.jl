md"""
An `S2Cap` is a circular patch on a sphere defined by a center and a angular radius.

The radius is expressed as a chord angle to make computations easier.
"""
struct S2Cap <: S2Region
    center::Vector
    radius::S1ChordAngle
    S2Cap(center::Vector, edge::Vector) = verifyPoint(center) && new(center, S1ChordAngle(center, edge))
    S2Cap(center::Vector, radius::S1ChordAngle) = verifyPoint(center) && new(center, radius)
    S2Cap(center::Vector, angle::S1Angle) = verifyPoint(center) && new(center, S1ChordAngle(angle))
end

equals(a::S2Cap, b::S2Cap) = a.center == b.center && a.radius == b.radius


    

