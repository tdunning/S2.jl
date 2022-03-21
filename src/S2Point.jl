#import Base:isapprox, abs, angle

# isapprox, norm, norm2, normalize, abs, ortho, rotate, angle

S2Point(x, y, z) = [x, y, z]
S2Point(lat::Number, lon::Number) = let
    x = cos(deg2rad(lat)) * cos(deg2rad(lon))
    y = cos(deg2rad(lat)) * sin(deg2rad(lon))
    z = sin(deg2rad(lat))
    S2Point(x, y, z)
end
S2Point(latlon::String) = let args = parse.(Float64, split(latlon, ':'))
    S2Point(args[1], args[2])
end

# Origin of the coordinate system, [0,0,0]
ORIGIN = S2Point(0, 0, 0)

# Direction of the x-axis.
X_POS = S2Point(1, 0, 0)

# Opposite direction of the x-axis.
X_NEG = S2Point(-1, 0, 0)

# Direction of the y-axis.
Y_POS = S2Point(0, 1, 0)

# Opposite direction of the y-axis.
Y_NEG = S2Point(0, -1, 0)

# Direction of the z-axis. 
Z_POS = S2Point(0, 0, 1)

# Opposite direction of the z-axis.
Z_NEG = S2Point(0, 0, -1)

md"""
Verify that a vector is a valid S2 point
"""
verifyPoint(a::Vector) = (length(a) == 3 && eltype(a) <: Number) || throw(ArgumentError("Value is not a valid S2Point: $a"))

coordinates(a::Vector) = verifyPoint(a) && (a[1], a[2], a[3])

getCapBound(a::Vector) = verifyPoint(a) && S2Cap(a, S1ChordAngle(0))

norm2(a::Vector) = sum(a.^2)

function distance2(a::Vector, b::Vector)
    verifyPoint(a)
    verifyPoint(b)
    dx = a - b
    return dx ⋅ dx
end

distance(a::Vector, b::Vector) = sqrt(distance2(a, b))

md"""
The scalar triple product is `a ⋅ (b × c)` but more efficiently computed.
"""
function scalarTripleProduct(a::Vector, b::Vector)
    x = b[2] * c[3] - b[3] * c[2]
    y = b[3] * c[1] - b[1] * c[3]
    z = b[1] * c[2] - b[2] * c[1]
    return a[1] * x + a[2] * y + a[3] * z
end

function ortho(a::Vector)
    k = largestAbsComponent(a)
    if k == 1
        normalize(a × Z_POS)
    elseif k == 2
        normalize(a × X_POS)
    elseif k == 3
        normalize(a × Y_POS)
    else
        throw(ArgumentError("Maximum component must be 1, 2, or 3"))
    end
end

largestAbsComponent(a::Vector) = findmax(coordinates(a))[2]



md"""
Rotates this point around an arbitrary axis. The result is normalized.

@param axis point around which rotation should be performed.
@param radians radians to rotate the point counterclockwise around the given axis.
"""
function rotate(p::Vector, axis::Vector, radians::Real)::Vector{Float64}
    verifyPoint(p)
    verifyPoint(axis)
    
    p = normalize(p)
    axis = normalize(axis)
    pointOnAxis = (p ⋅ axis) * axis
    axisToPoint = p - pointOnAxis
    axisToPointNormal = axis × axisToPoint
    axisToPoint = cos(radians) * axisToPoint 
    axisToPointNormal = sin(radians) * axisToPointNormal
    return normalize(axisToPoint + axisToPointNormal + pointOnAxis)
end

md"""
The angle in radians between two points
"""
Base.angle(a::Vector, b::Vector) = atan(norm(a × b), a ⋅ b)

