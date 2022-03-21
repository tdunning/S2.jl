module S2
abstract type S2Region end

using Markdown, LinearAlgebra

import Base:(*), (+), (-), (/), (<), (≤)
import Base:min, max

S2_M_PI = Float64(π)

include("S2Point.jl")
include("S1Interval.jl")
#include("S2Cell.jl")
include("S1Angle.jl")
include("S1ChordAngle.jl")
include("S2Cap.jl")
include("S2LatLng.jl")

export S1Angle, radians, degrees, distance, distance2, normalizeRange
export S2Point, norm2, ortho, rotate, full, isfull, isinverted
export S1ChordAngle, toAngle
export S1Interval, center, complement, complementcenter
export interiorcontains, intersects, interiorintersects, union, intersection
export addPoint, clampPoint, fromPointPair, fromPoint, expanded, getDirectedHausdorffDistance
export S2Cap, getCapBound
export S2LatLng

end
