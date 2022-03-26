module S2
abstract type S2Region end

using Markdown, LinearAlgebra

import Base:(*), (+), (-), (/), (<), (≤)
import Base:min, max


S2_M_PI = Float64(π)

include("S1Angle.jl")
include("S1ChordAngle.jl")
include("S2Cap.jl")
include("S2LatLng.jl")
include("R1Interval.jl")
include("S1Interval.jl")
include("S2LatLngRect.jl")
include("S2Point.jl")
include("S2Cell.jl")

export S1Angle, radians, degrees, distance, distance2, normalizeRange
export S2Point, norm2, ortho, rotate, full, isfull, isinverted
export S1ChordAngle, toAngle
export R1Interval, S1Interval, center, complement, complementcenter
export interiorcontains, intersects, interiorintersects, union, intersection
export addPoint, clampPoint, fromPointPair, fromPoint, expanded, getDirectedHausdorffDistance
export S2Cap, getCapBound
export S2LatLng, S2LatLngRect, bound
export S2Cell

end
