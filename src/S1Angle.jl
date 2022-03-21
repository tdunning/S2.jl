# Original Java version Copyright 2005 Google Inc.
# Julia translation Copyright 2022 Ted Dunning
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

struct S1Angle
    radians::Float64

    S1Angle(radians::Number) = new(radians)

    md"""
    Return the angle between two points, which is also equal to the distance between these points
    on the unit sphere. The points do not need to be normalized.
    """
    S1Angle(x::Vector, y::Vector) = new(angle(x, y))
end

degrees(degrees::Number) = S1Angle(deg2rad(degrees))
radians(radians::Number) = S1Angle(radians)

md"""
Returns the angle in radians
"""
radians(x::S1Angle) = x.radians

md"""
Returns the angle in degrees
"""
degrees(x::S1Angle) = rad2deg(x.radians)

md"""
An angle larger than any finite angle
"""
const S1Angle_INFINITY = S1Angle(Inf)

const S1Angle_ZERO = S1Angle(0)

(<)(a::S1Angle, b::S1Angle) = a.radians < b.radians

(≤)(a::S1Angle, b::S1Angle) = a.radians ≤ b.radians

max(a::S1Angle, b::S1Angle) = a > b ? a : b;
min(a::S1Angle, b::S1Angle) = a < b ? a : b;

md"""
Returns the distance along the surface of a sphere of the given radius
"""
distance(a::S1Angle, radius::Number) = a.radians * radius;

(-)(a::S1Angle) = S1Angle(-a.radians)
(+)(a::S1Angle, b::S1Angle) = S1Angle(a.radians + b.radians)
(-)(a::S1Angle, b::S1Angle) = S1Angle(a.radians - b.radians)
(*)(k::Number, a::S1Angle) = S1Angle(k * a.radians)
(/)(a::S1Angle, k::Number) = S1Angle(a.radians / k)
Base.cos(a::S1Angle) = cos(a.radians)
Base.sin(a::S1Angle) = sin(a.radians)
Base.tan(a::S1Angle) = tan(a.radians)


md"""
Returns the angle normalized to the range (-180, 180] degrees.
"""
normalizeRange(a::S1Angle) = S1Angle(normalizeRange(a.radians))

function normalizeRange(r:: Number)
    if -π < r ≤ π
        return r
    else
        b = r - 2π * floor(r / (2π))
        if b ≤ π
            return b
        else
            return b - 2π
        end
    end
end


