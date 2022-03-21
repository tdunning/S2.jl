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

# This class represents a point on the unit sphere as a pair of latitude-longitude coordinates.

struct S2LatLng
    latRadians::Float64
    lngRadians::Float64
end

S2LatLng(lat::S1Angle, lon::S1Angle) = S2LatLng(radians(lat), radians(lon))
S2LatLng(p::Vector) = verifyPoint(p) && S2LatLng(latitude(p), longitude(p))

radians2LL(lat::Number, lon::Number) = S2LatLng(lat, lon)
deg2LL(lat::Number, lon::Number) = S2LatLng(deg2rad(lat), deg2rad(lon))

latitude(p::Vector) = verifyPoint(p) && S1Angle(atan(p[3], hypot(p[1], p[2])))
latitude(p::S2LatLng) = S1Angle(p.latRadians)
longitude(p::Vector) = verifyPoint(p) && S1Angle(atan(p[2], p[1]))
longitude(p::S2LatLng) = S1Angle(p.lngRadians)

md"""
Returns a new S2LatLng based on this instance for which {@link #isValid()} will be {@code
true}.

* Latitude is clipped to the range {@code [-90, 90]}
* Longitude is normalized to be in the range {@code [-180, 180]}

f the current point is valid then the returned point will have the same coordinates.
"""
function LinearAlgebra.normalize(p::S2LatLng)::S2LatLng
    S2LatLng(max(-π, min(π, p.latRadians)),
             normalizeRange(p.lngRadians))
end

function S2Point(ll::S2LatLng) 
    ϕ = ll.latRadians
    θ = lngRadians
    cosphi = cos(ϕ)
    return S2Point(cos(θ) * cosphi, sin(θ) * cosphi, sin(ϕ))
end

md"""
Returns the distance (measured along the surface of the sphere) to the given point.

This code uses the Haversine formula which is very good for short distances but is
not as good for nearly antipodal points (but that only means errors of 10cm). An
alternative is to convert to S2Points and then to S1Angle to get 15 digits of accuracy
for all distances.  
"""
function distance(a::S2LatLng, b::S2LatLng)
    dlat = Math.sin(0.5 * (a.latRadians - b.latRadians))
    dlng = Math.sin(0.5 * (a.lngRadians - b.lngRadians))
    x = dlat * dlat + dlng * dlng * cos(a.latRadians) * cos(b.latRadians)
    return S1Angle.radians(2 * asin(sqrt(min(1.0, x))))
end

md"Returns the surface distance to the given point assuming a constant radius."
distance(a::S2LatLng, b::S2LatLng, radius) = distance(distance(a, b), radius)

md"""
Adds the given point to this point. Note that there is no guarantee that the new point will be
valid.
"""
(+)(a::S2LatLng, b::S2LatLng) = S2LatLng(a.latRadians + b.latRadians, a.lngRadians + o.lngRadians)

md"""
Subtracts the given point from this point. Note that there is no guarantee that the new point
will be <em>valid</em>.
"""
(-)(a::S2LatLng, b::S2LatLng) = S2LatLng(a.latRadians - b.latRadians, a.lngRadians - b.lngRadians)

(*)(a::S2LatLng, m) = S2LatLng(a.latRadians * m, a.lngRadians * m)

equals(a::S2LatLng, b::S2LatLng) = (a.latRadians == b.latRadians) && (a.lngRadians == b.lngRadians)

