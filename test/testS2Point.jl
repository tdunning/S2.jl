# Original Java version Copyright 2013 Google Inc.
# This Julia translation Copyright 2022 Ted Dunning
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

using Test
using S2

import LinearAlgebra:(⋅), (×)

a = S2Point(1, 2, 3);
b = S2Point(0, 1, 0);

@test S2Point(0, 0) ≈ S2Point("0:0")
@test S2Point("0:0") ≈ S2.X_POS
@test S2Point(90, 0) == S2Point("90:0")
@test S2Point("90:0") ≈ S2.Z_POS
@test S2Point(0, 90) ≈ S2Point("0:90")
@test S2Point("0:90") ≈ S2.Y_POS

@test S2Point(1, 3, 3) == (a + b)
@test S2Point(1, 1, 3) == (a - b)
@test S2Point(2, 4, 6) == (2 * a)
@test S2Point(0.5, 1, 1.5) == (a / 2)
@test 2.0 == (a ⋅ b)
@test S2Point(-3, 0, 1) == a × b
@test b == normalize(b)
@test 14 == norm2(a)
@test sqrt(14) == norm(a)
@test 1 == norm2(normalize(a))
@test a == abs.(-a)
@test S2.Y_NEG == ortho(S2.X_POS)

begin
    # Check simple axial cases.
    p = S2Point("0:0");
    @test S2Point("0:90") ≈ rotate(p, S2.Z_POS, deg2rad(90)) 
    @test S2Point("-90:0") ≈ rotate(p, S2.Y_POS, deg2rad(90))
    @test S2Point("0:0") ≈ rotate(p, S2.X_POS, deg2rad(90)) 

    # Verify rotations in the plane containing the origin and two random points.
    for i = 1:1000
        let a = normalize(randn(3))
            b = normalize(randn(3))
            begin
                axis = a × b
                θ = angle(a, b);
                
                @test b ≈ rotate(a, axis, θ)
                # Rotate 'b' onto 'a'.
                @test a ≈ rotate(b, axis, -θ)
                # Rotate 'a' to the midpoint of 'a' and 'b'.
                @test isapprox(normalize(a + b), rotate(a, axis, θ / 2), atol=1e-12)
                # Rotate 'b' to the antipodal point of 'a'.
                @test isapprox((-a), rotate(b, axis, π - θ), atol=1e-12)
            end
        end
    end
end

begin
    point = normalize(randn(3))
    expectedCapBound = S2Cap(point, S1ChordAngle(0))
    @test expectedCapBound == bound(S2Cap, point)
    
    ll = S2LatLng(point)
    @test fromPointPair(S2LatLngRect, ll, ll) == bound(S2LatLngRect, point)
    @test fromPoint(S2LatLngRect, ll) == bound(S2LatLngRect, point)
    
    # The leaf cell containing a point is still much larger than the point.
    cell = S2Cell(point)
    @test !(contains(point, cell))
    @test mayIntersect(point, cell)
end

# let points = [S2Point("0:0"), S2Point("1:1"), S2Point("2:2")]
#     for j = 0:(length(points)-1]
#         subset = points[0:(j-1)]
#         S2Shape shape = S2Point.Shape.fromList(subset);
#         @test !(containsOrigin(shape);
#         @test !(hasInterior(shape));
#         @test length(subset) == numEdges(shape)
#     end
# end

