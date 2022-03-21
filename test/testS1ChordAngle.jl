# Original Java version Copyright 2014 Google Inc.
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

using S2
using Test
using LinearAlgebra

for i in 1:100
    frame = qr(randn(3, 3)).Q
    x, y, z = [frame[:,i] for i in 1:3]

    @test S2.S1Angle_ZERO == toAngle(S1ChordAngle(z, z))
    @test isapprox(π, radians(toAngle(S1ChordAngle(-z, z))), atol=1e-7)
    @test π / 2 ≈ radians(toAngle(S1ChordAngle(x, z)))
    w = normalize(y + z)
    @test π / 4 ≈ radians(toAngle(S1ChordAngle(w, z)))
end

@test 0.0 == degrees(toAngle(S1ChordAngle(0)))
@test 60.0 ≈ degrees(toAngle(S1ChordAngle(1)))
@test 90.0 ≈ degrees(toAngle(S1ChordAngle(2)))
@test 180.0 ≈ degrees(toAngle(S1ChordAngle(4)))
@test 180.0 ≈ degrees(toAngle(S1ChordAngle(5)))


@test S2.S1Angle_ZERO == toAngle(S2.S1ChordAngle_ZERO)


@test S2.isZero(S2.S1ChordAngle_ZERO)
@test !S2.isNegative(S2.S1ChordAngle_ZERO)
@test !S2.isSpecial(S2.S1ChordAngle_ZERO)
@test !S2.isSpecial(S2.S1ChordAngle_STRAIGHT)
@test S2.S2.isNegative(S2.S1ChordAngle_NEGATIVE)
@test S2.isSpecial(S2.S1ChordAngle_NEGATIVE)
@test S2.isInfinity(S2.S1ChordAngle_INFINITY)
@test S2.isSpecial(S2.S1ChordAngle_INFINITY)


@test 0.0 == radians(toAngle(S1ChordAngle(S2.S1Angle_ZERO)))
@test 4.0 == S1ChordAngle(radians(π)).distance2
@test isapprox(π, radians(toAngle(S1ChordAngle(radians(π)))), atol=1e-7)
@test radians(toAngle(S1ChordAngle(radians(-1)))) < 0.0
@test 1.0 ≈ radians(toAngle(S1ChordAngle(radians(1.0))))


let zero = S1ChordAngle(0)
    degree30 = S1ChordAngle(degrees(30))
    degree60 = S1ChordAngle(degrees(60))
    degree90 = S1ChordAngle(degrees(90))
    degree120 = S1ChordAngle(degrees(120))
    degree180 = S1ChordAngle(degrees(180))
    @test 0.0 == degrees(toAngle(zero + zero))
    @test 0.0 == degrees(toAngle(zero - zero))
    @test 0.0 == degrees(toAngle(degree60 - degree60))
    @test 0.0 == degrees(toAngle(degree180 - degree180))
    @test 0.0 == degrees(toAngle(zero - degree60))
    @test 0.0 == degrees(toAngle(degree30 - degree90))
    @test 60.0 ≈ degrees(toAngle(degree60 + zero))
    @test 60.0 ≈ degrees(toAngle(degree60 - zero))
    @test 60.0 ≈ degrees(toAngle(zero + degree60))
    @test 90.0 ≈ degrees(toAngle(degree30 + degree60))
    @test 90.0 ≈ degrees(toAngle(degree60 + degree30))
    @test 60.0 ≈ degrees(toAngle(degree90 - degree30))
    @test 30.0 ≈ degrees(toAngle(degree90 - degree60))
    @test 180.0 == degrees(toAngle(degree180 + zero))
    @test 180.0 == degrees(toAngle(degree180 - zero))
    @test 180.0 == degrees(toAngle(degree90 + degree90))
    @test 180.0 == degrees(toAngle(degree120 + degree90))
    @test 180.0 == degrees(toAngle(degree120 + degree120))
    @test 180.0 == degrees(toAngle(degree30 + degree180))
    @test 180.0 == degrees(toAngle(degree180 + degree180))
end

let N = 20
    for i in 0:N
        θ = π * i / N
        angle = S1ChordAngle(radians(θ))
        @test isapprox(sin(θ), sin(angle), atol=1e-15)
        @test isapprox(cos(θ), cos(angle), atol=1e-15)
        # Since the tan(x) is unbounded near Pi/4, we map the result back to an
        # angle before comparing.  (The assertion is that the result is equal to
        # the tangent of a nearby angle.)
        @test isapprox(atan(tan(θ)), atan(tan(angle)), atol=1e-15)
    end
end

# Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
let
    angle90 = S1ChordAngle(2)
    angle180 = S1ChordAngle(4)
    begin
        @test 1.0 == sin(angle90)
        @test 0.0 == cos(angle90)
        @test Inf == tan(angle90)
        @test 0.0 == sin(angle180)
        @test -1.0 == cos(angle180)
        @test 0.0 == tan(angle180)
    end
end


#  /**
#   * Verifies that the error bound returned by {@link S1ChordAngle#getS2PointConstructorMaxError} is
#   * large enough.
#   */
#  public void testGetS2PointConstructorMaxError() {
#    for (int iter = 0; iter < 10000; ++iter) {
#      rand.setSeed(iter);
#      S2Point x = randomPoint();
#      S2Point y = randomPoint();
#      if (super.oneIn(10)) {
#        // Occasionally test a point pair that is nearly identical or antipodal.
#        S1Angle r = S1Angle.radians(1e-15 * rand.nextDouble());
#        y = S2EdgeUtil.interpolateAtDistance(r, x, y);
#        if (oneIn(2)) {
#          y = S2Point.neg(y);
#        }
#      }
#      S1ChordAngle dist = new S1ChordAngle(x, y);
#      double error = dist.getS2PointConstructorMaxError();
#      String msg = "angle=" + dist + ", iter=" + iter;
#      assertTrue(msg, S2Predicates.compareDistance(x, y, dist.plusError(error).getLength2()) <= 0);
#      assertTrue(msg, S2Predicates.compareDistance(x, y, dist.plusError(-error).getLength2()) >= 0);
#    }
#  }
#
#  public void testS1AngleConsistency() {
#    // This test checks that the error bounds in the S1ChordAngle constructors
#    // are consistent with the maximum error in S1Angle(x, y).
#    double maxS1AngleError = 3.25 * DBL_EPSILON;
#    for (int iter = 0; iter < 10000; ++iter) {
#      S2Point x = randomPoint();
#      S2Point y = randomPoint();
#      S1ChordAngle dist1 = S1ChordAngle.fromS1Angle(new S1Angle(x, y));
#      S1ChordAngle dist2 = new S1ChordAngle(x, y);
#      double maxError =
#          (maxS1AngleError
#              + dist1.getS1AngleConstructorMaxError()
#              + dist2.getS2PointConstructorMaxError());
#      assertTrue(dist1.compareTo(dist2.plusError(maxError)) <= 0);
#      assertTrue(dist1.compareTo(dist2.plusError(-maxError)) >= 0);
#    }
#  }
#
#  /**
#   * Assert that {@code actual} is almost equal to {@code expected}, within floating point error.
#   */
#  private static void assertNearlyEquals(double expected, double actual) {
#    assertEquals(expected, actual, DOUBLE_ERROR);
#  }
#
