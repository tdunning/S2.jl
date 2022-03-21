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

using S2
using LinearAlgebra
using Test


# Basics

# Check that the conversion between Pi radians and 180 degrees is exact.
@test Float64(π) == radians(S1Angle(π))
@test 180.0 == degrees(S1Angle(π))
@test Float64(π) == radians(degrees(180))
@test 180.0 == degrees(degrees(180))

@test 90.0 == degrees(radians(π / 2))
@test 90.0 == degrees(S1Angle(π / 2))

@test -90.0 == degrees(S1Angle(-π / 2))
@test Float64(-π / 4) == radians(degrees(-45))


@test 30.0 ≈ degrees(degrees(10) + degrees(20))
@test -10 ≈ degrees(degrees(10) - degrees(20))
@test 20 ≈ degrees(2 * degrees(10))
@test 5 ≈ degrees(degrees(10) / 2.0)
@test 1 ≈ cos(S1Angle(0))
@test 1 ≈ sin(degrees(90))
@test 1 ≈ tan(degrees(45))


@test 100.0 * π ≈ distance(S1Angle(π), 100.0)
@test 50.0 * π ≈ distance(S1Angle(π / 2), 100.0)
@test 25.0 * π ≈ distance(radians(π / 4), 100.0)


@test 0.0 ≈ radians(normalizeRange(degrees(360.0)))
@test -90.0 ≈ degrees(normalizeRange(degrees(-90.0)))
@test 180.0 ≈ degrees(normalizeRange(degrees(-180.0)))
@test 90.0 ≈ degrees(normalizeRange(degrees(90.0)))
@test 180.0 ≈ degrees(normalizeRange(degrees(180.0)))
@test -90.0 ≈ degrees(normalizeRange(degrees(270.0)))
@test 180.0 ≈ degrees(normalizeRange(degrees(540.0)))
@test 90.0 ≈ degrees(normalizeRange(degrees(-270.0)))

@test π ≈ radians(normalizeRange(radians(π)))

# -π maps to π.
@test π ≈ radians(normalizeRange(radians(-π)))

let a = degrees(90.0)
    @test a === normalizeRange(a)
end

@test 1.5 == radians(S1Angle(0.5) + S1Angle(0.75) + S1Angle(0.25))

@test 180 == degrees(normalizeRange(degrees(90) + degrees(45) + degrees(30) + degrees(15)))

@test 0.75 * π ≈ radians(S1Angle(π / 2) + degrees(45))

