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

begin
    # Create some standard intervals to use in the tests.  These include the
    # empty and full intervals, intervals containing a single point, and
    # intervals spanning one or more "quadrants" which are numbered as follows:
    #    quad1 == [0, π/2]
    #    quad2 == [π/2, π]
    #    quad3 == [-π, -π/2]
    #    quad4 == [-π/2, 0]
    emptyInterval = empty(S1Interval);
    fullInterval = full(S1Interval);

    # Single-point intervals:
    zero =  S1Interval(0, 0);
    pi2 =  S1Interval(π/2, π/2);
    pi = S1Interval(π, π);
    mipi = S1Interval(-π, -π) # Same as "pi" after normalization.
    mipi2 = S1Interval(-π/2, -π/2);

    # Single quadrants:
    quad1 = S1Interval(0, π/2);
    quad2 = S1Interval(π/2, -π);
    quad3 = S1Interval(π, -π/2);
    quad4 = S1Interval(-π/2, 0);
    # Quadrant pairs:
    quad12 = S1Interval(0, -π);
    quad23 = S1Interval(π/2, -π/2);
    quad34 = S1Interval(-π, 0);
    quad41 = S1Interval(-π/2, π/2);
    # Quadrant triples:
    quad123 = S1Interval(0, -π/2);
    quad234 = S1Interval(π/2, 0);
    quad341 = S1Interval(π, π/2);
    quad412 = S1Interval(-π/2, -π);
    # Small intervals around the midpoints between quadrants, such that
    # the center of each interval is offset slightly CCW from the midpoint.
    mid12 = S1Interval(π/2 - 0.01, π/2 + 0.02);
    mid23 = S1Interval(π - 0.01, -π + 0.02);
    mid34 = S1Interval(-π/2 - 0.01, -π/2 + 0.02);
    mid41 = S1Interval(-0.01, 0.02);

    # Spot-check the constructors and accessors.
    @test 0.0 == quad12.lo
    @test Float64(π) == quad12.hi
    @test Float64(π) == pi.lo
    @test Float64(π) == pi.hi

    # Check that [-π, -π] is normalized to [π, π].
    @test S2.S2_M_PI == mipi.lo
    @test S2.S2_M_PI == mipi.hi
    @test Float64(S2.S2_M_PI/2) == quad23.lo
    @test Float64(-S2.S2_M_PI/2) == quad23.hi

    # is_valid(), is_empty(), is_full(), is_inverted()
    @test isvalid(zero) && !isempty(zero) && !isfull(zero)
    @test isvalid(emptyInterval) && isempty(emptyInterval) && !isfull(emptyInterval)
    @test isinverted(emptyInterval)
    @test isvalid(fullInterval) && !isempty(fullInterval) && isfull(fullInterval)
    @test !isempty(quad12) && !isfull(quad12) && !isinverted(quad12)
    @test !isempty(quad23) && !isfull(quad23) && isinverted(quad23)
    @test isvalid(pi) && !isempty(pi) && !isinverted(pi)
    @test isvalid(mipi) && !isempty(mipi) && !isinverted(mipi)
    
    @test S2.S2_M_PI/2 == center(quad12)
    @test 3.0 - S2.S2_M_PI == center(S1Interval(3.1, 2.9))
    @test S2.S2_M_PI - 3.0 == center(S1Interval(-2.9, -3.1))
    @test S2.S2_M_PI == center(S1Interval(2.1, -2.1))
    @test S2.S2_M_PI == center(pi)
    @test S2.S2_M_PI == center(mipi)
    @test S2.S2_M_PI == abs(center(quad23))
    @test 0.75 * S2.S2_M_PI == center(quad123)
    
    
    @test S2.S2_M_PI == length(quad12)
    @test 0.0 == length(pi)
    @test 0.0 == length(mipi)
    @test 1.5 * S2.S2_M_PI == length(quad123)
    @test S2.S2_M_PI == abs(length(quad23))
    @test 2 * S2.S2_M_PI == length(fullInterval)
    @test length(emptyInterval) < 0
    
    
    
    @test isfull(complement(emptyInterval))
    @test isempty(complement(fullInterval))
    @test isfull(complement(pi))
    @test isfull(complement(mipi))
    @test isfull(complement(zero))
    @test complement(quad12) ≈ quad34
    @test complement(quad34) ≈ quad12
    @test complement(quad123) ≈ quad4
    
    
    @test !contains(emptyInterval, 0) && !contains(emptyInterval, S2.S2_M_PI) && !contains(emptyInterval, -S2.S2_M_PI)
    @test !interiorcontains(emptyInterval, S2.S2_M_PI) && !interiorcontains(emptyInterval, -S2.S2_M_PI)
    @test contains(fullInterval, 0) && contains(fullInterval, S2.S2_M_PI) && contains(fullInterval, -S2.S2_M_PI)
    @test interiorcontains(fullInterval, S2.S2_M_PI) && interiorcontains(fullInterval, -S2.S2_M_PI)
    @test contains(quad12, 0) && contains(quad12, S2.S2_M_PI) && contains(quad12, -S2.S2_M_PI)
    @test interiorcontains(quad12, S2.S2_M_PI/2) && !interiorcontains(quad12, 0)
    @test !interiorcontains(quad12, S2.S2_M_PI) && !interiorcontains(quad12, -S2.S2_M_PI)
    @test contains(quad23, S2.S2_M_PI/2) && contains(quad23, -S2.S2_M_PI/2)
    @test contains(quad23, S2.S2_M_PI) && contains(quad23, -S2.S2_M_PI)
    @test !contains(quad23, 0)
    @test !interiorcontains(quad23, S2.S2_M_PI/2) && !interiorcontains(quad23, -S2.S2_M_PI/2)
    @test interiorcontains(quad23, S2.S2_M_PI) && interiorcontains(quad23, -S2.S2_M_PI)
    @test !interiorcontains(quad23, 0)
    @test contains(pi, S2.S2_M_PI) && contains(pi, -S2.S2_M_PI) && !contains(pi, 0)
    @test !interiorcontains(pi, S2.S2_M_PI) && !interiorcontains(pi, -S2.S2_M_PI)
    @test contains(mipi, S2.S2_M_PI) && contains(mipi, -S2.S2_M_PI) && !contains(mipi, 0)
    @test !interiorcontains(mipi, S2.S2_M_PI) && !interiorcontains(mipi, -S2.S2_M_PI)
    @test contains(zero, 0) && !interiorcontains(zero, 0)


    function testIntervalOps(x::S1Interval, y::S1Interval,
                             expectedRelation::String,
                             expectedUnion::S1Interval, expectedIntersection::S1Interval) 
        # Test all of the interval operations on the given pair of intervals.
        # "expectedRelation" is a sequence of "T" and "F" characters corresponding
        # to the expected results of Contains(), InteriorContains(), Intersects(),
        # and InteriorIntersects() respectively.
        @test contains(x, y) == (expectedRelation[1] == 'T')
        @test interiorcontains(x, y) == (expectedRelation[2] == 'T')
        @test intersects(x, y) == (expectedRelation[3] == 'T')
        @test interiorintersects(x, y) == (expectedRelation[4] == 'T')

        # bounds() returns a const reference to a member variable, so we need to
        # make a copy when invoking it on a temporary object.
        @test expectedUnion == union(x, y)
        @test expectedIntersection == intersection(x, y)

        @test contains(x, y) == (union(x, y) == x)
        @test intersects(x, y) == !isempty(intersection(x, y))

        if y.lo == y.hi
            r = addPoint(x, y.lo);
            @test expectedUnion == r
        end
    end


    begin
        # Contains(S1Interval), InteriorContains(S1Interval),
        # Intersects(), InteriorIntersects(), Union(), Intersection()
        testIntervalOps(emptyInterval, emptyInterval, "TTFF", emptyInterval, emptyInterval);
        testIntervalOps(emptyInterval, fullInterval, "FFFF", fullInterval, emptyInterval);
        testIntervalOps(emptyInterval, zero, "FFFF", zero, emptyInterval);
        testIntervalOps(emptyInterval, pi, "FFFF", pi, emptyInterval);
        testIntervalOps(emptyInterval, mipi, "FFFF", mipi, emptyInterval);
        
        testIntervalOps(fullInterval, emptyInterval, "TTFF", fullInterval, emptyInterval);
        testIntervalOps(fullInterval, fullInterval, "TTTT", fullInterval, fullInterval);
        testIntervalOps(fullInterval, zero, "TTTT", fullInterval, zero);
        testIntervalOps(fullInterval, pi, "TTTT", fullInterval, pi);
        testIntervalOps(fullInterval, mipi, "TTTT", fullInterval, mipi);
        testIntervalOps(fullInterval, quad12, "TTTT", fullInterval, quad12);
        testIntervalOps(fullInterval, quad23, "TTTT", fullInterval, quad23);
        
        testIntervalOps(zero, emptyInterval, "TTFF", zero, emptyInterval);
        testIntervalOps(zero, fullInterval, "FFTF", fullInterval, zero);
        testIntervalOps(zero, zero, "TFTF", zero, zero);
        testIntervalOps(zero, pi, "FFFF", S1Interval(0, S2.S2_M_PI), emptyInterval);
        testIntervalOps(zero, pi2, "FFFF", quad1, emptyInterval);
        testIntervalOps(zero, mipi, "FFFF", quad12, emptyInterval);
        testIntervalOps(zero, mipi2, "FFFF", quad4, emptyInterval);
        testIntervalOps(zero, quad12, "FFTF", quad12, zero);
        testIntervalOps(zero, quad23, "FFFF", quad123, emptyInterval);
        
        testIntervalOps(pi2, emptyInterval, "TTFF", pi2, emptyInterval);
        testIntervalOps(pi2, fullInterval, "FFTF", fullInterval, pi2);
        testIntervalOps(pi2, zero, "FFFF", quad1, emptyInterval);
        testIntervalOps(pi2, pi, "FFFF", S1Interval(S2.S2_M_PI/2, S2.S2_M_PI), emptyInterval);
        testIntervalOps(pi2, pi2, "TFTF", pi2, pi2);
        testIntervalOps(pi2, mipi, "FFFF", quad2, emptyInterval);
        testIntervalOps(pi2, mipi2, "FFFF", quad23, emptyInterval);
        testIntervalOps(pi2, quad12, "FFTF", quad12, pi2);
        testIntervalOps(pi2, quad23, "FFTF", quad23, pi2);
        
        testIntervalOps(pi, emptyInterval, "TTFF", pi, emptyInterval);
        testIntervalOps(pi, fullInterval, "FFTF", fullInterval, pi);
        testIntervalOps(pi, zero, "FFFF", S1Interval(S2.S2_M_PI, 0), emptyInterval);
        testIntervalOps(pi, pi, "TFTF", pi, pi);
        testIntervalOps(pi, pi2, "FFFF", S1Interval(S2.S2_M_PI/2, S2.S2_M_PI), emptyInterval);
        testIntervalOps(pi, mipi, "TFTF", pi, pi);
        testIntervalOps(pi, mipi2, "FFFF", quad3, emptyInterval);
        testIntervalOps(pi, quad12, "FFTF", S1Interval(0, S2.S2_M_PI), pi);
        testIntervalOps(pi, quad23, "FFTF", quad23, pi);
        
        testIntervalOps(mipi, emptyInterval, "TTFF", mipi, emptyInterval);
        testIntervalOps(mipi, fullInterval, "FFTF", fullInterval, mipi);
        testIntervalOps(mipi, zero, "FFFF", quad34, emptyInterval);
        testIntervalOps(mipi, pi, "TFTF", mipi, mipi);
        testIntervalOps(mipi, pi2, "FFFF", quad2, emptyInterval);
        testIntervalOps(mipi, mipi, "TFTF", mipi, mipi);
        testIntervalOps(mipi, mipi2, "FFFF", S1Interval(-S2.S2_M_PI, -S2.S2_M_PI/2), emptyInterval);
        testIntervalOps(mipi, quad12, "FFTF", quad12, mipi);
        testIntervalOps(mipi, quad23, "FFTF", quad23, mipi);
        
        testIntervalOps(quad12, emptyInterval, "TTFF", quad12, emptyInterval);
        testIntervalOps(quad12, fullInterval, "FFTT", fullInterval, quad12);
        testIntervalOps(quad12, zero, "TFTF", quad12, zero);
        testIntervalOps(quad12, pi, "TFTF", quad12, pi);
        testIntervalOps(quad12, mipi, "TFTF", quad12, mipi);
        testIntervalOps(quad12, quad12, "TFTT", quad12, quad12);
        testIntervalOps(quad12, quad23, "FFTT", quad123, quad2);
        testIntervalOps(quad12, quad34, "FFTF", fullInterval, quad12);
        
        testIntervalOps(quad23, emptyInterval, "TTFF", quad23, emptyInterval);
        testIntervalOps(quad23, fullInterval, "FFTT", fullInterval, quad23);
        testIntervalOps(quad23, zero, "FFFF", quad234, emptyInterval);
        testIntervalOps(quad23, pi, "TTTT", quad23, pi);
        testIntervalOps(quad23, mipi, "TTTT", quad23, mipi);
        testIntervalOps(quad23, quad12, "FFTT", quad123, quad2);
        testIntervalOps(quad23, quad23, "TFTT", quad23, quad23);
        testIntervalOps(quad23, quad34, "FFTT", quad234, S1Interval(-S2.S2_M_PI, -S2.S2_M_PI/2));
        
        testIntervalOps(quad1, quad23, "FFTF", quad123, S1Interval(S2.S2_M_PI/2, S2.S2_M_PI/2));
        testIntervalOps(quad2, quad3, "FFTF", quad23, mipi);
        testIntervalOps(quad3, quad2, "FFTF", quad23, pi);
        testIntervalOps(quad2, pi, "TFTF", quad2, pi);
        testIntervalOps(quad2, mipi, "TFTF", quad2, mipi);
        testIntervalOps(quad3, pi, "TFTF", quad3, pi);
        testIntervalOps(quad3, mipi, "TFTF", quad3, mipi);
        
        testIntervalOps(quad12, mid12, "TTTT", quad12, mid12);
        testIntervalOps(mid12, quad12, "FFTT", quad12, mid12);
        
        quad12eps = S1Interval(quad12.lo, mid23.hi);
        quad2hi = S1Interval(mid23.lo, quad12.hi);
        testIntervalOps(quad12, mid23, "FFTT", quad12eps, quad2hi);
        testIntervalOps(mid23, quad12, "FFTT", quad12eps, quad2hi);
        
        # This test checks that the union of two disjoint intervals is the smallest
        # interval that contains both of them.  Note that the center of "mid34"
        # slightly CCW of -S2.S2_M_PI/2 so that there is no ambiguity about the result.
        quad412eps = S1Interval(mid34.lo, quad12.hi)
        testIntervalOps(quad12, mid34, "FFFF", quad412eps, emptyInterval);
        testIntervalOps(mid34, quad12, "FFFF", quad412eps, emptyInterval);
        
        quadeps12 = S1Interval(mid41.lo, quad12.hi);
        quad1lo = S1Interval(quad12.lo, mid41.hi);
        testIntervalOps(quad12, mid41, "FFTT", quadeps12, quad1lo);
        testIntervalOps(mid41, quad12, "FFTT", quadeps12, quad1lo);
        
        quad2lo = S1Interval(quad23.lo, mid12.hi);
        quad3hi = S1Interval(mid34.lo, quad23.hi);
        quadeps23 = S1Interval(mid12.lo, quad23.hi);
        quad23eps = S1Interval(quad23.lo, mid34.hi);
        quadeps123 = S1Interval(mid41.lo, quad23.hi);
        testIntervalOps(quad23, mid12, "FFTT", quadeps23, quad2lo);
        testIntervalOps(mid12, quad23, "FFTT", quadeps23, quad2lo);
        testIntervalOps(quad23, mid23, "TTTT", quad23, mid23);
        testIntervalOps(mid23, quad23, "FFTT", quad23, mid23);
        testIntervalOps(quad23, mid34, "FFTT", quad23eps, quad3hi);
        testIntervalOps(mid34, quad23, "FFTT", quad23eps, quad3hi);
        testIntervalOps(quad23, mid41, "FFFF", quadeps123, emptyInterval);
        testIntervalOps(mid41, quad23, "FFFF", quadeps123, emptyInterval);
    end

    begin
    @test addPoint(emptyInterval, 0) == zero
    @test addPoint(emptyInterval, S2.S2_M_PI) == pi
    @test addPoint(emptyInterval, -S2.S2_M_PI) == mipi
    @test addPoint(addPoint(emptyInterval, S2.S2_M_PI), -S2.S2_M_PI) == pi
    @test addPoint(addPoint(emptyInterval, -S2.S2_M_PI), S2.S2_M_PI) == pi
    @test addPoint(addPoint(emptyInterval, mid12.lo), mid12.hi) == mid12
    @test addPoint(addPoint(emptyInterval, mid23.lo), mid23.hi) == mid23
    @test addPoint(addPoint(quad1, -0.9 * S2.S2_M_PI), -S2.S2_M_PI/2) == quad123
    @test isfull(addPoint(fullInterval, 0))
    @test isfull(addPoint(fullInterval, S2.S2_M_PI))
    @test isfull(addPoint(fullInterval, -S2.S2_M_PI))
    end

    begin
        r = S1Interval(-S2.S2_M_PI, -S2.S2_M_PI);
        @test S2.S2_M_PI == clampPoint(r, -S2.S2_M_PI)
        @test S2.S2_M_PI == clampPoint(r, 0)
        r = S1Interval(0, S2.S2_M_PI);
        @test 0.1 == clampPoint(r, 0.1)
        @test 0.0 == clampPoint(r, -S2.S2_M_PI/2 + 1e-15)
        @test S2.S2_M_PI == clampPoint(r, -S2.S2_M_PI/2 - 1e-15)
        r = S1Interval(S2.S2_M_PI - 0.1, -S2.S2_M_PI + 0.1);
        @test S2.S2_M_PI == clampPoint(r, S2.S2_M_PI)
        @test S2.S2_M_PI - 0.1 == clampPoint(r, 1e-15)
        @test -S2.S2_M_PI + 0.1 == clampPoint(r, -1e-15)
        @test 0.0 == clampPoint(full(S1Interval), 0)
        @test S2.S2_M_PI == clampPoint(full(S1Interval), S2.S2_M_PI)
        @test S2.S2_M_PI == clampPoint(full(S1Interval), -S2.S2_M_PI)
    end

    @test fromPointPair(-S2.S2_M_PI, S2.S2_M_PI) == pi
    @test fromPointPair(S2.S2_M_PI, -S2.S2_M_PI) == pi
    @test fromPointPair(mid34.hi, mid34.lo) == mid34
    @test fromPointPair(mid23.lo, mid23.hi) == mid23


    @test expanded(emptyInterval, 1) == emptyInterval
    @test expanded(fullInterval, 1) == fullInterval
    @test expanded(zero, 1) == S1Interval(-1, 1)
    @test expanded(mipi, 0.01) == S1Interval(S2.S2_M_PI - 0.01, -S2.S2_M_PI + 0.01)
    @test expanded(pi, 27) == fullInterval
    @test expanded(pi, S2.S2_M_PI/2) == quad23
    @test expanded(pi2, S2.S2_M_PI/2) == quad12
    @test expanded(mipi2, S2.S2_M_PI/2) == quad34

    @test expanded(emptyInterval, -1) == emptyInterval
    @test expanded(fullInterval, -1) == fullInterval
    @test expanded(quad123, -27) == emptyInterval
    @test expanded(quad234, -27) == emptyInterval
    @test expanded(quad123, -S2.S2_M_PI/2) == quad2
    @test expanded(quad341, -S2.S2_M_PI/2) == quad4
    @test expanded(quad412, -S2.S2_M_PI/2) == quad1


    begin
        # Choose two values kLo and kHi such that it's okay to shift an endpoint by
        # kLo (i.e., the resulting interval is equivalent) but not by kHi.
        kLo = 4 * eps(Float64) # < maxError default
        kHi = 6 * eps(Float64) # > maxError default

        # Empty intervals.
        @test emptyInterval ≈ emptyInterval
        @test zero ≈ emptyInterval && emptyInterval ≈ zero
        @test pi ≈ emptyInterval && emptyInterval ≈ pi
        @test mipi ≈ emptyInterval && emptyInterval ≈ mipi
        @test !(emptyInterval ≈ fullInterval)
        @test emptyInterval ≈ S1Interval(1, 1 + 2 * kLo)
        @test !(emptyInterval ≈ S1Interval(1, 1 + 2 * kHi))
        @test S1Interval(S2.S2_M_PI - kLo, -S2.S2_M_PI + kLo) ≈ emptyInterval

        # Full intervals.
        @test fullInterval ≈ fullInterval
        @test !(fullInterval ≈ emptyInterval)
        @test !(fullInterval ≈ zero)
        @test !(fullInterval ≈ pi)
        @test fullInterval ≈ S1Interval(kLo, -kLo)
        @test !(fullInterval ≈ S1Interval(2 * kHi, 0))
        @test S1Interval(-S2.S2_M_PI + kLo, S2.S2_M_PI - kLo) ≈ fullInterval
        @test !(S1Interval(-S2.S2_M_PI, S2.S2_M_PI - 2 * kHi) ≈ fullInterval)

        # Singleton intervals.
        @test pi ≈ pi && mipi ≈ pi
        @test pi ≈ S1Interval(S2.S2_M_PI - kLo, S2.S2_M_PI - kLo)
        @test !(pi ≈ S1Interval(S2.S2_M_PI - kHi, S2.S2_M_PI - kHi))
        @test pi ≈ S1Interval(S2.S2_M_PI - kLo, -S2.S2_M_PI + kLo)
        @test !(pi ≈ S1Interval(S2.S2_M_PI - kHi, -S2.S2_M_PI))
        @test !(zero ≈ pi)
        @test union(union(pi, mid12), zero) ≈ quad12
        @test intersection(quad2, quad3) ≈ pi
        @test intersection(quad3, quad2) ≈ pi

        # Intervals whose corresponding endpoints are nearly the same but where the
        # endpoints are in opposite order (i.e., inverted intervals).
        @test !(S1Interval(0, kLo) ≈ S1Interval(kLo, 0))
        @test !(S1Interval(S2.S2_M_PI - 0.5 * kLo, -S2.S2_M_PI + 0.5 * kLo) ≈ S1Interval(-S2.S2_M_PI + 0.5 * kLo, S2.S2_M_PI - 0.5 * kLo))

        # Other intervals.
        @test S1Interval(1 - kLo, 2 + kLo) ≈ S1Interval(1, 2)
        @test S1Interval(1 + kLo, 2 - kLo) ≈ S1Interval(1, 2)
        @test S1Interval(2 - kLo, 1 + kLo) ≈ S1Interval(2, 1)
        @test S1Interval(2 + kLo, 1 - kLo) ≈ S1Interval(2, 1)
        @test !(S1Interval(1 - kHi, 2 + kLo) ≈ S1Interval(1, 2))
        @test !(S1Interval(1 + kHi, 2 - kLo) ≈ S1Interval(1, 2))
        @test !(S1Interval(2 - kHi, 1 + kLo) ≈ S1Interval(2, 1))
        @test !(S1Interval(2 + kHi, 1 - kLo) ≈ S1Interval(2, 1))
        @test !(S1Interval(1 - kLo, 2 + kHi) ≈ S1Interval(1, 2))
        @test !(S1Interval(1 + kLo, 2 - kHi) ≈ S1Interval(1, 2))
        @test !(S1Interval(2 - kLo, 1 + kHi) ≈ S1Interval(2, 1))
        @test !(S1Interval(2 + kLo, 1 - kHi) ≈ S1Interval(2, 1))
    end
    
    begin
        @test 0.0 == getDirectedHausdorffDistance(emptyInterval, emptyInterval)
        @test 0.0 == getDirectedHausdorffDistance(emptyInterval, mid12)
        @test S2.S2_M_PI == getDirectedHausdorffDistance(mid12, emptyInterval)

        @test 0.0 == getDirectedHausdorffDistance(quad12, quad123)
        in = S1Interval(3.0, -3.0) # an interval whose complement center is 0.
        @test 3.0 == getDirectedHausdorffDistance(S1Interval(-0.1, 0.2), in)
        @test 3.0 - 0.1 == getDirectedHausdorffDistance(S1Interval(0.1, 0.2), in)
        @test 3.0 - 0.1 == getDirectedHausdorffDistance(S1Interval(-0.2, -0.1), in)
    end
end
