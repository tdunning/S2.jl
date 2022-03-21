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

md"""
An S1Interval represents a closed interval on a unit circle (also known as a 1-dimensional
sphere). It is capable of representing the empty interval (containing no points), the full
interval (containing all points), and zero-length intervals (containing a single point).

Points are represented by the angle they make with the positive x-axis in the range [-π, π].
An interval is represented by its lower and upper bounds (both inclusive, since the interval is
closed). The lower bound may be greater than the upper bound, in which case the interval is
"inverted" (i.e. it passes through the point (-1, 0)).

Note that the point (-1, 0) has two valid representations, π and -π. The normalized
representation of this point internally is π, so that endpoints of normal intervals are in the
range (-π, π]. However, we take advantage of the point -π to construct two special intervals:
the full() interval is [-π, π], and the Empty() interval is [π, -π].
"""
struct S1Interval
    lo::Float64
    hi::Float64

    S1Interval(low::Number, high::Number) = S1Interval(low, high, false)
    S1Interval(a::S1Interval) = S1Interval(a.lo, a.hi)
    function S1Interval(lowx::Number, highx::Number, checked::Bool)
        low, high = Float64(lowx), Float64(highx)
        lo, hi = low, high
        if !checked
            if low == -S2_M_PI && high != S2_M_PI
                lo = S2_M_PI
            end
            if high == -S2_M_PI && low != S2_M_PI
                hi = S2_M_PI
            end
        end
        new(lo, hi)
    end
end



md"""The empty interval"""
Base.empty(T::Type{S1Interval}) = S1Interval(S2_M_PI, -S2_M_PI, true)

md"""The full interval"""
full(T::Type{S1Interval}) = S1Interval(-S2_M_PI, S2_M_PI, true)

md"Construct an interval containing a single point"
function fromPoint(p::Number) 
    if p == -S2_M_PI
        p = S2_M_PI
    end
    
    return S1Interval(p, p)
end

md"""
Convenience method to construct the minimal interval containing the two given points. This is
equivalent to starting with an empty interval and calling addPoint() twice, but it is more
efficient.
"""
function fromPointPair(p1::Number, p2::Number)::S1Interval 
    if (p1 == -S2_M_PI) 
        p1 = S2_M_PI
    end
    if (p2 == -S2_M_PI) 
        p2 = S2_M_PI
    end
    if (positiveDistance(p1, p2) <= S2_M_PI) 
        S1Interval(p1, p2)
    else 
        S1Interval(p2, p1)
    end
end

md"""
An interval is valid if neither bound exceeds π in absolute value, and the value -π appears
only in the empty() and full() intervals.
"""
Base.isvalid(a::S1Interval) = abs(a.lo) <= S2_M_PI && abs(a.hi) <= S2_M_PI &&
    !(a.lo == -S2_M_PI && a.hi != S2_M_PI) && !(a.hi == -S2_M_PI && a.lo != S2_M_PI)


md"Returns true if the interval contains all points on the unit circle"
isfull(a::S1Interval) = a.hi - a.lo == 2 * S2_M_PI

md"Returns true if the interval is empty, i.e. it contains no points"
Base.isempty(a::S1Interval) = a.lo - a.hi == 2 * S2_M_PI

md"Returns true if lo() > hi(). (This is also true for empty intervals)"
isinverted(a::S1Interval) = a.lo > a.hi

md"""
Returns the midpoint of the interval. For full and empty intervals, the result is arbitrary
"""
function center(a::S1Interval) 
    center = 0.5 * (a.lo + a.hi)
    if !isinverted(a)
        return center
    end
    # Return the center in the range (-π, π].
    return (center <= 0) ? (center + S2_M_PI) : (center - S2_M_PI)
end

md"Returns the length of the interval. The length of an empty interval is negative"
function Base.length(a::S1Interval)
    length = a.hi - a.lo
    if length >= 0
      return length
    end
    length += 2 * S2_M_PI
    # Empty intervals have a negative length.
    return (length > 0) ? length : -1
end

md"""
Return the complement of the interior of the interval. An interval and its complement have the
same boundary but do not share any interior values. The complement operator is not a bijection,
since the complement of a singleton interval (containing a single value) is the same as the
complement of an empty interval.
"""
function complement(a::S1Interval) 
    if (a.lo == a.hi) 
        return full(S1Interval)
    end
    # Handles empty and full.
    return S1Interval(a.hi, a.lo, true)
end

md"""
Return the midpoint of the complement of the interval. For full and empty intervals, the result
is arbitrary. For a singleton interval (containing a single point), the result is its antipodal
point on S1.
"""
function complementcenter(a::S1Interval) 
    if a.lo != a.hi 
        return center(complement(a))
    else 
        # Singleton
        return (a.hi <= 0) ? (a.hi + S2_M_PI) : (a.hi - S2_M_PI)
    end
end

md"Returns true if the interval (which is closed) contains the point `p`"
function Base.contains(a::S1Interval, px::Number)
    p = Float64(px)
    # Works for empty, full, and singleton intervals.
    if p == -S2_M_PI
        p = S2_M_PI
    end
    if (isinverted(a)) 
        return (p >= a.lo || p <= a.hi) && !isempty(a)
    else 
        return a.lo <= p <= a.hi
    end
end

md"Returns true if the interior of the interval contains the point `p`"
function interiorcontains(a::S1Interval, px::Number)
    p = Float64(px)
    # Works for empty, full, and singleton intervals.
    if p == -S2_M_PI
        p = S2_M_PI
    end

    if isinverted(a) 
        return p > a.lo || p < a.hi
    else
        return (a.lo < p < a.hi) || isfull(a)
    end
end

md"""
Returns true if the interval `a` contains the interval `b`. Works for empty, full, and
singleton intervals.
"""
function Base.contains(a::S1Interval, b::S1Interval) 
    # It might be helpful to compare the structure of these tests to
    # the simpler Contains(double) method above.
    if isinverted(a)
        if isinverted(b)
            return b.lo >= a.lo && b.hi <= a.hi
        end
        return (b.lo >= a.lo || b.hi <= a.hi) && !isempty(a)
    else 
        if isinverted(b)
            return isfull(a) || isempty(b)
        end
    end
    return b.lo >= a.lo && b.hi <= a.hi
end


md"""
Returns true if the interior of interval `a` contains the entire interval `b`. Note that
`interiorcontains(x, x)` is true only when `x` is the empty or full interval, and
`interiorcontains(x, S1Interval(p,p))` is equivalent to `InteriorContains(x, p)`.
"""
function interiorcontains(a::S1Interval, b::S1Interval)
    if isinverted(a)
        if !isinverted(b)
            return b.lo > a.lo || b.hi < a.hi
        end
        return (b.lo > a.lo && b.hi < a.hi) || isempty(b)
    else 
        if isinverted(b)
            return isfull(a) || isempty(b)
        end
        return (b.lo > a.lo && b.hi < a.hi) || isfull(a)
    end
end

md"""
Returns true if the two intervals contain any points in common. Note that the point +/-π has
two representations, so the intervals [-π,-3] and [2,π] intersect, for example.
"""
function intersects(a::S1Interval, b::S1Interval)
    if (isempty(a) || isempty(b)) 
        return false
    end
    if isinverted(a) 
        # Every non-empty inverted interval contains π
        return isinverted(b) || b.lo <= a.hi || b.hi >= a.lo
    else 
      if isinverted(b)
          return b.lo <= a.hi || b.hi >= a.lo
      end
      return b.lo <= a.hi && b.hi >= a.lo
    end
end

md"""
Returns true if the interior of this interval contains any point of the interval {@code y}
(including its boundary). Works for empty, full, and singleton intervals.
"""
function interiorintersects(a::S1Interval, b::S1Interval)
    if isempty(a) || isempty(b) || a.lo == a.hi
        return false
    end
    if isinverted(a)
        return isinverted(b) || b.lo < a.hi || b.hi > a.lo
    else
        if isinverted(b)
            return b.lo < a.hi || b.hi > a.lo
        end
      return (b.lo < a.hi && b.hi > a.lo) || isfull(a)
    end
end

md"""
Return the Hausdorff distance from `a` to interval `b`. For two S1Intervals `a` and `b`,
this distance is defined by `h(x, y) = max_{p in x} min_{q in y} d(p, q)`, where `d(.,.)` 
is measured along `S1`.
"""
function getDirectedHausdorffDistance(a::S1Interval, b::S1Interval)
    if contains(b, a)
        # this accounts for empty `a`
        return 0.0
    end
    if isempty(b)
        # max distance in `S1`
        return S2_M_PI
    end

    bComplementCenter = complementcenter(b)
    if contains(a, bComplementCenter)
        return positiveDistance(b.hi, bComplementCenter)
    else
        # The Hausdorff distance is realized by either two hi endpoints or two
        # lo endpoints, whichever is farther apart.
        hiHi = contains(S1Interval(b.hi, bComplementCenter), a.hi) ? positiveDistance(b.hi, a.hi) : 0
        loLo = contains(S1Interval(bComplementCenter, b.lo), a.lo) ? positiveDistance(a.lo, b.lo) : 0
        return max(hiHi, loLo)
    end
end

md"""
Expands the interval by the minimum amount necessary so that it contains the point {@code p}
(an angle in the range [-π, π]).
"""
function addPoint(a::S1Interval, p::Number)
    if p == -S2_M_PI
        p = S2_M_PI
    end

    if contains(a, p)
        return a
    end

    if isempty(a)
        return fromPoint(p)
    else
        # Compute distance from p to each endpoint.
        dlo = positiveDistance(p, a.lo)
        dhi = positiveDistance(a.hi, p)
        if (dlo < dhi) 
            return S1Interval(p, a.hi)
        else
            return S1Interval(a.lo, p)
        end
        # Adding a point can never turn a non-full interval into a full one.
    end
end

md"""
Returns the closest point in the interval to the point {@code p}. The interval must be
non-empty.
"""
function clampPoint(a::S1Interval, p::Number) 
    if p == -S2_M_PI
        p = S2_M_PI
    end

    if contains(a, p)
        return p
    end
    
    # Compute distance from p to each endpoint.
    dlo = positiveDistance(p, a.lo)
    dhi = positiveDistance(a.hi, p)
    return (dlo < dhi) ? a.lo : a.hi
end

md"""
Returns a new interval that has been expanded on each side by the distance {@code margin}. If
"margin" is negative, then shrink the interval on each side by "margin" instead. The resulting
interval may be empty or full. Any expansion (positive or negative) of a full interval remains
full, and any expansion of an empty interval remains empty.
"""
function expanded(a::S1Interval, margin::Number)
    if margin >= 0
        if isempty(a)
            return a
        end
        # Check whether this interval will be full after expansion, allowing
        # for a 1-bit rounding error when computing each endpoint.
        if length(a) + 2 * margin + 2 * eps(Float64) >= 2 * S2_M_PI 
            return full(S1Interval)
        end
    else
        if isfull(a)
            return a
        end
        # Check whether this interval will be empty after expansion, allowing
        # for a 1-bit rounding error when computing each endpoint.
        if length(a) + 2 * margin - 2 * eps(Float64) <= 0 
            return empty(S1Interval)
        end
    end

    lo = rem(a.lo - margin, 2 * S2_M_PI, RoundNearest)
    hi = rem(a.hi + margin, 2 * S2_M_PI, RoundNearest)
    
    return S1Interval(lo, hi)
end

md"""
Returns the smallest interval that contains interval `a` and the interval `b`.
"""
function Base.union(a::S1Interval, b::S1Interval)
    lo = a.lo
    hi = a.hi
    
    # The `isfull(b)` case is handled correctly in all cases by the code
    # below, but can follow three separate code paths depending on whether
    # `a` is inverted, is non-inverted but contains π, or neither.
    if !isempty(b)
        if contains(a, b.lo)
            if contains(a, b.hi)
                # Either this interval contains y, or the union of the two
                # intervals is the full interval.
                if (!contains(a, b)) 
                    return full(S1Interval)
                end
            else
                hi = b.hi
            end
        elseif contains(a, b.hi)
            lo = b.lo
        elseif isempty(a) || contains(b, lo)
            # Interval `a` contains neither endpoint of `b`. This means that either `b`
            # contains all of `a`, or the two intervals are disjoint.
            lo = b.lo
            hi = b.hi
        else 
            # Check which pair of endpoints are closer together.
            dlo = positiveDistance(b.hi, lo)
            dhi = positiveDistance(hi, b.lo)
            if dlo < dhi 
                lo = b.lo
            else
                hi = b.hi
            end
        end
    end
    return S1Interval(lo, hi)
end


md"""
Returns the smallest interval that contains the intersection of interval `a` with `b`
Note that the region of intersection may consist of two disjoint intervals.
"""
function intersection(a::S1Interval, b::S1Interval)
    lo, hi = a.lo, a.hi
    # The `isfull(b)` case is handled correctly in all cases by the code below, but can follow three
    # separate code paths depending on whether this interval is inverted, is non-inverted but
    # contains π, or neither.

    if isempty(b) 
        return empty(S1Interval)
    elseif contains(a, b.lo)
        if contains(a, b.hi)
            # Either `a` contains `b`, or the region of intersection consists of two disjoint
            # subintervals. In either case, we want to set the interval to the shorter of the two
            # original intervals.
            if length(b) < length(a)
                return S1Interval(b.lo, b.hi, true) # isfull(a) code path
            end
        else
            return S1Interval(b.lo, hi, true)
        end
    elseif contains(a, b.hi)
        return S1Interval(lo, b.hi, true)
    else
        # Interval `a` contains neither endpoint of `b`. This means that either `b`
        # contains all of this interval, or the two intervals are disjoint.
        if !contains(b, lo)
            return empty(S1Interval)
        end
    end
    return S1Interval(lo, hi)
end

md"""
Returns true if this interval can be transformed into the interval {@code y} by moving each
endpoint by at most "maxError" (and without the endpoints crossing, which would invert the
interval). Empty and full intervals are considered to start at an arbitrary point on the unit
circle, thus any interval with (length <= 2*maxError) matches the empty interval, and any
interval with (length >= 2*π - 2*maxError) matches the full interval.
"""
function Base.isapprox(a::S1Interval, b::S1Interval; atol=1e-15)
    # Full and empty intervals require special cases because the "endpoints"
    # are considered to be positioned arbitrarily.
    if isempty(a)
        return length(b) <= 2 * atol
    end
    if isempty(b)
        return length(a) <= 2 * atol
    end
    if isfull(a)
        return length(b) >= 2 * (S2_M_PI - atol)
    end
    if isfull(b)
        return length(a) >= 2 * (S2_M_PI - atol)
    end

    # The purpose of the last test below is to verify that moving the endpoints
    # does not invert the interval, e.g. [-1e20, 1e20] vs. [1e20, -1e20].
    return abs(rem(b.lo - a.lo, 2 * S2_M_PI, RoundNearest)) <= atol &&
        abs(rem(b.hi - a.hi, 2 * S2_M_PI, RoundNearest)) <= atol &&
        abs(length(a) - length(b)) <= 2 * atol
end

function isequals(a::S1Interval, b::S1Interval) 
    return a.lo == b.lo && a.hi == b.hi
end

md"""
Computes the distance from {@code a} to {@code b} in the range [0, 2*π). This is equivalent to
{@code drem(b - a - π, 2 * π) + π}, except that it is more numerically stable
(it does not lose precision for very small positive distances).
"""
function positiveDistance(a::Number, b::Number)
    d = Float64(b - a);
    if d >= 0
        return d
    else
        # We want to ensure that if b == π and a == (-π + eps),
        # the return result is approximately 2π and not zero.
        return (b + S2_M_PI) - (a - S2_M_PI)
    end
end

