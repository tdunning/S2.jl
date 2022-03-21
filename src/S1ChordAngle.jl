# The original Java version is Copyright 2014 Google Inc.
# This Julia version is Copyright 2022 Ted Dunning
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
An `S1ChordAngle` expresses an angle on a sphere but expresses this angle in terms
of the square of the length of the chord between the two points that define the angle.

This is equivalent of the `haversine` which was used before computers to simplify
spherical trigonometry computations and allows many common operations to be done
without trigonometry. For instance, the sine of a `S1ChordAngle` can be computed 
using a single square root operation. An `S1ChordAngle` is, however, limited to 
angles in the range from `[0, π]` which corresponds to the range `[0, 4]`.

For large angles, there can be some loss of accuracy. Near the maximum value, the 
representation of an `S1ChordAngle` can be expected to have an error of about
`(10^{-15} / x)`.


This file is an idiomatic translation of the Java implementation of S2
"""
struct S1ChordAngle
    distance2::Float64

    S1ChordAngle(distance2::Number) = S1ChordAngle(distance2, false)

    function S1ChordAngle(distance2::Number, force::Bool)
        if force
            new(distance2)
        else
            new(min(S1ChordAngle_MAX_LENGTH2, distance2))
        end
    end
    
    function S1ChordAngle(x::Vector, y::Vector)
        # assert (norm2(x) == 1 && norm2(y) == 1)
        # The distance may slightly exceed 4.0 due to roundoff errors
        S1ChordAngle(distance(x, y)^2)
    end

    # Returns a new chord angle approximated from {@code angle} (see 
    # getS1AngleConstructorMaxError() for the max magnitude of the error).
    #
    # Angles outside the range [0, Pi] are handled as follows:
    #
    # <ul>
    #   <li>{@link S1Angle#INFINITY} is mapped to {@link #INFINITY}
    #   <li>negative angles are mapped to {@link #NEGATIVE}
    #   <li>finite angles larger than Pi are mapped to {@link #STRAIGHT}
    # </ul>
    #
    # Note that this operation is relatively expensive and should be avoided. To use {@link
    # S1ChordAngle} effectively, you should structure your code so that input arguments are converted
    # to S1ChordAngles at the beginning of your algorithm, and results are converted back to `S1Angle`s
    # only at the end.
    
    function S1ChordAngle(angle::S1Angle) 
        if (radians(angle) < 0) 
            return S1ChordAngle_NEGATIVE;
        elseif (angle == S1Angle_INFINITY) 
            return S1ChordAngle_INFINITY
        else
            # The chord length is 2 * sin(angle / 2)
            length = 2 * sin(0.5 * min(π, radians(angle)))
            return S1ChordAngle(length * length);
        end
    end
end

const S1ChordAngle_MAX_LENGTH2 = 4.0
const S1ChordAngle_ZERO = S1ChordAngle(0)
const S1ChordAngle_RIGHT = S1ChordAngle(2)
const S1ChordAngle_STRAIGHT = S1ChordAngle(4)
const S1ChordAngle_INFINITY = S1ChordAngle(Inf, true)    
const S1ChordAngle_NEGATIVE = S1ChordAngle(-1)

isNegative(a::S1ChordAngle) = a.distance2 < 0
isInfinity(a::S1ChordAngle) = isinf(a.distance2)
isSpecial(a::S1ChordAngle) = isNegative(a) || isInfinity(a)
isZero(a::S1ChordAngle) = a.distance2 == 0


# Returns true if getLength2() is within the normal range of 0 to 4 (inclusive) or the angle is
# special.

function isValid(a::S1ChordAngle)::Bool
    return (a.distance2 >= 0 && a.distance2 <= MAX_LENGTH2) || isNegative(a) || isInfinity(a)
end


# Convert the chord angle to an {@link S1Angle}. {@link #INFINITY} is converted to {@link
# S1Angle#INFINITY}, and {@link #NEGATIVE} is converted to a negative {@link S1Angle}. This
# operation is relatively expensive.

function toAngle(a::S1ChordAngle)
    if (isNegative(a)) 
        return radians(-1);
    elseif (isInfinity(a)) 
        return S1Angle_INFINITY;
    else
        return radians(2 * asin(0.5 * sqrt(a.distance2)))
    end
end



# Returns a new S1ChordAngle whose chord distance represents the sum of the angular distances
# represented by the 'a' and 'b' chord angles.
#
# Note that this method is much more efficient than converting the chord angles to S1Angles
# and adding those. It requires only one square root plus a few additions and multiplications.

function (+)(a::S1ChordAngle, b::S1ChordAngle) 
    if isSpecial(a) || isSpecial(b)
        throw(ArgumentError("Can't add special chord angles"))
    end

    let a2 = a.distance2
        b2 = b.distance2
        begin
            # Optimization for the common case where "b" is an error tolerance parameter that happens to be
            # set to zero.
            if (b2 == 0) 
                return a
            end

            # Clamp the angle sum to at most 180 degrees.
            if (a2 + b2 >= S1ChordAngle_MAX_LENGTH2) 
                return S1ChordAngle_STRAIGHT
            end

            # Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
            # Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
            # Then the formula below can be derived from c = 2 * sin(A+B) and the relationships
            #   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
            #   cos(X) = sqrt(1 - sin^2(X)) .
            x = a2 * (1 - 0.25 * b2); # isValid() => non-negative
            y = b2 * (1 - 0.25 * a2); # isValid() => non-negative
            return S1ChordAngle(min(S1ChordAngle_MAX_LENGTH2, x + y + 2 * sqrt(x * y)))
        end
    end
end

md"""
Subtract one S1ChordAngle from another.

Note that this method is much more efficient than converting the chord angles to S1Angles
and adding those. It requires only one square root plus a few additions and multiplications.
"""
function (-)(a::S1ChordAngle, b::S1ChordAngle) :: S1ChordAngle
    # See comments in add(S1ChordAngle, S1ChordAngle).
    if isSpecial(a) || isSpecial(b)
        throw(ArgumentError("Can't subtract special chord angles"))
    end
    let a2 = a.distance2
        b2 = b.distance2
        begin
            if (b2 == 0) 
                return a
            end
            if (a2 <= b2) 
                return S1ChordAngle_ZERO;
            end
            let x = a2 * (1 - 0.25 * b2)
                y = b2 * (1 - 0.25 * a2)
                return S1ChordAngle(max(0.0, x + y - 2 * sqrt(x * y)))
            end
        end
    end
end

md"""
Returns the smaller of the given instances
"""
Base.min(a::S1ChordAngle, b::S1ChordAngle) = a.distance2 <= b.distance2 ? a : b


md"""
Returns the larger of the given instances
"""
Base.max(a::S1ChordAngle, b::S1ChordAngle) = a.distance2 > b.distance2 ? a : b

md"""
Returns the square of Math.sin(toAngle().radians()), but computed more efficiently
"""
function sin2(a::S1ChordAngle) 
    if isSpecial(a)
        throw(ArgumentError("Can't take the sin of a special angle"))
    end
    # Let "a" be the (non-squared) chord length, and let A be the corresponding half-angle
    # (a = 2*sin(A)). The formula below can be derived from:
    #   sin(2*A) = 2 * sin(A) * cos(A)
    #   cos^2(A) = 1 - sin^2(A)
    # This is much faster than converting to an angle and computing its sine.
    return a.distance2 * (1 - 0.25 * a.distance2)
end

md"""
Returns Math.sin(toAngle().radians()), but computed more efficiently
"""
Base.sin(a::S1ChordAngle) = sqrt(sin2(a))

md"""
Returns Math.cos(toAngle().radians()), but computed more efficiently
"""
function Base.cos(a::S1ChordAngle) 
    # cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
    if isSpecial(a)
        throw(ArgumentError("Can't take the cosine of a special angle"))
    end
              
    return 1 - 0.5 * a.distance2
end

md"""
Returns Math.tan(toAngle().radians()), but computed more efficiently.
"""
Base.tan(a::S1ChordAngle) = sin(a) / cos(a);

md"""
Returns a new S1ChordAngle that has been adjusted by the given error bound (which can be
positive or negative). {@code error} should be the value returned by one of the error bound
methods below. For example:

```
   {@code S1ChordAngle a = new S1ChordAngle(x, y);}
   {@code S1ChordAngle a1 = a.plusError(a.getS2PointConstructorMaxError());}
```

If this {@link #isSpecial}, we return {@code this}.
"""
function plusError(error::Float64) 
    return isSpecial() ? this : fromLength2(max(0.0, min(S1ChordAngle_MAX_LENGTH2, length2 + error)))
end

md"""
Returns the error in {@link #fromS1Angle}
"""
getS1AngleConstructorMaxError() = S2.DBL_EPSILON * length2

md"""
There is a relative error of {@code 2.5 * DBL_EPSILON} when computing the squared distance,
plus a relative error of {@code 2 * DBL_EPSILON} and an absolute error of {@code 16 *
DBL_EPSILON^2} because the lengths of the input points may differ from 1 by up to {@code 2 *
DBL_EPSILON} each. (This is the maximum length error in {@link S2Point#normalize}).
"""
getS2PointConstructorMaxError() = (4.5 * S2.DBL_EPSILON * length2) + (16 * S2.DBL_EPSILON * S2.DBL_EPSILON)

