# S2.jl
A Julia port of the S2 library.

S2 provides a wide variety of spherical geometry primitives in addition
to an efficient form of geo-hashing.

This port provides the geo-hashing and a subset of the spherical
geometry, translated and slightly mutated into native Julia with
Juliesque idioms for key functions. The primary reference for the port is
the Java version of S2.

The changes in this version largely have to do with the switch from
object orientation to a function overloading style. For instance, where
have library would have multiple methods like `getCapBound` or
`getRectBound` in the `S2Region` interface, this version has a single
`bounds` function that dispatches on desired result type as well as the
type of the object for which bounds are needed. Thus you would have
`bounds(S2Cap, x)` instead of `x.getCapBound`.

A major impact of this is that we can use a native vector instead of the
`S2Point` type so that we can get + and * for free and we can use all of
the power of the existing `LinearAlgebra` module in base Julia. 
