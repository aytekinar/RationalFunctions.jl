# RationalFunctions

A modest attempt to extend [`Polynomials`][poly-link] package to support rational
functions in [Julia][julia-link].

[poly-link]:  https://github.com/Keno/Polynomials.jl
[julia-link]: http://julialang.org/

## Build Status and Code Coverage

-  Build status: [![Build Status][build-img]][build-link]
-  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/aytekinar/RationalFunctions.jl.svg?branch=master
[build-link]: https://travis-ci.org/aytekinar/RationalFunctions.jl
[ca-img]: https://coveralls.io/repos/github/aytekinar/RationalFunctions.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/aytekinar/RationalFunctions.jl?branch=master
[cc-img]: https://codecov.io/gh/aytekinar/RationalFunctions.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/aytekinar/RationalFunctions.jl

## Description

`RationalFunctions` aims at supporting rational functions of single variable. The
package extends `Polynomials`, and tries to provide a set of basic mathematical
functonality between `Number`s, `Poly`s and `RationalFunction` objects.

The package is in its early phase and provides the below functionality:

- [x] **Construction.** `RationalFunction` object creation from `Number`s, `Vector`s
  and `Poly`s, and `Poly` object creation from `RationalFunction`s (where appropriate),
- [x] **Conversion and promotion.** Conversion and promotion rules among `Number`s,
  `Poly`s and `RationalFunction`s (where appropriate),
- [x] **Display.** `MIME("text/plain")`- and `MIME("text/latex")`-aware object display,
- [x] **Convenience.** Some convenience functions such as `coeffs`, `degree`, `roots`,
  `variable`, `num`, `den`, `zeros` and `poles`,
- [x] **Identity.** Methematical identities `one` and `zero`,
- [x] **Comparison.** Comparison related functions `hash`, `==`, `isequal` and
  `isapprox`,
- [x] **Evaluation.** Function evaluation at a given `Number` or collections of
  `Number`s,
- [x] **Math.** Basic operations such as `inv`, `transpose`, `conj`, `derivative`,
  `reduce`, `+`, `-`, `*`, and `/`,
- [x] **Solution.** Solution of `lhs = rhs` through `solve(lhs, rhs)` where one
  of `lhs` and `rhs` is a `RationalFunction` whereas the other is a `Number`, `Poly`
  or `RationalFunction`,
- [ ] **Function fitting.** Function fitting (of corresponding degrees) to a given
  collection of `(input,output)` pairs, and,
- [ ] **Plotting.** Function plotting for a given collection of input values (through
  [`RecipesBase`][recipes-link]).

Any comments are appreciated, especially on

- [ ] how to implement `NaN`s and `Inf`s, and `isnan` and `isinf`, and,
- [ ] **Zero functions.** When `num = zero(num)` or `den = Inf`.

[recipes-link]: https://github.com/JuliaPlots/RecipesBase.jl

## Basic Usage

### Construction

The easiest way to construct a `RationalFunction` object is to call its constructor
with `Poly` objects. If a `Poly` object is not provided, you need to provide the
variable name (defaults to `:x`). In either case, you should also provide the
conjugation property of the `RationalFunction` object (defaults to
`RationalFunctions.Conj{false}`).

As it is implemented now, `RationalFunction` objects do not allow for mixing
different variables --- a `RationalFunction` object's variable is the combination
of both its variable name (*i.e.*, `:x`) and its conjugation property (*i.e.*,
`RationalFunctions.Conj{false}`).

#### Example

```julia
julia> using Polynomials
julia> using RationalFunctions

julia> numpoly = poly([1,1,2], :s) # construct (s-1)(s-1)(s-2)
Poly(-2 + 5⋅s - 4⋅s^2 + s^3)
julia> denpoly = poly([1,2,3], :s) # construct (s-1)(s-2)(s-3)
Poly(-6 + 11⋅s - 6⋅s^2 + s^3)

julia> # Construct a rational function from numpoly, denpoly
julia> r1 = RationalFunction(numpoly, denpoly)
f(s) = num(s)/den(s), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> # Construct a rational function from numpoly, denpoly (conjugated input)
julia> r2 = RationalFunction(numpoly, denpoly, RationalFunctions.Conj{true})
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> # Construct a rational function from coefficients
julia> r3 = RationalFunction([1,2,3],[2,4,6],:s)
f(s) = num(s)/den(s), where,
num(s) is Poly(1 + 2⋅s + 3⋅s^2), and,
den(s) is Poly(2 + 4⋅s + 6⋅s^2).

julia> # Construct a rational function from coefficients (conjugated input)
julia> r4 = RationalFunction([1,2,3],[2,4,6],:s,RationalFunctions.Conj{true})
f(s) = num(s)/den(s), where,
num(s) is Poly(1 + 2⋅s + 3⋅s^2), and,
den(s) is Poly(2 + 4⋅s + 6⋅s^2).

julia> # Construct a rational function from coefficients (dropping the variable)
julia> r5 = RationalFunction([1,2,3],[2,4,6])
f(s) = num(s)/den(s), where(s) is Poly(1 + 2⋅s + 3⋅s^2), and,
den(s) is Poly(2 + 4⋅s + 6⋅s^2).

julia> # Construct a directly from polynomial division
julia> r6 = numpoly/denpoly
f(s) = num(s)/den(s), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).
```

### Convenience functions

Some functions exist for your convenience when working with `RationalFunction`s.

#### Example

```julia
julia> r1
f(s) = num(s)/den(s), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> coeffs(r1) # Tuple of vector of coefficients of num(r1) and den(r1)
([-2,5,-4,1],[-6,11,-6,1])

julia> degree(r1) # Tuple of degrees of num(r1) and den(r1)
(3,3)

julia> roots(r1) # Tuple of zeros and poles of reduce(r1)
([1.0],[3.0])

julia> variable(r1) # Tuple of variables of num(r1), den(r1) and conjugation property
(Poly(s),Poly(s),RationalFunctions.Conj{false})

julia> num(r1) # Numerator in `Poly`
Poly(-2 + 5⋅s - 4⋅s^2 + s^3)

julia> den(r1) # Denominator in `Poly`
Poly(-6 + 11⋅s - 6⋅s^2 + s^3)

julia> zeros(r1) # Values which make reduce(r1) zero
1-element Array{Float64,1}:
 1.0

julia> poles(r1) # Values which make reduce(r1) ∞
1-element Array{Float64,1}:
 3.0
```

### Function evaluation

Using the usual call notation, you can evaluate function values at given input(s).

#### Example

```julia
julia> p1 = Poly(1+2*rand(2))
Poly(1.6809852898721749 + 2.613098491401878⋅x)
julia> p2 = Poly(1+2*rand(3))
Poly(1.6629832340509254 + 2.9921125048432287⋅x + 2.8500993637891843⋅x^2)

julia> r7 = p1/p2
f(x) = num(x)/den(x), where,
num(x) is Poly(1.6809852898721749 + 2.613098491401878⋅x), and,
den(x) is Poly(1.6629832340509254 + 2.9921125048432287⋅x + 2.8500993637891843⋅x^2).

julia> r8 = r7' # r8 = conj(r7)
f(x̄) = num(x̄)/den(x̄), where,
num(x) is Poly(1.6809852898721749 + 2.613098491401878⋅x), and,
den(x) is Poly(1.6629832340509254 + 2.9921125048432287⋅x + 2.8500993637891843⋅x^2).

julia> r7(1+1im)
0.4392153604476556 - 0.25879126601068836im

julia> r8(1+1im)
0.4392153604476556 + 0.25879126601068836im

julia> r7(randn(8)*1im)
8-element Array{Complex{Float64},1}:
 1.13511-0.379719im
 1.12115+0.444117im
 0.323077-0.693061im
 1.10553+0.153815im
 1.01852-0.625716im
 0.562162+0.784794im
 1.13827-0.288414im
 1.06107+0.0758621im

julia> r8(randn(5,5))
5×5 Array{Float64,2}:
 -0.484168   0.497893  -0.666764  -0.528562  0.646693
  0.503046   0.898892   0.50783   -0.668359  0.596057
  0.89041    0.810045   0.543107   0.309938  0.823568
  0.839318  -0.668416   0.979665   0.101942  0.389897
  1.00446    0.734795   0.991642  -0.668412  0.863274

julia> [r(x) for x in 1+9rand(5) + randn(5)*1im, r in [r7, r8]]
5×2 Array{Complex{Float64},2}:
 0.330996+0.0544501im    0.330996-0.0544501im
 0.342919-0.251789im     0.342919+0.251789im
 0.224866-0.000359367im  0.224866+0.000359367im
 0.102199-0.0100334im    0.102199+0.0100334im
 0.245507+0.137555im     0.245507-0.137555im
```

### Mathematical operations

You can combine `Number`s, `Poly`s and `RationalFunction`s (where appropriate) to
form mathematical expressions. Also, basic operations on `RationalFunction`s are
also defined.

#### Example

```julia
julia> r1
f(s) = num(s)/den(s), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> r2
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> r3
f(s) = num(s)/den(s), where,
num(s) is Poly(1 + 2⋅s + 3⋅s^2), and,
den(s) is Poly(2 + 4⋅s + 6⋅s^2).

julia> r4
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(1 + 2⋅s + 3⋅s^2), and,
den(s) is Poly(2 + 4⋅s + 6⋅s^2).

julia> r5
f(x) = num(x)/den(x), where,
num(x) is Poly(1 + 2⋅x + 3⋅x^2), and,
den(x) is Poly(2 + 4⋅x + 6⋅x^2).

julia> numpoly
Poly(-2 + 5⋅s - 4⋅s^2 + s^3)

julia> denpoly
Poly(-6 + 11⋅s - 6⋅s^2 + s^3)

julia> r1+r2
ERROR: r1+r2: r1 and r2 have different variables ((s,RationalFunctions.Conj{false}) vs (s,RationalFunctions.Conj{true}))

julia> r1+r5
ERROR: r1+r2: r1 and r2 have different variables ((s,RationalFunctions.Conj{false}) vs (x,RationalFunctions.Conj{false}))

julia> r1+r3
f(s) = num(s)/den(s), where,
num(s) is Poly(-10.0 + 1.0⋅s - 2.0⋅s^2 + 38.0⋅s^3 - 36.0⋅s^4 + 9.0⋅s^5), and,
den(s) is Poly(-12.0 - 2.0⋅s - 4.0⋅s^2 + 44.0⋅s^3 - 32.0⋅s^4 + 6.0⋅s^5).

julia> r2+r4
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(-10.0 + 1.0⋅s - 2.0⋅s^2 + 38.0⋅s^3 - 36.0⋅s^4 + 9.0⋅s^5), and,
den(s) is Poly(-12.0 - 2.0⋅s - 4.0⋅s^2 + 44.0⋅s^3 - 32.0⋅s^4 + 6.0⋅s^5).
julia> r1 == r6
true

julia> r1 * denpoly == numpoly
true

julia> r1 * denpoly / numpoly == 1
true

julia> reduce(r1)
f(s) = num(s)/den(s), where,
num(s) is Poly(-0.5 + 0.5⋅s), and,
den(s) is Poly(-1.5 + 0.5⋅s).

julia> reduce(r3)
f(s) = num(s)/den(s), where,
num(s) is Poly(0.5), and,
den(s) is Poly(1.0).

julia> reduce(r1+r3)
f(s) = num(s)/den(s), where,
num(s) is Poly(-1.25 + 0.75⋅s), and,
den(s) is Poly(-1.5 + 0.5⋅s).

julia> derivative(r1)
f(s) = num(s)/den(s), where,
num(s) is Poly(-8 + 24⋅s - 26⋅s^2 + 12⋅s^3 - 2⋅s^4), and,
den(s) is Poly(36 - 132⋅s + 193⋅s^2 - 144⋅s^3 + 58⋅s^4 - 12⋅s^5 + s^6).

julia> derivative(r1) |> reduce
f(s) = num(s)/den(s), where,
num(s) is Poly(1.0), and,
den(s) is Poly(-4.5 + 3.0⋅s - 0.5⋅s^2).

julia> derivative(r3)
f(s) = num(s)/den(s), where,
num(s) is Poly(0), and,
den(s) is Poly(4 + 16⋅s + 40⋅s^2 + 48⋅s^3 + 36⋅s^4).

julia> derivative(r3) |> reduce
f(s) = num(s)/den(s), where,
num(s) is Poly(0.0), and,
den(s) is Poly(1.0).
```

### Solution of rational function equalities

You can solve for the variable in `RationalFunction` equalities.

#### Example

```julia
julia> r1
f(s) = num(s)/den(s), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> solve(r1)
1-element Array{Float64,1}:
 1.0

julia> solve(r1, 5)
1-element Array{Float64,1}:
 3.5

julia> solve(r1, Inf)
1-element Array{Float64,1}:
 3.0

julia> solve(r1, 5*num(r1))
4-element Array{Float64,1}:
 1.0
 1.12111
 1.79085
 3.08803

julia> solve(r1, 4*r1+2)
1-element Array{Float64,1}:
 1.8

julia> solve(r1 * denpoly, 5*numpoly)
3-element Array{Complex{Float64},1}:
 1.0+3.46305e-8im
 1.0-3.46305e-8im
          2.0+0.0im
```

## Note

The package (unfortunately) lacks

- Proper documentation, and,
- Unit testing.

However, a proper documentation (at least as `?command` in Julia) will be provided
soon.
