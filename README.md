# RationalFunctions

[![Unix][unix-img]][unix-link] [![Windows][win-img]][win-link]
[![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

A modest attempt to extend [`Polynomials`][poly-link] package to support rational
functions in [Julia][julia-link].

[unix-img]: https://img.shields.io/travis/aytekinar/RationalFunctions.jl/master.svg?label=unix
[unix-link]: https://travis-ci.org/aytekinar/RationalFunctions.jl
[win-img]: https://img.shields.io/appveyor/ci/aytekinar/rationalfunctions-jl/master.svg?label=windows
[win-link]: https://ci.appveyor.com/project/aytekinar/rationalfunctions-jl/branch/master
[ca-img]: https://img.shields.io/coveralls/aytekinar/RationalFunctions.jl/master.svg?label=coveralls
[ca-link]: https://coveralls.io/github/aytekinar/RationalFunctions.jl?branch=master
[cc-img]: https://img.shields.io/codecov/c/github/aytekinar/RationalFunctions.jl/master.svg?label=codecov
[cc-link]: https://codecov.io/gh/aytekinar/RationalFunctions.jl?branch=master

[poly-link]:  https://github.com/Keno/Polynomials.jl
[julia-link]: http://julialang.org/

## Description

`RationalFunctions` aims at supporting rational functions of single variable. The
package extends `Polynomials`, and tries to provide a set of basic mathematical
functionality between `Number`s, `Poly`s and `RationalFunction` objects.

## Basic Usage

The easiest way to construct a `RationalFunction` object is to call its constructor
with `Poly` objects. If a `Poly` object is not provided, you need to provide the
variable name (defaults to `:x`). In either case, you should also provide the
conjugation property of the `RationalFunction` object (defaults to
`RationalFunctions.Conj{false}`).

As it is implemented now, `RationalFunction` objects do not allow for mixing
different variables --- a `RationalFunction` object's variable is the combination
of both its variable name (*i.e.*, `:x`) and its conjugation property (*i.e.*,
`RationalFunctions.Conj{false}`).

For more information, check the documentation with `?RationalFunction` command.
You will also see the interface, *i.e.*, exported methods of the package under
`See also:` section. You can read their documentation in a similar way.

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

julia> r1(1+1im)
0.2 - 0.39999999999999997im

julia> # Construct a rational function from numpoly, denpoly (conjugated input)
julia> r2 = RationalFunction(numpoly, denpoly, RationalFunctions.Conj{true})
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(-2 + 5⋅s - 4⋅s^2 + s^3), and,
den(s) is Poly(-6 + 11⋅s - 6⋅s^2 + s^3).

julia> r2(1+1im)
0.2 + 0.39999999999999997im

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
```

### Convenience functions

Some functions exist for your convenience when working with `RationalFunction`s.

Please read the corresponding documentation in Julia by issuing `?coeffs`, `?degree`,
`?roots`, `?variable`, `?num`, `?den`, `?zeros`, `?poles`, `?funcfit`, `?derivative`,
`?reduce`, `?residue`, or `?solve`.

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

julia> roots(r1) # Tuple of roots of num(r1) and den(r1)
(Complex{Float64}[2.0+0.0im,1.0+2.83263e-8im,1.0-2.83263e-8im],[3.0,2.0,1.0])

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

julia> r7 = RationalFunction(p1, p2)
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
WARNING: r1+r2: `r1` (s,Conj{false}) and `r2` (s,Conj{true}) have different variables
ERROR: DomainError: ...

julia> r1+r5
WARNING: r1+r2: `r1` (s,Conj{false}) and `r2` (x,Conj{false}) have different variables
ERROR: DomainError: ...

julia> r1+r3
f(s) = num(s)/den(s), where,
num(s) is Poly(-10.0 + 1.0⋅s - 2.0⋅s^2 + 38.0⋅s^3 - 36.0⋅s^4 + 9.0⋅s^5), and,
den(s) is Poly(-12.0 - 2.0⋅s - 4.0⋅s^2 + 44.0⋅s^3 - 32.0⋅s^4 + 6.0⋅s^5).

julia> r2+r4
f(s̄) = num(s̄)/den(s̄), where,
num(s) is Poly(-10.0 + 1.0⋅s - 2.0⋅s^2 + 38.0⋅s^3 - 36.0⋅s^4 + 9.0⋅s^5), and,
den(s) is Poly(-12.0 - 2.0⋅s - 4.0⋅s^2 + 44.0⋅s^3 - 32.0⋅s^4 + 6.0⋅s^5).

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

### Plotting of rational function outputs for given inputs

`RationalFunctions` provides plotting recipes through [`RecipesBase`][recipes-link]
for `Real`-coefficient `RationalFuncion` objects for `Real` inputs.

#### Example

The example below is run on a Linux machine with `Julia 0.5` having [`Plots`][plots-link]
as the frontend and [`GR`][gr-link] as the backend.

If you would like to use another (set of) package(s) for plotting, just disregard
the below code excerpt, and rely on the corresponding package's documentation and
the function evaluation implemented in `RationalFunctions`.

```julia
julia> using Plots; # Plots as the frontend, which is compatible with RecipesBase

julia> gr(); # GR as the backend

julia> x = -2:1E-1:2;

julia> xinit = x[5:5:end]; # x-values for function fitting

julia> yinit1 = map(x->x^2+3x+5, xinit); # y-values for function fitting

julia> yinit2 = map(x->x^2+2, xinit); # y-values for function fitting

julia> init1 = map((x,y)->(x,y), xinit, yinit1);

julia> init2 = map((x,y)->(x,y), xinit, yinit2);

julia> r1 = funcfit(xinit, yinit1, 2);

julia> r2 = funcfit(xinit, yinit2, 2);

julia> plot(r1, x, label = "r1(x)"); # plot r1 vs x

julia> plot!(r2, x, init2, label = "r2(x)"); # plot r2 vs x with given (x,y)-pairs scattered

julia> # plot!(r2, x, xinit, yinit2, label = "r2(x)") # same as above
```

![Figure obtained from the above code excerpt.][test-plot]

[recipes-link]: https://github.com/JuliaPlots/RecipesBase.jl
[plots-link]: https://github.com/tbreloff/Plots.jl
[gr-link]: https://github.com/jheinen/GR.jl
[test-plot]: https://github.com/aytekinar/RationalFunctions.jl/releases/download/v0.0.2/test-plot-linux-0.5.png
