# RationalFunctions

A modest attempt to extend [`Polynomials`][poly-link] package to support rational
functions in [Julia][julia-link].

[poly-link]:  https://github.com/Keno/Polynomials.jl
[julia-link]: http://julialang.org/

### Build Status and Code Coverage

-  Build status: [![Build Status][build-img]][build-link]
-  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/aytekinar/RationalFunctions.jl.svg?branch=master
[build-link]: https://travis-ci.org/aytekinar/RationalFunctions.jl
[ca-img]: https://coveralls.io/repos/github/aytekinar/RationalFunctions.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/aytekinar/RationalFunctions.jl?branch=master
[cc-img]: https://codecov.io/gh/aytekinar/RationalFunctions.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/aytekinar/RationalFunctions.jl

### Description

The package implements rational functions (of single variable) which support
the following basic operations:

- [x] **Identities.** Methematical identities `one` and `zero`,
- [x] **Comparison.** Comparison related functions `hash`, `==`, `!=`, `isequal`,
- [x] **Operations.** Basic operations such as `inv`, `transpose`, `conj`, `derivative`,
  `+`, `-`, `*`, and `/`,
- [x] **Evaluation.** Function evaluation at a given `Number` or collections of
  `Number`s,
- [ ] **Function fitting.** Function fitting (of corresponding degrees) to a given
  set of `(input,output)` pairs, and,
- [ ] **Plotting.** Function plotting for a given set of input values (through
  [`RecipesBase`][recipes-link]).

The package is in its early phase. Any comments are appreciated, especially on

- [ ] how to implement `NaN`s and `Inf`s, and,
- [ ] **Zero functions.** When `num = zero(num)` or `den = Inf`.

[recipes-link]: https://github.com/JuliaPlots/RecipesBase.jl

### Note

The package (unfortunately) lacks

- Proper documentation, and,
- Unit testing.
