# Inner and outer constructors

"""
Alias for `Symbol`-like types: `AbstractString`, `Char`, `Symbol`.

See also: `RationalFunctions.PolynomialLike`.
"""
const SymbolLike = Union{AbstractString,Char,Symbol}

"""
Alias for `Polynomial-like types: `Number`, `Polynomial`.

See also: `RationalFunctions.SymbolLike`.
"""
const PolynomialLike = Union{Number,Polynomial}

"""
Constructor for `RationalFunction` objects.

Construct `RationalFunction` objects from `Polynomial` objects:

    RationalFunction(num, den = one(num), conj = Val{:notc})

where,
  * at least one of `num` and `den` objects is `Polynomial`, while the other can be either
    a `Number` or a `Vector`, and,
  * `conj` is a type which indicates whether the variable will be conjugated
    (`Val{:conj}`) in function evaluations, or not (`Val{:notc}`).

Construct `RationalFunction` objects by providing the coefficients (in ascending order):

    RationalFunction(num, den = one(eltype(num)), var = :x, conj = Val{:notc})

where,
  * `num` and `den` objects can be either a `Number` or a `Vector`, and,
  * `var` is either an `AbstractString`, a `Char` or a `Symbol`.

# Examples
```julia
julia> r1 = RationalFunction(poly([1,2,3]));
julia> r2 = RationalFunction(poly([1,2,3]), Polynomial([1,2,3]));
julia> r3 = RationalFunction(poly([1,2,3]), RationalFunctions.Val{:conj});
julia> r4 = RationalFunction([1,2,3]);
julia> r5 = RationalFunction(1, [1, 2, 3], "s");
julia> r6 = RationalFunction([1,2,3], 't', RationalFunctions.Val{:conj});
```

See also: `RationalFunctions.SymbolLike`, `RationalFunctions.PolynomialLike`, `coeffs`,
`degree`, `roots`, `variable`, `num`, `den`, `zeros`, `poles`, `funcfit`,
`derivative`, `reduce`, `solve`.
"""
struct RationalFunction{T,S,U,V}
  num::Polynomial{U}
  den::Polynomial{V}

  # Full construction (from numerator and denominator polynomials)
  @compat function (::Type{RationalFunction})(num::Polynomial{U}, den::Polynomial{V},
    ::Type{Val{:conj}}) where {U<:Number,V<:Number}
    if num.var ≠ den.var
      throw(DomainError((num, den), "RationalFunction(num,den): num and den `Polynomial`s have different variables"))
    end
    new{Val{num.var},Val{:conj},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction})(num::Polynomial{U}, den::Polynomial{V},
    ::Type{Val{:notc}}) where {U<:Number,V<:Number}
    if num.var ≠ den.var
      throw(DomainError((num, den), "RationalFunction(num,den): num and den `Polynomial`s have different variables"))
    end
    new{Val{num.var},Val{:notc},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction})(num::Polynomial{U}, den::Polynomial{V}) where
    {U<:Number,V<:Number}
    if num.var ≠ den.var
      throw(DomainError((num, den), "RationalFunction(num,den): num and den `Polynomial`s have different variables"))
    end
    new{Val{num.var},Val{:notc},U,V}(num,den)
  end
end

# Partial construction from `Polynomial`s
RationalFunction(num::Polynomial, conj::Type{Val{S}} = Val{:notc}) where {S} =
  RationalFunction(num, one(num), conj)

# Full construction from `Number`s, `Vector`s and `Polynomial`s
RationalFunction(num::Number, den::Polynomial, conj::Type{Val{S}} = Val{:notc}) where {S} =
  RationalFunction(Polynomial([num], den.var), den, conj)

RationalFunction(num::Polynomial, den::Number, conj::Type{Val{S}} = Val{:notc}) where {S} =
  RationalFunction(num, Polynomial([den], num.var), conj)

RationalFunction(num::Vector, den::Polynomial, conj::Type{Val{S}} = Val{:notc}) where {S} =
  RationalFunction(Polynomial(num, den.var), den, conj)

RationalFunction(num::Polynomial, den::Vector, conj::Type{Val{S}} = Val{:notc}) where {S} =
  RationalFunction(num, Polynomial(den, num.var), conj)

# Full construction from numbers and vectors
RationalFunction(num::Vector, den::Vector, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S} = RationalFunction(Polynomial(num, var), Polynomial(den, var), conj)

RationalFunction(num::Number, den::Number, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S} = RationalFunction(Polynomial([num], var), Polynomial([den], var), conj)

RationalFunction(num::Vector, den::Number, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S} = RationalFunction(Polynomial(num, var), Polynomial([den], var), conj)

RationalFunction(num::Number, den::Vector, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S} = RationalFunction(Polynomial([num], var), Polynomial(den, var), conj)

# Partial construction from numbers and vectors
RationalFunction(num::Vector{U}, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S,U<:Number} = RationalFunction(Polynomial(num, var),
  Polynomial([one(U)], var), conj)

RationalFunction(num::U, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) where {S,U<:Number} = RationalFunction(Polynomial([num], var),
  Polynomial([one(U)], var), conj)

"""
    Polynomial(r::RationalFunction)

Create a `Polynomial` object from `r` if `degree(den(r)) == 0`.
"""
function Polynomial(r::RationalFunction)
  if degree(r.den) ≠ 0
    throw(DomainError(r.den, "Polynomial(r): r.den is not constant"))
  end
  r.num / r.den[0]
end
