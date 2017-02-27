# Inner and outer constructors

"""
Alias for `Symbol`-like types: `AbstractString`, `Char`, `Symbol`.

See also: `RationalFunctions.PolyLike`.
"""
const SymbolLike = Union{AbstractString,Char,Symbol}

"""
Alias for `Poly`-like types: `Number`, `Poly`.

See also: `RationalFunctions.SymbolLike`.
"""
const PolyLike = Union{Number,Poly}

"""
Constructor for `RationalFunction` objects.

Construct `RationalFunction` objects from `Poly` objects:

    RationalFunction(num, den = one(num), conj = Val{:notc})

where,
  * at least one of `num` and `den` objects is `Poly`, while the other can be either
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
julia> r2 = RationalFunction(poly([1,2,3]), Poly([1,2,3]));
julia> r3 = RationalFunction(poly([1,2,3]), RationalFunctions.Val{:conj});
julia> r4 = RationalFunction([1,2,3]);
julia> r5 = RationalFunction(1, [1, 2, 3], "s");
julia> r6 = RationalFunction([1,2,3], 't', RationalFunctions.Val{:conj});
```

See also: `RationalFunctions.SymbolLike`, `RationalFunctions.PolyLike`, `coeffs`,
`degree`, `roots`, `variable`, `num`, `den`, `zeros`, `poles`, `funcfit`,
`derivative`, `reduce`, `solve`.
"""
immutable RationalFunction{T,S,U,V}
  num::Poly{U}
  den::Poly{V}

  # Full construction (from numerator and denominator polynomials)
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V}, ::Type{Val{:conj}})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Val{num.var},Val{:conj},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V}, ::Type{Val{:notc}})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Val{num.var},Val{:notc},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Val{num.var},Val{:notc},U,V}(num,den)
  end
end

# Partial construction from `Poly`s
RationalFunction{S}(num::Poly, conj::Type{Val{S}} = Val{:notc}) =
  RationalFunction(num, one(num), conj)

# Full construction from `Number`s, `Vector`s and `Poly`s
RationalFunction{S}(num::Number, den::Poly, conj::Type{Val{S}} = Val{:notc}) =
  RationalFunction(Poly([num], den.var), den, conj)

RationalFunction{S}(num::Poly, den::Number, conj::Type{Val{S}} = Val{:notc}) =
  RationalFunction(num, Poly([den], num.var), conj)

RationalFunction{S}(num::Vector, den::Poly, conj::Type{Val{S}} = Val{:notc}) =
  RationalFunction(Poly(num, den.var), den, conj)

RationalFunction{S}(num::Poly, den::Vector, conj::Type{Val{S}} = Val{:notc}) =
  RationalFunction(num, Poly(den, num.var), conj)

# Full construction from numbers and vectors
RationalFunction{S}(num::Vector, den::Vector, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly(num, var), Poly(den, var), conj)

RationalFunction{S}(num::Number, den::Number, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly([num], var), Poly([den], var), conj)

RationalFunction{S}(num::Vector, den::Number, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly(num, var), Poly([den], var), conj)

RationalFunction{S}(num::Number, den::Vector, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly([num], var), Poly(den, var), conj)

# Partial construction from numbers and vectors
RationalFunction{S,U<:Number}(num::Vector{U}, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly(num, var),
  Poly([one(U)], var), conj)

RationalFunction{S,U<:Number}(num::U, var::SymbolLike = :x,
  conj::Type{Val{S}} = Val{:notc}) = RationalFunction(Poly([num], var),
  Poly([one(U)], var), conj)

"""
    Poly(r::RationalFunction)

Create a `Poly` object from `r` if `degree(den(r)) == 0`.
"""
function Poly(r::RationalFunction)
  if degree(r.den) ≠ 0
    warn("Poly(r): r.den is not constant")
    throw(DomainError())
  end
  r.num / r.den[0]
end
