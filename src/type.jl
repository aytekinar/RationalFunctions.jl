# Inner and outer constructors

"""
Value type for encoding the variable used in constructing `RationalFunction`s.

See also: `RationalFunctions.Conj`.
"""
immutable Var{T}
end

"""
Value type for encoding the conjugation property in constructing `RationalFunction`s.

See also: `RationalFunctions.Var`.
"""
immutable Conj{T}
end

"""
Alias for `Symbol`-like types: `Symbol`, `AbstractString`, `Char`.

See also: `RationalFunctions.PolyLike`.
"""
typealias SymbolLike  Union{Symbol,AbstractString,Char}

"""
Alias for `Poly`-like types: `Number`, `Poly`.

See also: `RationalFunctions.SymbolLike`.
"""
typealias PolyLike    Union{Poly,Number}

"""
Constructor for `RationalFunction` objects.

Construct `RationalFunction` objects from `Poly` objects:

    RationalFunction(num, den = one(num), conj = RationalFunctions.Conj{false})

where,
  * at least one of `num` and `den` objects is `Poly`, while the other can be either
    a `Number` or a `Vector`, and,
  * `conj` is a type which indicates whether the variable will be conjugated
    (`RationalFunctions.Conj{true}`) in function evaluations, or not (`RationalFunctions.Conj{false}`).

Construct `RationalFunction` objects by providing the coefficients (in ascending order):

    RationalFunction(num, den = one(eltype(num)), var = :x, conj = RationalFunctions.Conj{false})

where,
  * `num` and `den` objects can be either a `Number` or a `Vector`, and,
  * `var` is either a `Symbol`, a `Char` or an `AbstractString`.

# Examples
```julia
julia> r1 = RationalFunction(poly([1,2,3]));
julia> r2 = RationalFunction(poly([1,2,3]), Poly([1,2,3]));
julia> r3 = RationalFunction(poly([1,2,3]), RationalFunctions.Conj{true});
julia> r4 = RationalFunction([1,2,3]);
julia> r5 = RationalFunction(1, [1, 2, 3], "s");
julia> r6 = RationalFunction([1,2,3], 't', RationalFunctions.Conj{true});
```

See also: `RationalFunctions.Var`, `RationalFunctions.Conj`, `RationalFunctions.SymbolLike`,
`RationalFunctions.PolyLike`, `coeffs`, `degree`, `roots`, `variable`, `num`, `den`,
`zeros`, `poles`, `funcfit`, `derivative`, `reduce`, `solve`.
"""
immutable RationalFunction{T,S,U,V}
  num::Poly{U}
  den::Poly{V}

  # Full construction (from numerator and denominator polynomials)
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V}, ::Type{Conj{true}})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Var{num.var},Conj{true},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V}, ::Type{Conj{false}})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Var{num.var},Conj{false},U,V}(num,den)
  end
  @compat function (::Type{RationalFunction}){U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V})
    if num.var ≠ den.var
      warn("RationalFunction(num,den): num and den `Poly`s have different variables")
      throw(DomainError())
    end
    new{Var{num.var},Conj{false},U,V}(num,den)
  end
end

# Partial construction from `Poly`s
RationalFunction{S}(num::Poly, conj::Type{Conj{S}} = Conj{false}) =
  RationalFunction(num, one(num), conj)

# Full construction from `Number`s, `Vector`s and `Poly`s
RationalFunction{S}(num::Number, den::Poly, conj::Type{Conj{S}} = Conj{false}) =
  RationalFunction(Poly([num], den.var), den, conj)

RationalFunction{S}(num::Poly, den::Number, conj::Type{Conj{S}} = Conj{false}) =
  RationalFunction(num, Poly([den], num.var), conj)

RationalFunction{S}(num::Vector, den::Poly, conj::Type{Conj{S}} = Conj{false}) =
  RationalFunction(Poly(num, den.var), den, conj)

RationalFunction{S}(num::Poly, den::Vector, conj::Type{Conj{S}} = Conj{false}) =
  RationalFunction(num, Poly(den, num.var), conj)

# Full construction from numbers and vectors
RationalFunction{S}(num::Vector, den::Vector, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var), Poly(den, var), conj)

RationalFunction{S}(num::Number, den::Number, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var), Poly([den], var), conj)

RationalFunction{S}(num::Vector, den::Number, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var), Poly([den], var), conj)

RationalFunction{S}(num::Number, den::Vector, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var), Poly(den, var), conj)

# Partial construction from numbers and vectors
RationalFunction{S,U<:Number}(num::Vector{U}, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var),
  Poly([one(U)], var), conj)

RationalFunction{S,U<:Number}(num::U, var::SymbolLike = :x,
  conj::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var),
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
