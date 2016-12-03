# Inner and outer constructors

immutable Var{T}
end

immutable Conj{T}
end

typealias SymbolLike  Union{Symbol,AbstractString,Char}
typealias PolyLike    Union{Poly,Number}

immutable RationalFunction{T,S,U,V}
  num::Poly{U}
  den::Poly{V}

  # Full construction (from numerator and denominator polynomials)
  @compat function (::Type{RationalFunction}){S,U<:Number,V<:Number}(num::Poly{U},
    den::Poly{V}, ::Type{Conj{S}} = Conj{false})
    @assert isa(S, Bool)        "RationalFunction: S can be true or false in Conj{S}"
    @assert num.var == den.var  "RationalFunction: numerator and denominator have different variables"
    new{Var{den.var},Conj{S},U,V}(num,den)
  end

  # Partial construction (only from numerator polynomial)
  @compat function (::Type{RationalFunction}){S,U<:Number,V<:Number}(num::Poly{U},
    ::Type{V} = Float64, ::Type{Conj{S}} = Conj{false})
    @assert isa(S, Bool)        "RationalFunction: S can be true or false in Conj{S}"
    den = Poly([one(V)], num.var)
    new{Var{num.var},Conj{S},U,V}(num, den)
  end
end

# Full construction from vectors
RationalFunction{S}(num::Vector, den::Vector, var::SymbolLike = :x,
  t::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var), Poly(den, var), t)

# Partial construction from vectors
RationalFunction{S,V<:Number}(num::Vector, var::SymbolLike = :x,
  t1::Type{V} = Float64, t2::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var), t1, t2)

# Full construction from numbers
RationalFunction{S}(num::Number, den::Number, var::SymbolLike = :x,
  t::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var), Poly([den], var), t)

# Partial construction from numbers
RationalFunction{S,V<:Number}(num::Number, var::SymbolLike = :x,
  t1::Type{V} = Float64, t2::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var), t1, t2)

# Full construction from numbers and vectors
RationalFunction{S}(num::Vector, den::Number, var::SymbolLike = :x,
  t::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly(num, var), Poly([den], var), t)

RationalFunction{S}(num::Number, den::Vector, var::SymbolLike = :x,
  t::Type{Conj{S}} = Conj{false}) = RationalFunction(Poly([num], var), Poly(den, var), t)

# Full construction from numbers and `Poly`s
RationalFunction{S}(num::Number, den::Poly, t::Type{Conj{S}} = Conj{false}) =
  RationalFunction(Poly([num], den.var), den, t)

RationalFunction{S}(num::Poly, den::Number, t::Type{Conj{S}} = Conj{false}) =
  RationalFunction(num, Poly([den], num.var), t)

# Full construction from vectors and `Poly`s
RationalFunction{S}(num::Vector, den::Poly, t::Type{Conj{S}} = Conj{false}) =
  RationalFunction(Poly(num, den.var), den, t)

RationalFunction{S}(num::Poly, den::Vector, t::Type{Conj{S}} = Conj{false}) =
  RationalFunction(num, Poly(den, num.var), t)

# Construction from Poly division
function /(p1::Poly, p2::Poly)
  @assert p1.var == p2.var "/(p1,p2): Polynomials do not have the same variables"
  RationalFunction(p1, p2, Conj{false})
end

./(p1::Poly, p2::Poly) = /(p1,p2)

# Convert to Poly
function Poly(r::RationalFunction)
  @assert degree(r.den) == 0 "Poly(r): r.den is not constant"
  r.num / r.den[0]
end
