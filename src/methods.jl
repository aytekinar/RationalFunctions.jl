# Copy operation
copy{T,S}(r::RationalFunction{T,S}) = RationalFunction(copy(r.num), copy(r.den), S)

# Iteration interface
eltype{T,S,U,V}(::Type{RationalFunction{T,S,U,V}}) = (U, V)
eltype{T,S,U,V}(r::RationalFunction{T,S,U,V})      = (U, V)

# Convenience functions
coeffs(r::RationalFunction) = (coeffs(r.num), coeffs(r.den))

degree(r::RationalFunction) = (degree(r.num), degree(r.den))

roots(r::RationalFunction)  = (rnew = reduce(r); (roots(rnew.num), roots(rnew.den)))

variable{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}}) =
  (variable(U, T), variable(V, T), Conj{S})
variable{T,S,U,V}(r::RationalFunction{Var{T},S,U,V})            =
  (variable(U, T), variable(V, T), S)

num(r::RationalFunction)    = r.num
den(r::RationalFunction)    = r.den
zeros(r::RationalFunction)  = (rnew = reduce(r); roots(rnew.num))
poles(r::RationalFunction)  = (rnew = reduce(r); roots(rnew.den))

## Identities
one{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}})  =
  RationalFunction(Poly([one(U)], T), Poly([one(V)], T), Conj{S})
one{T,S}(r::RationalFunction{T,S})                          =
  RationalFunction(one(r.num), one(r.den), S)
zero{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}}) =
  RationalFunction(Poly(U[], T), Poly([one(V)], T), Conj{S})
zero{T,S}(r::RationalFunction{T,S})                         =
  RationalFunction(zero(r.num), one(r.den), S)

## Comparison
hash{T,S}(r::RationalFunction{Var{T},Conj{S}}, h::UInt)             =
  hash(T, hash(S, hash(coeffs(r.num), hash(coeffs(r.den), h))))

=={T,S}(r1::RationalFunction{T,S}, r2::RationalFunction{T,S})       =
  r1.num * r2.den == r1.den * r2.num
==(r1::RationalFunction, r2::RationalFunction)                      = false

isequal{T,S}(r1::RationalFunction{T,S}, r2::RationalFunction{T,S})  = hash(r1) == hash(r2)
isequal(r1::RationalFunction, r2::RationalFunction)                 = false

# Auxiliary definitions for `isapprox`
eps{T<:AbstractFloat}(::Type{T})          = Base.eps(T)
eps{T<:AbstractFloat}(::Type{Complex{T}}) = Base.eps(T)
eps{T}(::Type{T})                         = zero(T)

function isapprox{T,S,U1,V1,U2,V2}(r1::RationalFunction{Var{T},Conj{S},U1,V1},
  r2::RationalFunction{Var{T},Conj{S},U2,V2};
  rtol::Real = sqrt(eps(promote_type(U1,V2,U2,V2))), atol::Real = 0)
  p1 = r1.num * r2.den
  p2 = r1.den * d2.num

  isapprox(coeffs(p1), coeffs(p2); rtol = rtol, atol = atol)
end

function isapprox{T1,S1,T2,S2}(r1::RationalFunction{Var{T1},Conj{S1}},
  r2::RationalFunction{Var{T2},Conj{S2}}; rtol::Real = 0, atol::Real = 0)
  warn("r1≈r2: `r1` ($T1,Conj{$S1}) and `r2` ($T2,Conj{$S2}) have different variables")
  throw(DomainError())
end

# Function evaluation
function _funceval(r::RationalFunction, x::Number)
  mindegree = min(degree(r)...)
  result    = r.num(x)/r.den(x)
  !isnan(result) && return result

  k         = 1
  num       = r.num
  den       = r.den
  while isnan(result) && k ≤ mindegree
    # Apply L'Hospital
    num     = polyder(num)
    den     = polyder(den)
    result  = num(x)/den(x)
    k      += 1
  end
  result
end

_funceval(r::RationalFunction, X) = map(x->_funceval(r, x), X)

@compat (r::RationalFunction)(x::Number)                    = _funceval(r, x)
@compat (r::RationalFunction{T,Conj{true}}){T}(x::Number)   = _funceval(r, conj(x))
@compat (r::RationalFunction{T,Conj{true}}){T}(x::Real)     = _funceval(r, x)
@compat (r::RationalFunction)(X)                            = _funceval(r, X)

# Mathematical operations (always return temporaries for correctness of the results)
## Inversion
inv{T,S}(r::RationalFunction{T,S})  = RationalFunction(copy(r.den), copy(r.num), S)

## Transposition
transpose(r::RationalFunction)      = copy(r)

## Conjugation
function conj{T,S}(r::RationalFunction{Var{T},Conj{S}})
  numcoeff, dencoeff = coeffs(r)
  RationalFunction(Poly(conj(copy(numcoeff)), T), Poly(conj(copy(dencoeff)), T),
    Conj{!S})
end

## Derivative
function derivative{T,S}(r::RationalFunction{T,S}, n::Int = 1)
  if n < 0
    warn("derivative(r, n): `n` must be non-negative")
    throw(DomainError())
  end
  n == 0 && return copy(r)
  num   = polyder(r.num)*r.den - r.num*polyder(r.den)
  den   = r.den*r.den
  temp  = RationalFunction(num, den, S)
  for count in 2:n
    num   = polyder(temp.num)*temp.den - temp.num*polyder(temp.den)
    den   = temp.den*temp.den
    temp  = RationalFunction(num, den, S)
  end
  temp
end

## Reduction
function reduce{T,S}(r::RationalFunction{T,S})
  g       = gcd(r.num, r.den)
  common  = degree(g) == 0 ? one(g) : g

  num, _  = divrem(r.num, common)
  den, _  = divrem(r.den, common)
  RationalFunction(num, den, S)
end

## Basic operations between rational functions
function +{T,S}(r1::RationalFunction{Var{T},Conj{S}}, r2::RationalFunction{Var{T},Conj{S}})
  g       = gcd(r1.den, r2.den)
  common  = Polynomials.degree(g) == 0 ? one(g) : g

  den1, _ = divrem(r1.den, common)
  den2, _ = divrem(r2.den, common)
  num     = r1.num * den2 + den1 * r2.num
  den     = den1*den2*common
  RationalFunction(num, den, Conj{S})
end

function +{T1,S1,T2,S2}(r1::RationalFunction{Var{T1},Conj{S1}}, r2::RationalFunction{Var{T2},Conj{S2}})
  warn("r1+r2: `r1` ($T1,Conj{$S1}) and `r2` ($T2,Conj{$S2}) have different variables")
  throw(DomainError())
end

function *{T,S}(r1::RationalFunction{Var{T},Conj{S}}, r2::RationalFunction{Var{T},Conj{S}})
  num = r1.num * r2.num
  den = r1.den * r2.den
  RationalFunction(num, den, Conj{S})
end

function *{T1,S1,T2,S2}(r1::RationalFunction{Var{T1},Conj{S1}}, r2::RationalFunction{Var{T2},Conj{S2}})
  warn("r1*r2: `r1` ($T1,Conj{$S1}) and `r2` ($T2,Conj{$S2}) have different variables")
  throw(DomainError())
end

dot(r1::RationalFunction, r2::RationalFunction) = *(r1, r2)

function /{T,S}(r1::RationalFunction{Var{T},Conj{S}}, r2::RationalFunction{Var{T},Conj{S}})
  num = r1.num * r2.den
  den = r1.den * r2.num
  RationalFunction(num, den, Conj{S})
end

function /{T1,S1,T2,S2}(r1::RationalFunction{Var{T1},Conj{S1}}, r2::RationalFunction{Var{T2},Conj{S2}})
  warn("r1/r2: `r1` ($T1,Conj{$S1}) and `r2` ($T2,Conj{$S2}) have different variables")
  throw(DomainError())
end

-{T,S}(r::RationalFunction{T,S})                = RationalFunction(-r.num, copy(r.den), S)

-(r1::RationalFunction, r2::RationalFunction)   = +(r1, -r2)

.+(r1::RationalFunction, r2::RationalFunction)  = +(r1, r2)
.*(r1::RationalFunction, r2::RationalFunction)  = *(r1, r2)
./(r1::RationalFunction, r2::RationalFunction)  = /(r1, r2)
.-(r1::RationalFunction, r2::RationalFunction)  = -(r1, r2)

## Basic operations between `Number`s
==(r::RationalFunction, n::Number)              = ==(r.num, n*r.den)
isapprox{T,S,U,V,Z<:Number}(r::RationalFunction{T,S,U,V}, n::Z;
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) =
  isapprox(coeffs(r.num), coeffs(n*r.den); rtol = rtol, atol = atol)

+{T,S}(r::RationalFunction{T,S}, n::Number) = RationalFunction(r.num + n*r.den, copy(r.den), S)
+(n::Number, r::RationalFunction)           = +(r, n)
*{T,S}(r::RationalFunction{T,S}, n::Number) = RationalFunction(n*r.num, copy(r.den), S)
*(n::Number, r::RationalFunction)           = *(r, n)
dot(r::RationalFunction, n::Number)         = *(r, n)
dot(n::Number, r::RationalFunction)         = *(r, n)
/{T,S}(r::RationalFunction{T,S}, n::Number) = RationalFunction(copy(r.num), n*r.den, S)
/{T,S}(n::Number, r::RationalFunction{T,S}) = RationalFunction(n*r.den, copy(r.num), S)
-(r::RationalFunction, n::Number)           = +(r, -n)
-(n::Number, r::RationalFunction)           = +(-r, n)

.+(r::RationalFunction, n::Number)  = +(r, n)
.*(r::RationalFunction, n::Number)  = *(r, n)
./(r::RationalFunction, n::Number)  = /(r, n)
.-(r::RationalFunction, n::Number)  = +(r, -n)

.+(n::Number, r::RationalFunction)  = +(r, n)
.*(n::Number, r::RationalFunction)  = *(r, n)
./(n::Number, r::RationalFunction)  = /(n, r)
.-(n::Number, r::RationalFunction)  = +(-r, n)

## Basic operations between `Poly`s
function =={T,S}(r::RationalFunction{Var{T},S}, p::Poly)
  T ≠ p.var && return false
  return r.num == p*r.den
end

function isapprox{T,S,U,V,Z<:Number}(r::RationalFunction{Var{T},S,U,V}, p::Poly{Z};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0)
  if T ≠ p.var
    warn("r≈p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  isapprox(coeffs(r.num), coeffs(p*r.den); rtol = rtol, atol = atol)
end

function +{T,S}(r::RationalFunction{Var{T},S}, p::Poly)
  if T ≠ p.var
    warn("r+p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(r.num + p*r.den, r.den, S)
end
+(p::Poly, r::RationalFunction) = +(r, p)

function *{T,S}(r::RationalFunction{Var{T},S}, p::Poly)
  if T ≠ p.var
    warn("r*p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(p*r.num, r.den, S)
end
*(p::Poly, r::RationalFunction)           = *(r, p)
dot(r::RationalFunction, p::Poly)         = *(r, p)
dot(p::Poly, r::RationalFunction)         = *(r, p)

function /{T,S}(r::RationalFunction{Var{T},S}, p::Poly)
  if T ≠ p.var
    warn("r/p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(r.num, p*r.den, S)
end

function /{T,S}(p::Poly, r::RationalFunction{Var{T},S})
  if T ≠ p.var
    warn("p/r: `p` ($(p.var)) and `r` ($T) have different variables")
    throw(DomainError())
  end
  RationalFunction(p*r.den, r.num, S)
end

-(r::RationalFunction, p::Poly) = +(r, -p)
-(p::Poly, r::RationalFunction) = +(-r, p)

.+(r::RationalFunction, p::Poly)  = +(r, p)
.*(r::RationalFunction, p::Poly)  = *(r, p)
./(r::RationalFunction, p::Poly)  = /(r, p)
.-(r::RationalFunction, p::Poly)  = +(r, -p)

.+(p::Poly, r::RationalFunction)  = +(r, p)
.*(p::Poly, r::RationalFunction)  = *(r, p)
./(p::Poly, r::RationalFunction)  = /(p, r)
.-(p::Poly, r::RationalFunction)  = +(-r, p)

## Solve
function solve(lhs::RationalFunction, rhs::Number = 0)
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  rhs == zero(rhs)  && return zeros(lhs)
  isinf(rhs)        && return poles(lhs)
  zeros(lhs - rhs)
end

function solve(lhs::RationalFunction, rhs::Poly)
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  rhs == zero(rhs)  && return zeros(lhs)
  zeros(lhs - rhs)
end

solve(lhs::PolyLike, rhs::RationalFunction) = solve(rhs, lhs)

function solve{T,S}(lhs::RationalFunction{Var{T},Conj{S}},
  rhs::RationalFunction{Var{T},Conj{S}})
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  zeros(lhs - rhs)
end

function solve{T1,S1,T2,S2}(lhs::RationalFunction{Var{T1},Conj{S1}},
  rhs::RationalFunction{Var{T2},Conj{S2}})
  warn("solve(lhs,rhs): `lhs` ($T1,Conj{$S1}) and `rhs` ($T2,Conj{$S2}) have different variables")
  throw(DomainError())
end
