# Copy operation
copy{T,S}(r::RationalFunction{T,S}) = RationalFunction(copy(r.num), copy(r.den), S)

# Iteration interface
eltype{T,S,U,V}(::Type{RationalFunction{T,S,U,V}}) = (U, V)
eltype{T,S,U,V}(r::RationalFunction{T,S,U,V})      = (U, V)

# Convenience functions
"""
    coeffs(r::RationalFunction) -> Tuple{Vector,Vector}

Return the coefficient vectors of the numerator and denominator polynomials of `r`.

See also: `num`, `den`, `degree`.
"""
coeffs(r::RationalFunction) = (coeffs(r.num), coeffs(r.den))

"""
    degree(r::RationalFunction) -> Tuple{Int,Int}

Return the degrees of the numerator and denominator polynomials of `r`.

See also: `num`, `den`, `coeffs`.
"""
degree(r::RationalFunction) = (degree(r.num), degree(r.den))

"""
    roots(r::RationalFunction) -> Tuple{Vector,Vector}

Return the roots of the numerator and denominator polynomials of `r`.

See also: `num`, `den`, `zeros`, `poles`, `solve`.
"""
roots(r::RationalFunction)  = (roots(r.num), roots(r.den))

"""
    variable(r::RationalFunction) -> Tuple{Poly,Poly,Val}

Return the variables of the numerator and denominator polynomials of `r` as well
as the conjugation property.
"""
variable{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}}) =
  (variable(U, T), variable(V, T), Val{S})
variable{T,S,U,V}(r::RationalFunction{Val{T},Val{S},U,V})      =
  (variable(U, T), variable(V, T), Val{S})

"""
    num(r::RationalFunction) -> Poly

Return the numerator polynomial of `r`.

See also: `den`.
"""
num(r::RationalFunction)    = r.num

"""
    den(r::RationalFunction) -> Poly

Return the denominator polynomial of `r`.

See also: `num`.
"""
den(r::RationalFunction)    = r.den

"""
    zeros(r::RationalFunction) -> Vector

Return the values which make `r` zero.

See also: `roots`, `poles`, `solve`.
"""
zeros(r::RationalFunction)  = (rnew = reduce(r); roots(rnew.num))

"""
    poles(r::RationalFunction) -> Vector

Return the values which make `r` unbounded, *i.e.*, tend to infinity.

See also: `roots`, `zeros`, `solve`.
"""
poles(r::RationalFunction)  = (rnew = reduce(r); roots(rnew.den))

## Identities
one{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}})  =
  RationalFunction(Poly([one(U)], T), Poly([one(V)], T), Val{S})
one{T,S}(r::RationalFunction{T,S})                          =
  RationalFunction(one(r.num), one(r.den), S)
zero{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}}) =
  RationalFunction(Poly(U[], T), Poly([one(V)], T), Val{S})
zero{T,S}(r::RationalFunction{T,S})                         =
  RationalFunction(zero(r.num), one(r.den), S)

## Comparison
hash{T,S}(r::RationalFunction{Val{T},Val{S}}, h::UInt)             =
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

function isapprox{T,S,U1,V1,U2,V2}(r1::RationalFunction{Val{T},Val{S},U1,V1},
  r2::RationalFunction{Val{T},Val{S},U2,V2};
  rtol::Real = sqrt(eps(promote_type(U1,V2,U2,V2))), atol::Real = 0)
  p1 = r1.num * r2.den
  p2 = r1.den * r2.num

  isapprox(coeffs(p1), coeffs(p2); rtol = rtol, atol = atol)
end

function isapprox{T1,S1,T2,S2}(r1::RationalFunction{Val{T1},Val{S1}},
  r2::RationalFunction{Val{T2},Val{S2}}; rtol::Real = 0, atol::Real = 0)
  warn("r1≈r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables")
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

_funceval(r::RationalFunction, X) = map(r, X)

@compat (r::RationalFunction)(x::Number)                    = _funceval(r, x)
@compat (r::RationalFunction{T,Val{:conj}}){T}(x::Number)   = _funceval(r, conj(x))
@compat (r::RationalFunction{T,Val{:conj}}){T}(x::Real)     = _funceval(r, x)
@compat (r::RationalFunction)(X)                            = _funceval(r, X)

"""
    funcfit(x, y, m::Int, n::Int = 0, var = :x) -> RationalFunction

Fit a rational function of variable `var`, which
  * has a numerator and a denominator of degrees `m` and `n`, respectively, and,
  * which will output the values in `y` to the corresponding inputs in `x`.

When `m` is also dropped, *i.e.*,

    funcfit(x, y, var = :x) -> RationalFunction

fit a rational function of variable `var`, which
  * has a numerator and a denominator of degrees `m` and `n`, respectively, and,
  * which will output the values in `y` to the corresponding inputs in `x`,
where `n, r = divrem(length(x), 2)` and `m = n - 1 + r`.
"""
function funcfit(x, y, m::Int, n::Int, var::SymbolLike = :x)
  if length(x) ≠ length(y)
    warn("funcfit(x, y, m, n, var): length(x) ≠ length(y)")
    throw(DomainError())
  elseif m < 0 || n < 0
    warn("funcfit(x, y, m, n, var): `m` and `n` must be non-negative")
    throw(DomainError())
  elseif m + n + 1 > length(x)
    warn("funcfit(x, y, m, n, var): not enough data points for given degrees")
    throw(DomainError())
  end

  T   = promote_type(eltype(x), eltype(y))

  A1  = T[ x[i]^j        for i in 1:length(x), j in 0:m]
  A2  = T[-x[i]^j * y[i] for i in 1:length(x), j in 1:n]
  A   = [A1 A2]

  b   = A \ y

  RationalFunction(b[1:m+1], [1; b[m+2:end]], var)
end

funcfit(x, y, m::Int, var::SymbolLike = :x) = funcfit(x, y, m, 0, var)
function funcfit(x, y, var::SymbolLike = :x)
  n, r  = divrem(length(x), 2)
  m     = n - 1 + r
  funcfit(x, y, m, n, var)
end

# Mathematical operations (always return temporaries for correctness of the results)
## Inversion
inv{T,S}(r::RationalFunction{T,S})  = RationalFunction(copy(r.den), copy(r.num), S)

## Transposition
transpose(r::RationalFunction)      = copy(r)

## Conjugation
function conj{T,S}(r::RationalFunction{Val{T},Val{S}})
  numcoeff, dencoeff = coeffs(r)
  RationalFunction(Poly(conj(copy(numcoeff)), T), Poly(conj(copy(dencoeff)), T),
    Val{ifelse(S == :conj, :notc, :conj)})
end
# Related to #2. This is the solution in Julia v0.6 in `arraymath.jl`
conj{T,S,U,V}(m::AbstractArray{RationalFunction{Val{T},Val{S},U,V}}) = map(conj, m)

## Derivative
"""
    derivative(r::RationalFunction, n::Int = 1) -> RationalFunction

Take the `n`th order derivative of `r`.
"""
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
"""
    reduce(r::RationalFunction) -> RationalFunction

Do the pole-zero cancellation, *i.e.*, remove the common roots of the numerator
and denominator polynomials of `r`.

See also: `zeros`, `poles`, `roots`.
"""
function reduce{T,S}(r::RationalFunction{T,S})
  g       = gcd(r.num, r.den)
  common  = degree(g) == 0 ? one(g) : g

  num, _  = divrem(r.num, common)
  den, _  = divrem(r.den, common)
  RationalFunction(num, den, S)
end

## Basic operations between rational functions
function +{T,S}(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}})
  g       = gcd(r1.den, r2.den)
  common  = Polynomials.degree(g) == 0 ? one(g) : g

  den1, _ = divrem(r1.den, common)
  den2, _ = divrem(r2.den, common)
  num     = r1.num * den2 + den1 * r2.num
  den     = den1*den2*common
  RationalFunction(num, den, Val{S})
end

function +{T1,S1,T2,S2}(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}})
  warn("r1+r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables")
  throw(DomainError())
end

function *{T,S}(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}})
  num = r1.num * r2.num
  den = r1.den * r2.den
  RationalFunction(num, den, Val{S})
end

function *{T1,S1,T2,S2}(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}})
  warn("r1*r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables")
  throw(DomainError())
end

dot(r1::RationalFunction, r2::RationalFunction) = *(r1, r2)

function /{T,S}(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}})
  num = r1.num * r2.den
  den = r1.den * r2.num
  RationalFunction(num, den, Val{S})
end

function /{T1,S1,T2,S2}(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}})
  warn("r1/r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables")
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
==(n::Number, r::RationalFunction)              = ==(r, n)
isapprox{T,S,U,V,Z<:Number}(r::RationalFunction{T,S,U,V}, n::Z;
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) =
  isapprox(coeffs(r.num), coeffs(n*r.den); rtol = rtol, atol = atol)
isapprox{T,S,U,V,Z<:Number}(n::Z, r::RationalFunction{T,S,U,V};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) =
  isapprox(r, n; rtol = rtol, atol = atol)

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
function =={T,S}(r::RationalFunction{Val{T},S}, p::Poly)
  T ≠ p.var && return false
  return r.num == p*r.den
end
==(p::Poly, r::RationalFunction) = ==(r, p)

function isapprox{T,S,U,V,Z<:Number}(r::RationalFunction{Val{T},S,U,V}, p::Poly{Z};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0)
  if T ≠ p.var
    warn("r≈p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  isapprox(coeffs(r.num), coeffs(p*r.den); rtol = rtol, atol = atol)
end
isapprox{T,S,U,V,Z<:Number}(p::Poly{Z}, r::RationalFunction{Val{T},S,U,V};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) =
  isapprox(r, p; rtol = rtol, atol = atol)

function +{T,S}(r::RationalFunction{Val{T},S}, p::Poly)
  if T ≠ p.var
    warn("r+p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(r.num + p*r.den, r.den, S)
end
+(p::Poly, r::RationalFunction) = +(r, p)

function *{T,S}(r::RationalFunction{Val{T},S}, p::Poly)
  if T ≠ p.var
    warn("r*p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(p*r.num, r.den, S)
end
*(p::Poly, r::RationalFunction)           = *(r, p)
dot(r::RationalFunction, p::Poly)         = *(r, p)
dot(p::Poly, r::RationalFunction)         = *(r, p)

function /{T,S}(r::RationalFunction{Val{T},S}, p::Poly)
  if T ≠ p.var
    warn("r/p: `r` ($T) and `p` ($(p.var)) have different variables")
    throw(DomainError())
  end
  RationalFunction(r.num, p*r.den, S)
end

function /{T,S}(p::Poly, r::RationalFunction{Val{T},S})
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
"""
    solve(lhs, rhs = 0) -> Vector

Solve for the values which make `lhs = rhs`, where at least one of `lhs` and `rhs`
is a `RationalFunction` while the other can be either a `Number` or a `Poly`.

See also: `zeros`, `poles`, `roots`.
"""
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

function solve{T,S}(lhs::RationalFunction{Val{T},Val{S}},
  rhs::RationalFunction{Val{T},Val{S}})
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  zeros(lhs - rhs)
end

function solve{T1,S1,T2,S2}(lhs::RationalFunction{Val{T1},Val{S1}},
  rhs::RationalFunction{Val{T2},Val{S2}})
  warn("solve(lhs,rhs): `lhs` ($T1,Val{$S1}) and `rhs` ($T2,Val{$S2}) have different variables")
  throw(DomainError())
end
