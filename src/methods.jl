# Copy operation
copy(r::RationalFunction{T,S}) where {T,S} = RationalFunction(copy(r.num), copy(r.den), S)

# Iteration interface
eltype(::Type{RationalFunction{T,S,U,V}}) where {T,S,U,V} = (U, V)
eltype(r::RationalFunction{T,S,U,V}) where {T,S,U,V} = (U, V)

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
    variable(r::RationalFunction) -> Tuple{Polynomial,Polynomial,Val}

Return the variables of the numerator and denominator polynomials of `r` as well
as the conjugation property.
"""
variable(::Type{RationalFunction{Val{T},Val{S},U,V}}) where {T,S,U,V} =
  (variable(Polynomial{U}, T), variable(Polynomial{V}, T), Val{S})
variable(r::RationalFunction{Val{T},Val{S},U,V}) where {T,S,U,V} =
  (variable(Polynomial{U}, T), variable(Polynomial{V}, T), Val{S})

"""
    numerator(r::RationalFunction) -> Polynomial

Return the numerator polynomial of `r`.

See also: `den`.
"""
numerator(r::RationalFunction)    = r.num

"""
    denominator(r::RationalFunction) -> Polynomial

Return the denominator polynomial of `r`.

See also: `num`.
"""
denominator(r::RationalFunction)    = r.den

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
one(::Type{RationalFunction{Val{T},Val{S},U,V}}) where {T,S,U,V}  =
  RationalFunction(Polynomial([one(U)], T), Polynomial([one(V)], T), Val{S})
one(r::RationalFunction{T,S}) where {T,S}                         =
  RationalFunction(one(r.num), one(r.den), S)
zero(::Type{RationalFunction{Val{T},Val{S},U,V}}) where {T,S,U,V} =
  RationalFunction(Polynomial(U[], T), Polynomial([one(V)], T), Val{S})
zero(r::RationalFunction{T,S}) where {T,S}                        =
  RationalFunction(zero(r.num), one(r.den), S)

## Comparison
hash(r::RationalFunction{Val{T},Val{S}}, h::UInt) where {T,S}       =
  hash(T, hash(S, hash(coeffs(r.num), hash(coeffs(r.den), h))))

==(r1::RationalFunction{T,S}, r2::RationalFunction{T,S}) where {T,S} =
  r1.num * r2.den == r1.den * r2.num
==(r1::RationalFunction, r2::RationalFunction)                      = false

isequal(r1::RationalFunction{T,S}, r2::RationalFunction{T,S}) where {T,S} =
  hash(r1) == hash(r2)
isequal(r1::RationalFunction, r2::RationalFunction)                 = false

# Auxiliary definitions for `isapprox`
eps(::Type{T}) where {T<:AbstractFloat}          = Base.eps(T)
eps(::Type{Complex{T}}) where {T<:AbstractFloat} = Base.eps(T)
eps(::Type{T}) where {T}                         = zero(T)

function isapprox(r1::RationalFunction{Val{T},Val{S},U1,V1},
  r2::RationalFunction{Val{T},Val{S},U2,V2};
  rtol::Real = sqrt(eps(promote_type(U1,V2,U2,V2))), atol::Real = 0) where {T,S,U1,V1,U2,V2}
  p1 = r1.num * r2.den
  p2 = r1.den * r2.num

  isapprox(coeffs(p1), coeffs(p2); rtol = rtol, atol = atol)
end

function isapprox(r1::RationalFunction{Val{T1},Val{S1}},
  r2::RationalFunction{Val{T2},Val{S2}}; rtol::Real = 0, atol::Real = 0) where {T1,S1,T2,S2}
  throw(DomainError((r1,r2), "r1≈r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables"))
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
    num     = derivative(num)
    den     = derivative(den)
    result  = num(x)/den(x)
    k      += 1
  end
  result
end

_funceval(r::RationalFunction, X) = map(r, X)

@compat (r::RationalFunction)(x::Number)                         = _funceval(r, x)
@compat (r::RationalFunction{T,Val{:conj}})(x::Number) where {T} = _funceval(r, conj(x))
@compat (r::RationalFunction{T,Val{:conj}})(x::Real) where {T}   = _funceval(r, x)
@compat (r::RationalFunction)(X)                                 = _funceval(r, X)

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
    throw(DomainError((x,y,m,n,var), "funcfit(x, y, m, n, var): length(x) ≠ length(y)"))
  elseif m < 0 || n < 0
    throw(DomainError((x,y,m,n,var), "funcfit(x, y, m, n, var): `m` and `n` must be non-negative"))
  elseif m + n + 1 > length(x)
    throw(DomainError((x,y,m,n,var), "funcfit(x, y, m, n, var): not enough data points for given degrees"))
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
inv(r::RationalFunction{T,S}) where {T,S} = RationalFunction(copy(r.den), copy(r.num), S)

## Transposition
transpose(r::RationalFunction)      = copy(r)

## Conjugation
function conj(r::RationalFunction{Val{T},Val{S}}) where {T,S}
  numcoeff, dencoeff = coeffs(r)
  RationalFunction(Polynomial(conj(copy(numcoeff)), T), Polynomial(conj(copy(dencoeff)), T),
    Val{ifelse(S == :conj, :notc, :conj)})
end
# Related to #2. This is the solution in Julia v0.6 in `arraymath.jl`
conj(m::AbstractArray{RationalFunction{Val{T},Val{S},U,V}}) where {T,S,U,V} = map(conj, m)

## Derivative
"""
    derivative(r::RationalFunction, n::Int = 1) -> RationalFunction

Take the `n`th order derivative of `r`.
"""
function derivative(r::RationalFunction{T,S}, n::Int = 1) where {T,S}
  if n < 0
    throw(DomainError((r, n),"derivative(r, n): `n` must be non-negative" ))
  end
  n == 0 && return copy(r)
  num   = derivative(r.num)*r.den - r.num*derivative(r.den)
  den   = r.den*r.den
  temp  = RationalFunction(num, den, S)
  for count in 2:n
    num   = derivative(temp.num)*temp.den - temp.num*derivative(temp.den)
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
function reduce(r::RationalFunction{T,S}) where {T,S}
  g       = gcd(r.num, r.den)
  common  = degree(g) == 0 ? one(g) : g

  num, _  = divrem(r.num, common)
  den, _  = divrem(r.den, common)
  RationalFunction(num, den, S)
end

## Basic operations between rational functions
function +(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}}) where {T,S}
  g       = gcd(r1.den, r2.den)
  common  = Polynomials.degree(g) == 0 ? one(g) : g

  den1, _ = divrem(r1.den, common)
  den2, _ = divrem(r2.den, common)
  num     = r1.num * den2 + den1 * r2.num
  den     = den1*den2*common
  RationalFunction(num, den, Val{S})
end

function +(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}}) where {T1,S1,T2,S2}
  throw(DomainError((r1, r2), "r1+r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables"))
end

function *(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}}) where {T,S}
  num = r1.num * r2.num
  den = r1.den * r2.den
  RationalFunction(num, den, Val{S})
end

function *(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}}) where {T1,S1,T2,S2}
  throw(DomainError((r1, r2), "r1*r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables"))
end

dot(r1::RationalFunction, r2::RationalFunction) = *(r1, r2)

^(r::RationalFunction, n::Integer) = Base.power_by_squaring(r,n)

function /(r1::RationalFunction{Val{T},Val{S}}, r2::RationalFunction{Val{T},Val{S}}) where {T,S}
  num = r1.num * r2.den
  den = r1.den * r2.num
  RationalFunction(num, den, Val{S})
end

function /(r1::RationalFunction{Val{T1},Val{S1}}, r2::RationalFunction{Val{T2},Val{S2}}) where {T1,S1,T2,S2}
  throw(DomainError((r1, r2), "r1/r2: `r1` ($T1,Val{$S1}) and `r2` ($T2,Val{$S2}) have different variables"))
end

-(r::RationalFunction{T,S}) where {T,S}         = RationalFunction(-r.num, copy(r.den), S)

-(r1::RationalFunction, r2::RationalFunction)   = +(r1, -r2)

## Basic operations between `Number`s
==(r::RationalFunction, n::Number)              = ==(r.num, n*r.den)
==(n::Number, r::RationalFunction)              = ==(r, n)
isapprox(r::RationalFunction{T,S,U,V}, n::Z;
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) where {T,S,U,V,Z<:Number} =
  isapprox(coeffs(r.num), coeffs(n*r.den); rtol = rtol, atol = atol)
isapprox(n::Z, r::RationalFunction{T,S,U,V};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) where {T,S,U,V,Z<:Number} =
  isapprox(r, n; rtol = rtol, atol = atol)

+(r::RationalFunction{T,S}, n::Number) where {T,S} = RationalFunction(r.num + n*r.den, copy(r.den), S)
+(n::Number, r::RationalFunction)                  = +(r, n)
*(r::RationalFunction{T,S}, n::Number) where {T,S} = RationalFunction(n*r.num, copy(r.den), S)
*(n::Number, r::RationalFunction)                  = *(r, n)
dot(r::RationalFunction, n::Number)                = *(r, n)
dot(n::Number, r::RationalFunction)                = *(r, n)
/(r::RationalFunction{T,S}, n::Number) where {T,S} = RationalFunction(copy(r.num), n*r.den, S)
/(n::Number, r::RationalFunction{T,S}) where {T,S} = RationalFunction(n*r.den, copy(r.num), S)
-(r::RationalFunction, n::Number)                  = +(r, -n)
-(n::Number, r::RationalFunction)                  = +(-r, n)

## Basic operations between `Polynomial`s
function ==(r::RationalFunction{Val{T},S}, p::Polynomial) where {T,S}
  T ≠ p.var && return false
  return r.num == p*r.den
end
==(p::Polynomial, r::RationalFunction) = ==(r, p)

function isapprox(r::RationalFunction{Val{T},S,U,V}, p::Polynomial{Z};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) where {T,S,U,V,Z<:Number}
  if T ≠ p.var
    throw(DomainError((r,p), "r≈p: `r` ($T) and `p` ($(p.var)) have different variables"))
  end
  isapprox(coeffs(r.num), coeffs(p*r.den); rtol = rtol, atol = atol)
end
isapprox(p::Polynomial{Z}, r::RationalFunction{Val{T},S,U,V};
  rtol::Real = sqrt(eps(promote_type(U,V,Z))), atol::Real = 0) where {T,S,U,V,Z<:Number} =
  isapprox(r, p; rtol = rtol, atol = atol)

function +(r::RationalFunction{Val{T},S}, p::Polynomial) where {T,S}
  if T ≠ p.var
    throw(DomainError((r,p), "r+p: `r` ($T) and `p` ($(p.var)) have different variables"))
  end
  RationalFunction(r.num + p*r.den, r.den, S)
end
+(p::Polynomial, r::RationalFunction) = +(r, p)

function *(r::RationalFunction{Val{T},S}, p::Polynomial) where {T,S}
  if T ≠ p.var
    throw(DomainError((r,p), "r*p: `r` ($T) and `p` ($(p.var)) have different variables"))
  end
  RationalFunction(p*r.num, r.den, S)
end
*(p::Polynomial, r::RationalFunction)           = *(r, p)
dot(r::RationalFunction, p::Polynomial)         = *(r, p)
dot(p::Polynomial, r::RationalFunction)         = *(r, p)

function /(r::RationalFunction{Val{T},S}, p::Polynomial) where {T,S}
  if T ≠ p.var
    throw(DomainError((r,p), "r/p: `r` ($T) and `p` ($(p.var)) have different variables"))
  end
  RationalFunction(r.num, p*r.den, S)
end

function /(p::Polynomial, r::RationalFunction{Val{T},S}) where {T,S}
  if T ≠ p.var
    throw(DomainError((r,p), "p/r: `p` ($(p.var)) and `r` ($T) have different variables"))
  end
  RationalFunction(p*r.den, r.num, S)
end

-(r::RationalFunction, p::Polynomial) = +(r, -p)
-(p::Polynomial, r::RationalFunction) = +(-r, p)

## Solve
"""
    solve(lhs, rhs = 0) -> Vector

Solve for the values which make `lhs = rhs`, where at least one of `lhs` and `rhs`
is a `RationalFunction` while the other can be either a `Number` or a `Polynomial`.

See also: `zeros`, `poles`, `roots`.
"""
function solve(lhs::RationalFunction, rhs::Number = 0)
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  rhs == zero(rhs)  && return zeros(lhs)
  isinf(rhs)        && return poles(lhs)
  zeros(lhs - rhs)
end

function solve(lhs::RationalFunction, rhs::Polynomial)
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  rhs == zero(rhs)  && return zeros(lhs)
  zeros(lhs - rhs)
end

solve(lhs::PolynomialLike, rhs::RationalFunction) = solve(rhs, lhs)

function solve(lhs::RationalFunction{Val{T},Val{S}},
  rhs::RationalFunction{Val{T},Val{S}}) where {T,S}
  lhs == rhs        && error("solve(lhs,rhs): `lhs` and `rhs` are equal")
  zeros(lhs - rhs)
end

function solve(lhs::RationalFunction{Val{T1},Val{S1}},
  rhs::RationalFunction{Val{T2},Val{S2}}) where {T1,S1,T2,S2}
  throw(DomainError((lhs,rhs), "solve(lhs,rhs): `lhs` ($T1,Val{$S1}) and `rhs` ($T2,Val{$S2}) have different variables"))
end

### Given rational function, return partial fraction decomposition
"""
    r, p, k = residue(num, den)
    r, p, k = residue(x::RationalFunction)

Given a numerator and denominator of a polynomial fraction or a RationalFunction,
return the partial fraction decomposition. The residues are returned as `r`, the
poles as `p`, and the direct terms as `k`. Both numerator and denominator can be
given as an array of coefficients or as Polynomial types.

Duplicate poles are handled by each subsequent duplicate having the
denominator's power raised by one. See package README for further explanation.

# Examples
``` julia
# Real poles with direct terms
julia> num1 = Polynomial([6, 9, 16, 8, 1])
julia> den1 = Polynomial([6, 11, 6, 1])
julia> residue(num1, den1)
([-6.0, -4.0, 3.0], [-3.0, -2.0, -1.0], [2.0, 1.0])

# Complex poles
julia> residue([10, 2], [0, 10, 2, 1])
(Complex{Float64}[-0.5-0.166667im, -0.5+0.166667im, 1.0+0.0im], Complex{Float64}[-1.0+3.0im, -1.0-3.0im, 0.0+0.0im], [0.0])

# Duplicate poles
julia> r_func = RationalFunction(Polynomial([1,0,1]),Polynomial([0,1])*Polynomial([-1,1])^2))
julia> residue(r_func)
([-0.0, 2.0, 1.0], [1.0, 1.0, 0.0], [0.0])
```
"""
function residue(num::Polynomial, den::Polynomial)
    p = roots(den)
    num_p = length(p)
    div_poly, rem_poly = divrem(num, den)
    k = coeffs(div_poly)
    rem_coeffs = coeffs(rem_poly)
    # Account for repeated poles. The dictionary created is used for pole powers
    unique_p = unique(p)
    dup_p_dict = Dict()
    for pole in unique_p
        temp_count = count((i)->i==pole, p)
        if temp_count > 1
          dup_p_dict[pole] = temp_count-1
        end
    end
    # Capture case in which the poles create terms with greater
    # powers than those in the numerator
    if isempty(dup_p_dict)
        max_power = num_p
    else
        max_power = max(num_p, maximum(values(dup_p_dict))+1)
    end
    if length(rem_coeffs)<max_power
        temp_rem_coeffs = zeros(max_power)
        temp_rem_coeffs[1:length(rem_coeffs)] = rem_coeffs
        rem_coeffs = temp_rem_coeffs
    end
    # Construct matrix representing the system of equations for
    # the residues
    residue_mtx = zeros(eltype(p),max_power, num_p)
    for col_idx in 1:num_p
        temp_p = p[col_idx]
        if temp_p in keys(dup_p_dict)
            temp_poly = Polynomial([-temp_p, 1])^dup_p_dict[temp_p]
            dup_p_dict[temp_p]-=1
        else
            temp_poly = Polynomial(1)
        end
        for poly_idx in 1:num_p
            if p[poly_idx] != temp_p
                temp_poly *= Polynomial([-p[poly_idx], 1])
            end
        end
        temp_coeffs = coeffs(temp_poly)
        zero_pad = zeros(eltype(p),1,max_power)
        zero_pad[1:length(temp_coeffs)] = temp_coeffs
        residue_mtx[:,col_idx] = zero_pad
    end
    r = residue_mtx\rem_coeffs
    return r, p, k
end
residue(num::Vector{T}, den::Vector{T}) where {T<:Number} = residue(Polynomial(num), Polynomial(den))
residue(rfunc::RationalFunction) = residue(numerator(rfunc), denominator(rfunc))

### Conversely, given r, p, and k terms, return polynomial fraction
"""
========================================

    num, dem = residue(r, p, k)

Given a set of residues, poles, and direct terms, calculate the corresponding
polynomial fraction.

# Examples
``` julia
# Real poles with direct terms
julia> r = [-6.0, -4.0, 3.0]
julia> p = [-3.0, -2.0, -1.0]
julia> k = [2.0, 1.0]
julia> residue(r, p, k)
([6, 9, 16, 8, 1], [6, 11, 6, 1])

# Complex poles
julia> r = [-0.5-0.166667im, -0.5+0.166667im, 1.0+0.0im]
julia> p = [-1.0+3.0im, -1.0-3.0im, 0.0+0.0im]
julia> k = [0.0]
julia> residue(r, p, k)
(Complex{Float64}[10.0+0.0im, 2.0+0.0im], Complex{Float64}[0.0+0.0im, 10.0+0.0im, 2.0+0.0im, 1.0+0.0im])

# Duplicate poles
julia> r = [-0.0, 2.0, 1.0]
julia> p = [1.0, 1.0, 0.0]
julia> k = [0.0]
julia> residue(r, p, k)
([1.0, 0.0, 1.0], [0.0, 1.0, -2.0, 1.0])
```
"""
function residue(r::Vector{T}, p::Vector{S}, k::Vector{U}) where {T<:Number,S<:Number,U<:Number}
    den_poly = fromroots(p)
    unique_p = unique(p)
    dup_p_dict = Dict()
    p_count_dict = Dict()
    for pole in unique_p
        temp_count = count((i)->i==pole, p)
        p_count_dict[pole] = temp_count
        if temp_count > 1
          dup_p_dict[pole] = temp_count-1
        end
    end
    r_poly = Polynomial(0)
    for resid_idx in 1:length(r)
        temp_num_term = Polynomial(1)
        if p[resid_idx] in keys(dup_p_dict)
            temp_num_term *= Polynomial([-p[resid_idx],1])^dup_p_dict[p[resid_idx]]
            dup_p_dict[p[resid_idx]]-=1
        end
        for pole in keys(p_count_dict)
            if pole != p[resid_idx]
                temp_num_term *= Polynomial([-pole, 1])^p_count_dict[pole]
            end
        end
        r_poly += temp_num_term*r[resid_idx]
    end
    k_poly = Polynomial(k)*den_poly
    num_poly = k_poly+r_poly
    return coeffs(num_poly), coeffs(den_poly)
end
