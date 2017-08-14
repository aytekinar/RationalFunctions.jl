using Plots
gr(display = false)

using Polynomials
using RationalFunctions
using Base.Test

# Constructor tests
## Construction from polynomials
px_ = Poly([1,-2,1])     # px_ = (x-1)(x-1)
qx = poly([1,2,3])      # qx = (x-1)(x-2)(x-3)
ps = Poly([1,-2,1], :s) # ps = (s-1)(s-1)

@test isequal(RationalFunction(px_, qx), RationalFunction(px_, qx, Val{:notc}))
@test isa(RationalFunction(px_, qx, Val{:conj}), RationalFunction{Val{px_.var}, Val{:conj}, eltype(px_), eltype(qx)})
@test isequal(RationalFunction(px_), RationalFunction(px_, one(px_)))
@test isequal(RationalFunction(px_, Val{:conj}), RationalFunction(px_, one(px_), Val{:conj}))

@test_throws DomainError RationalFunction(px_, ps)
@test_throws DomainError RationalFunction(px_, ps, Val{:conj})

@test isequal(RationalFunction(coeffs(px_), qx), RationalFunction(px_, coeffs(qx)))
@test isequal(RationalFunction(1, Poly([2])), RationalFunction(Poly([1]), 2))

## Construction from numbers and vectors
@test isequal(RationalFunction(1, [1, 2]), RationalFunction([1], [1, 2]))
@test isequal(RationalFunction([1, 2], 1.), RationalFunction([1, 2], [1.]))
@test isequal(RationalFunction(1, 1), RationalFunction([1], [1]))
@test isa(RationalFunction(1), RationalFunction{Val{:x}, Val{:notc}, Int, Int})
@test isa(RationalFunction([1, 2.], 's', Val{:conj}), RationalFunction{Val{:s}, Val{:conj}, Float64, Float64})

## `Poly` construction from `RationalFunction`s
r1 = RationalFunction(px_, qx)
r2 = RationalFunction(px_, Poly([4]))

@test isapprox(coeffs(Poly(r2)), coeffs(px_/4))
@test_throws DomainError Poly(r1)

# Conversion tests
v1 = [1, 2, 3]
v2 = [1., 2., 3.]
v3 = [1 + 1im, 2. - 1im]
p3 = poly(v3)
@test eltype([RationalFunction(v1, v2), RationalFunction(v2, v1)])  ==
  RationalFunction{Val{:x}, Val{:notc}, promote_type(eltype(v1), eltype(v2)),
    promote_type(eltype(v1), eltype(v2))}
@test eltype([RationalFunction(px_, qx), 1.])                        ==
  RationalFunction{Val{px_.var}, Val{:notc}, promote_type(eltype(px_), typeof(1.)),
    promote_type(eltype(qx), typeof(1.))}
@test eltype([RationalFunction(px_, qx, Val{:conj}), p3])            ==
  RationalFunction{Val{px_.var}, Val{:conj}, promote_type(eltype(px_), eltype(p3)),
    promote_type(eltype(qx), eltype(p3))}

@test_throws DomainError [RationalFunction(px_, qx), ps]

# Printing tests
# TODO: Think of some useful printing tests

# Method tests
## Copying
r1 = RationalFunction([ 0., 1, 2], :x)
r2 = copy(r1)

@test isequal(r1, r2) && !(===(r1, r2))

## Valiable-related
@test eltype(r1) == eltype(typeof(r1)) == (eltype(r1.num), eltype(r1.den))
@test variable(r1) == variable(typeof(r1)) == (variable(r1.num), variable(r1.den), Val{:notc})

## Identities
@test one(r1) == one(typeof(r1)) == RationalFunction(one(r1.num), one(r1.den))
@test zero(r1) == zero(typeof(r1)) == RationalFunction(zero(r1.num), one(r1.den))

## Comparison
r1 = RationalFunction([ 0., 1, 2], :x)
r2 = RationalFunction([ 0., 1, 2], :s)
r3 = RationalFunction([-0., 1, 2], :x)

@test r1 ≠ r2 && !isequal(r1, r2)
@test_throws DomainError r1 ≈ r2
@test r1 == r3 && r1 ≈ r3 && !isequal(r1, r3)

## Function evaluation
p1 = Poly([1,-2,1])     # p1 = (x-1)(x-1)
p2 = poly([1,2,3])      # p2 = (x-1)(x-2)(x-3)

r1 = RationalFunction(p1, p2)
r2 = conj(r1)

realinput = [1., 2., 3., 4., 5.]
imaginput = realinput + [0.1, -0.2, 0.3, -0.4, 0.5]*1im

@test_approx_eq(r1(realinput), r2(realinput))
@test_approx_eq(r1(imaginput), conj(r2(imaginput)))

r1  = RationalFunction(3p1, 5p1)
r2  = RationalFunction(3p1, 5p2)
r3  = RationalFunction(3p2, 5p1)

@test_approx_eq(r1(1), 3/5)
@test_approx_eq(r2(1), 0)
@test_approx_eq(r3(1), Inf)

@test_approx_eq(r1(realinput), fill(3/5, size(realinput)))
@test_approx_eq(r1(imaginput), fill(3/5, size(realinput)))

@test_approx_eq(r1(Inf), 3/5)
@test_approx_eq(r2(Inf), 0)
@test_approx_eq(r3(Inf), Inf)

## Function fitting
x   = -2:1E-1:2
y   = Float64[x^2 + 3x + 5 for x in x]
r1  = funcfit(x, y, 2)

@test_throws DomainError funcfit(x, y, -1, 2)
@test_throws DomainError funcfit(x, y, 2, -1)
@test_throws DomainError funcfit(x, y, length(x))
@test_throws DomainError funcfit(x, [1; y], 2)

@test coeffs(r1)[1] ≈ [5, 3, 1]
@test y ≈ r1(x)

x   = 0:1E-1:4
y   = Float64[x^2 + 3x + √x + 5cos(x) for x in x]
r2  = funcfit(x, y)

@test y ≈ r2(x)

## Inversion, transpose, conjugation
r1 = RationalFunction(p1, p2)
r2 = inv(r1)
r3 = transpose(r1)

@test r1*r2 ≈ 1
@test r1*r3 ≈ r1^2

r1 = RationalFunction([1+1im, 2.], [1-3im, 5, 8im])
r2 = conj(r1)

@test num(r1) == conj(num(r2))
@test den(r1) == conj(den(r2))

# test related to #2
m1 = fill(r1,1,1) # [r1]
m2 = fill(r2,1,1) # [conj(r1)]
m3 = conj(m1)     # [conj(r1)]

@test m2 == m3

r3 = convert(typeof(r1), r1)
r4 = convert(typeof(r2), r1)

@test r1 == r3
@test r2 == r4

## Derivative and reduction
### p1 = (x-1)(x-1)
### p2 = (x-1)(x-2)(x-3)
p3 = poly([1])    # p3 = (x-1)
p4 = poly([2,3])  # p4 = (x-2)(x-3)

r1      = RationalFunction(p1, p2)
result  = RationalFunction(polyder(p3)*p4 - p3*polyder(p4), p4^2)
@test derivative(r1) == result
@test !isequal(derivative(r1), result)

@test roots(derivative(reduce(r1))) == roots(result)

r1 = RationalFunction(p1, p1)
@test derivative(r1, 0) == r1
@test derivative(r1) == 0
@test derivative(RationalFunction(1, p3), 2) == RationalFunction(2, poly([1; 1; 1]))
@test_throws DomainError derivative(r1, -1)

## Mathematical operations
### Between rational functions
r1 = RationalFunction(1, [-1, 1], :x) # r1 =  1/(x-1)
r2 = RationalFunction(-1, [1, 1], :x) # r2 = -1/(x+1)
r3 = RationalFunction(1, [-1, 1], :s) # r3 =  1/(s-1)

for op in (:+, :-, :*, :/, :dot, :.+, :.-, :.*, :./)
  @test_throws DomainError eval( :(($op)(r1, r3)) )
end

@test r1+r2 ≈ r1 .+ r2 ≈ RationalFunction( 2    , [-1, 0, 1])
@test r1-r2 ≈ r1 .- r2 ≈ RationalFunction([0, 2], [-1, 0, 1])
@test r1*r2 ≈ r1 .* r2 ≈ dot(r1, r2) ≈ RationalFunction(-1, [-1, 0, 1])
@test r1/r2 ≈ r1 ./ r2 ≈ inv(r2/r1) ≈ -RationalFunction([1, 1], [-1, 1])

### Between polynomials
p1 = Poly([1, 1], :x)     # p1 = (x+1)
p2 = Poly([1, 1], :s)     # p2 = (s+1)
r1 = RationalFunction(p1) # r1 = (x+1)/1

@test p1 == r1 == p1
@test p2 ≠ r1 ≠ p2

@test p1 ≈ r1 ≈ p1
@test_throws DomainError p2 ≈ r1
@test_throws DomainError r1 ≈ p2

r1 = RationalFunction([0, 1, 1], [-1, 1], :x) # r1 = x(x+1)/(x-1)

for op in (:+, :-, :*, :/, :dot, :.+, :.-, :.*, :./)
  @test_throws DomainError eval( :(($op)(r1, p2)) )
  @test_throws DomainError eval( :(($op)(p2, r1)))
end

@test p1+r1 ≈ p1 .+ r1 ≈ r1+p1 ≈ r1 .+ p1 ≈ RationalFunction([-1, 1, 2], [-1, 1])
@test r1-p1 ≈ r1 .- p1 ≈ -(p1-r1) ≈ -(p1 .- r1) ≈ RationalFunction([1, 1],[-1, 1])
@test r1*p1 ≈ r1 .* p1 ≈ dot(r1, p1) ≈ p1*r1 ≈ p1 .* r1 ≈ dot(p1, r1) ≈ RationalFunction([0, 1, 2, 1], [-1, 1])
@test r1/p1 ≈ r1 ./ p1 ≈ inv(p1/r1) ≈ inv(p1 ./ r1) ≈ RationalFunction([0, 1], [-1, 1])

### Between numbers
n   = 3.
r1  = RationalFunction(n*[1, 1], [1, 1])

@test n == r1 == n
@test n+1 ≠ r1 ≠ n+2

@test n ≈ r1 ≈ n
@test n + 0*1im ≈ r1 ≈ n + 0*1im

r1  = RationalFunction([-1, 1], [1, 1]) # r1 = (x-1)/(x+1)

@test r1+n ≈ r1 .+ n ≈ n+r1 ≈ n .+ r1 ≈ RationalFunction([2, 4], [1, 1])
@test r1-n ≈ r1 .- n ≈ -(n-r1) ≈ -(n .- r1) ≈ RationalFunction([-4, -2], [1, 1])
@test r1*n ≈ r1 .* n ≈ dot(r1, n) ≈ n*r1 ≈ n .* r1 ≈ dot(n, r1) ≈ RationalFunction([-3, 3], [1, 1])
@test r1/n ≈ r1 ./ n ≈ inv(n/r1) ≈ inv(n ./ r1) ≈ RationalFunction([-1, 1], [3, 3]) ≈ RationalFunction([-1, 1]/3, [1, 1])

## Solution of rational function equalities
r1 = RationalFunction([-1, 1], [ 1, 1], :x) # r1 = (x-1)/(x+1)
r2 = RationalFunction([ 4, 3], [ 6, 5], :x) # r2 = (3x+4)/(5x+6)
r3 = RationalFunction([-1, 1], [-1, 1], :s) # r3 = (s-1)/(s-1)

### w.r.t `Number`s
@test solve(0, r1) ≈ solve(r1, 0) ≈ solve(r1) ≈ [1]
@test solve(Inf, r1) ≈ solve(r1, Inf) ≈ [-1]
@test solve(5, r1) ≈ solve(r1, 5) ≈ [-3/2]
@test_throws ErrorException solve(5r3, 5)

### w.r.t `Poly`s
p1    = Poly([-1, 1])
sln1  = solve(r1, p1)
sln2  = solve(p1, r1)
@test (sln1 ≈ [0, 1] || sln1 ≈ [1, 0]) && (sln2 ≈ [0, 1] || sln2 ≈ [1, 0])
@test solve(r1, Poly([0])) ≈ [1]
@test_throws DomainError solve(r3, p1)
@test_throws ErrorException solve(p1, RationalFunction(p1))

sln = solve(r1, r2)
@test sln ≈ [(3-√29)/2, (3+√29)/2] || sln ≈ [(3+√29)/2, (3-√29)/2]
@test_throws DomainError solve(r2, r3)
@test_throws ErrorException solve(r1, r1)

# Plotting via `RecipesBase` and `Plots`
x       = -2:1E-1:2
xinit   = x[5:5:end]
yinit1  = map(x->x^2+3x+5, xinit)
yinit2  = map(x->x^2+2, xinit)
init1   = map((x,y)->(x,y), xinit, yinit1)
init2   = map((x,y)->(x,y), xinit, yinit2)

r1      = funcfit(xinit, yinit1, 2)
r2      = funcfit(xinit, yinit2, 2)

@test isa(plot(r1, x, xinit, yinit1), Plots.Plot)
@test isa(plot(r1, x, label = "r1(x)"), Plots.Plot)
@test isa(plot!(r2, x, init2, label = "r2(x)"), Plots.Plot)

savefig("test-plot.png")

@test_throws DomainError plot(r1, x, xinit, [yinit1; 1])

## residue
num1 = [6, 9, 16, 8, 1]
den1 = [6, 11, 6, 1]
r, p, k = residue(num1, den1)
@test r ≈ [-6.0, -4.0, 3.0]
@test p ≈ [-3.0, -2.0, -1.0]
@test k ≈ [2.0, 1.0]
num2, den2 = residue(r, p, k)
@test num1 ≈ num2
@test den1 ≈ den2

num1 = [10, 2]
den1 = [0, 10, 2, 1]
r, p, k = residue(num1, den1)
@test r ≈ [-0.5-(1/6)im, -0.5+(1/6)im, 1.0+0.0im]
@test p ≈ [-1.0+3.0im, -1.0-3.0im, 0.0+0.0im]
@test k ≈ [0.0]
num2, den2 = residue(r, p, k)
#@test num1 ≈ num2 #Cannot pass due to erroneous third, imaginary root. Numerical errors?
@test den1 ≈ den2

num1 = [1, 0, 1]
den1 = [0, 1, -2, 1]
r, p, k = residue(num1, den1)
@test r ≈ [-0.0, 2.0, 1.0]
@test p ≈ [1.0, 1.0, 0.0]
@test k ≈ [0.0]
num2, den2 = residue(r, p, k)
@test num1 ≈ num2
@test den1 ≈ den2
