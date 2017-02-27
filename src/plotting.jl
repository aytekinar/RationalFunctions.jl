@recipe function f{T,S,U<:Real,V<:Real,W<:Real,Z<:Real}(r::RationalFunction{Val{T},Val{S},U,V},
  x::AbstractVector, xinit::AbstractVector{W}, yinit::AbstractVector{Z})
  if length(xinit) ≠ length(yinit)
    warn("plot(r, x, xinit, yinit): length(xinit) ≠ length(yinit)")
    throw(DomainError())
  end

  # Some defaults
  xguide  --> "Input ($T)"
  yguide  --> "Output (y)"

  @series begin
    primary :=  true
    label   --> "r($T)"
    x, r(x)
  end

  if !isempty(xinit)
    @series begin
      primary     :=  false
      seriestype  :=  :scatter
      label       --> "Initial Points"
      xinit, yinit
    end
  end
end

@recipe f{T,S,U<:Real,V<:Real}(r::RationalFunction{Val{T},Val{S},U,V},
  x::AbstractVector) = (r, x, Int[], Int[])

@recipe function f{T,S,U<:Real,V<:Real,Z<:Tuple}(r::RationalFunction{Val{T},Val{S},U,V},
  x::AbstractVector, init::AbstractArray{Z})
  xinit = map(x->x[1], init)
  yinit = map(x->x[2], init)
  (r, x, xinit, yinit)
end
