# Relation with itself
promote_rule{T,S,U1,V1,U2,V2}(::Type{RationalFunction{Var{T},Conj{S},U1,V1}},
  ::Type{RationalFunction{Var{T},Conj{S},U2,V2}}) =
  RationalFunction{Var{T},Conj{S},promote_type(U1,U2),promote_type(V1,V2)}
convert{T,S,U1,V1,U2,V2}(::Type{RationalFunction{Var{T},Conj{S},U1,V1}},
  r::RationalFunction{Var{T},Conj{S},U2,V2}) =
  RationalFunction(convert(Poly{U1}, r.num), convert(Poly{V1}, r.den), Conj{S})
convert{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}},
  r::RationalFunction{Var{T},Conj{S},U,V}) = r

# Relation with numbers
promote_rule{T,S,U,V,Y<:Number}(::Type{RationalFunction{Var{T},Conj{S},U,V}}, ::Type{Y}) =
  RationalFunction{Var{T},Conj{S},promote_type(U,Y),V}
convert{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}}, n::Number) =
  RationalFunction(Poly(U[n], T), Poly([one(V)], T), Conj{S})

# Relation with polynomials
promote_rule{T,S,U,V,Y<:Number}(::Type{RationalFunction{Var{T},Conj{S},U,V}},
  ::Type{Poly{Y}}) = RationalFunction{Var{T},Conj{S},promote_type(U,Y),V}
function convert{T,S,U,V}(::Type{RationalFunction{Var{T},Conj{S},U,V}}, p::Poly)
  @assert T == p.var "RationalFunction{Var{$T}}(p::Poly): p.var = $(p.var) (≠ $T)"
  newpoly = convert(Poly{U}, p)
  RationalFunction(newpoly, Poly([one(V)], T), Conj{S})
end
