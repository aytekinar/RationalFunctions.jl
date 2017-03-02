# Relation with itself
promote_rule{T,S,U1,V1,U2,V2}(::Type{RationalFunction{Val{T},Val{S},U1,V1}},
  ::Type{RationalFunction{Val{T},Val{S},U2,V2}}) =
  RationalFunction{Val{T},Val{S},promote_type(U1,U2),promote_type(V1,V2)}
convert{T,S,U1,V1,U2,V2}(::Type{RationalFunction{Val{T},Val{S},U1,V1}},
  r::RationalFunction{Val{T},Val{S},U2,V2}) =
  RationalFunction(convert(Poly{U1}, r.num), convert(Poly{V1}, r.den), Val{S})
convert{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}},
  r::RationalFunction{Val{T},Val{S},U,V}) = r
# Conversion between :conj and :nonc RF's (related to #2)
convert{T,S1,S2,U1,V1,U2,V2}(::Type{RationalFunction{Val{T},Val{S1},U1,V1}},
  r::RationalFunction{Val{T},Val{S2},U2,V2}) =
  RationalFunction(convert(Poly{U1}, conj(r.num)), convert(Poly{V1}, conj(r.den)),
  Val{S1})

# Relation with numbers
promote_rule{T,S,U,V,Y<:Number}(::Type{RationalFunction{Val{T},Val{S},U,V}}, ::Type{Y}) =
  RationalFunction{Val{T},Val{S},promote_type(U,Y),promote_type(V,Y)}
convert{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}}, n::Number) =
  RationalFunction(Poly(U[n], T), Poly([one(V)], T), Val{S})

# Relation with polynomials
promote_rule{T,S,U,V,Y<:Number}(::Type{RationalFunction{Val{T},Val{S},U,V}},
  ::Type{Poly{Y}}) = RationalFunction{Val{T},Val{S},promote_type(U,Y),promote_type(V,Y)}
function convert{T,S,U,V}(::Type{RationalFunction{Val{T},Val{S},U,V}}, p::Poly)
  newpoly = convert(Poly{U}, p)
  RationalFunction(newpoly, Poly([one(V)], T), Val{S})
end
