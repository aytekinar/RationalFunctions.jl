summary{T,S,U,V}(::RationalFunction{Val{T},Val{S},U,V})  = "RF{Val{$T},Val{$S},$U,$V}"

# Compact representations
function _compact{T,S}(stream, ::MIME"text/plain", r::RationalFunction{Val{T},Val{S}})
  var = ifelse(S == :conj, "$(T)̄", "$(T)")
  # print(stream, "num($(var))/den($(var))")
  print(stream, "n($(var))/d($(var))")
end

function _compact{T,S}(stream, ::MIME"text/latex", r::RationalFunction{Val{T},Val{S}})
  var = ifelse(S == :conj, "\\bar{$(T)}", "$(T)")
  print(stream, "\$")
  # print(stream, "\\tfrac{\\mathrm{num}($(var))}{\\mathrm{den}($(var))}")
  print(stream, "\\tfrac{n($(var))}{d($(var))}")
  print(stream, "\$")
end

# TODO: Think about text/html

# Full representations
function _full{T,S}(stream, m::MIME"text/plain", r::RationalFunction{Val{T},Val{S}})
  var = ifelse(S == :conj, "$(T)̄", "$(T)")
  println(stream, "f($(var)) = num($(var))/den($(var)), where,")
  print(stream, "num($(T)) is ")
  show(stream, m, r.num)
  println(stream, ", and,")
  print(stream, "den($(T)) is ")
  show(stream, m, r.den)
  print(stream, ".")
end

function _full{T,S}(stream, m::MIME"text/latex", r::RationalFunction{Val{T},Val{S}})
  var = ifelse(S == :conj, "\\bar{$(T)}", "$(T)")
  print(stream, "\$\$")
  print(stream, "f($(var)) = \\frac{\\mathrm{num}($(var))}{\\mathrm{den}($(var))}\\,,")
  print(stream, "\$\$ where \$\\mathrm{num}($(T))\$ is ")
  show(stream, m, r.num)
  print(stream, ", and \$\\mathrm{den}($(T))\$ is ")
  show(stream, m, r.den)
  print(stream, ".")
end

# TODO: Think about text/html

# `show` function
@compat Base.show(stream::IO, r::RationalFunction)                          =
  Base.show(stream, MIME("text/plain"), r)
@compat Base.show(stream::IO, mime::MIME"text/plain", r::RationalFunction)  =
  get(stream, :compact, false) ? _compact(stream, mime, r) : _full(stream, mime, r)
@compat Base.show(stream::IO, mime::MIME"text/latex", r::RationalFunction)  =
  get(stream, :compact, false) ? _compact(stream, mime, r) : _full(stream, mime, r)
