__precompile__(true)

module RationalFunctions

using Compat
using Polynomials
using RecipesBase

# Import conversion and promotion functions for overloading
import Base: convert, promote_rule

# Enable copy operations
import Base: copy

# Import iteration interface functions
import Base: eltype

# Import printing functions
import Base: summary
@compat import Base.show

# Import identities for overloading
import Base: one, zero

# Import comparison methods for overloading
import Base: isapprox, hash, isequal, transpose, conj

# Import mathematical operations for overloading
import Base: +, .+, -, .-, *, .*, dot, inv, /, ./, ==

# Import functions from Polynomials
import Polynomials: coeffs, degree, roots, variable, Poly

# Import num/den-type of functions
import Base: num, den, zeros, reduce

# Export only the useful functions
export  RationalFunction,
        # Polynomials-related functions
        coeffs,
        degree,
        roots,
        variable,
        residue,
        # Other functions
        derivative,
        funcfit,
        poles,
        solve

# Include files
include("type.jl")
include("printing.jl")
include("conversions.jl")
include("methods.jl")
include("plotting.jl")

end # module
