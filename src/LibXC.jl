__precompile__()
module LibXC
export description, kind, family, flags, citations, spin, energy, energy!
export potential, potential!, second_energy_derivative, third_energy_derivative
export energy_and_potential, energy_and_potential!, lda!, lda, XCFunctional, gga, gga!
export libxc_functionals

using DocStringExtensions
using Unitful: Quantity, @u_str
using DFTShims
export @u_str

include("constants.jl")
using .Constants

include("internals.jl")
using .Internals

# include("Checks.jl")
# using .Checks: *

# include("named_tuples.jl")
# include("lda.jl")
# include("gga.jl")
# include("overloads.jl")
end # module
