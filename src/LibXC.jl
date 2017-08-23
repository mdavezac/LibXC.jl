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

if Cdouble !== Float64
    error("No idea what happens if Cdouble ≠ Float64")
end

include("constants.jl")
using .Constants

include("internals.jl")
using .Internals

include("Checks.jl")
using .Checks

include("OutputTuples.jl")
using .OutputTuples

function lda end
function lda! end
function gga end
function gga! end
function energy end
function energy! end
function potential end
function potential! end
function energy_and_potential! end
function energy_and_potential end
function second_energy_derivative end
function second_energy_derivative! end
function third_energy_derivative end
function third_energy_derivative! end

include("FunctionalMacros.jl")
using .FunctionalMacros
include("lda.jl")
using .LDA
include("gga.jl")
using .GGA
# include("overloads.jl")
end # module
