__precompile__()
module LibXC
export description, kind, family, flags, citations, spin, energy, energy!,
       first_energy_derivative, first_energy_derivative!, second_energy_derivative,
       third_energy_derivative,  energy_and_first_derivative, energy_and_first_derivative!,
       lda!, lda, XCFunctional, gga, gga!, libxc_functionals

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
DFTShims.SpinCategory(f::XCFunctional) =
    spin(f) == Constants.polarized ? ColinearSpinFirst() : SpinDegenerate()

include("OutputTuples.jl")
using .OutputTuples

function lda end
function lda! end
function gga end
function gga! end
function energy end
function energy! end
function first_energy_derivative end
function first_energy_derivative! end
function energy_and_first_derivative! end
function energy_and_first_derivative end
function second_energy_derivative end
function second_energy_derivative! end
function third_energy_derivative end
function third_energy_derivative! end

include("ArrayManips.jl")
using .ArrayManips

include("FunctionalMacros.jl")
using .FunctionalMacros


include("lda.jl")
using .LDA
include("gga.jl")
using .GGA

""" setups docs, useful for docs and doctests """
_doc_pages() = Any["Home" => "index.md"]
end # module
