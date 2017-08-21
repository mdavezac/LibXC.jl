module Checks
using ..Internals: family, flags, spin, AbstractLibXCFunctional
using ..Constants
using DFTShims: SpinCategory, SpinDegenerate, ColinearSpinFirst, Dispatch, components,
                add_spin_axis, ColinearSpin
import DFTShims: SpinCategory
using ...LibXC
const DD = Dispatch.Dimensions

export @check_functional, @check_availability, @check_size

const SizeTuple = Tuple{Integer, Vararg{Integer}}

output_size(T::Type{<: DD.Scalars.All}, S::SpinCategory, d::SizeTuple) =
    add_spin_axis(S, d, length(components(T, S)))

SpinCategory(func::AbstractLibXCFunctional) =
    spin(func) == Constants.polarized ? ColinearSpinFirst(): SpinDegenerate()

output_size(T::Type{<: DD.Scalars.All}, sp::Bool, d::SizeTuple) =
    output_size(T, sp ? ColinearSpinFirst(): SpinDegenerate(), d)
output_size(T::Type{<: DD.Scalars.All}, s::Constants.SPIN, d::SizeTuple) =
    output_size(T, s == Constants.polarized, d)
output_size(T::Type{<: DD.Scalars.All}, func::AbstractLibXCFunctional, d::SizeTuple) =
    output_size(T, spin(func), d)

spinless_size(::SpinDegenerate, d::SizeTuple) = d
spinless_size(::ColinearSpin, d::SizeTuple) = Base.last(Base.IteratorsMD.split(d, Val{:1}))

""" Adds check for functional type """
macro check_functional(funcname, functype)
    msg = "Incorrect number of arguments: input is not an $functype functional"
    quote
        if LibXC.family($(esc(funcname))) ≠ LibXC.Constants.$functype
            throw(ArgumentError($msg))
        end
    end
end

""" Adds argument check for energy derivative availability """
macro check_availability(funcname, functype)
    msg = "Functional does not implement energy."
    quote
        LibXC.Constants.$functype ∉ LibXC.flags($(esc(funcname))) && error($msg)
    end
end

valid_size(S::SpinCategory, checker::AbstractArray, checkee::AbstractArray) =
    Base.size(checkee) == output_size(eltype(checkee), S, checker.size(checker))

""" Adds argument check for size compatibility """
macro check_size(spin_type, rhoname, outname)
    msg = "sizes of $rhoname and $outname are incompatible"
    quote
        if LibXC.Checks.valid_size($spin_type, $(esc(rhoname)), $(esc(outname)))
            throw(ArgumentError($msg))
        end
    end
end
end
