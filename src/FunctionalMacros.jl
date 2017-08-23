module FunctionalMacros
import LibXC: energy, energy!, potential, potential!, energy_and_potential,
              energy_and_potential!, second_energy_derivative, second_energy_derivative!,
              third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Constants
using ..Internals
using ..Internals: AbstractLibXCFunctional, CFuncType, libxc
using ..OutputTuples 

using AxisArrays

using DFTShims: ColinearSpinFirst, Dispatch, is_spin_polarized, components
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

""" Conversion between AxisArray of different unit/types to something LibXC understands """
_conversion(etype::Type, array::AxisArray) = begin
    eltype(array) === etype && return array
    result = similar(array, etype)
    result .= array
    result
end

_requested_outputs(o::Vararg{DD.AxisArrays.All}) = begin
    @lintpragma("Ignore unused o")
    check = (x, y) -> eltype(x) <: y
    result = Constants.FLAGS[]
    eltypes = eltype.(o)
    findfirst(x -> x <: DD.Scalars.ϵ, eltypes) ≠ 0 && push!(result, Constants.exc)
    findfirst(x -> x <: DD.Scalars.∂ϵ_∂ρ, eltypes) ≠ 0 && push!(result, Constants.vxc)
    findfirst(x -> x <: DD.Scalars.∂²ϵ_∂ρ², eltypes) ≠ 0 && push!(result, Constants.fxc)
    findfirst(x -> x <: DD.Scalars.∂³ϵ_∂ρ³, eltypes) ≠ 0 && push!(result, Constants.kxc)
    unique(result)
end

""" Constructs mutating functions taking AxisArray arguments """
_mutating_wrapper_functionals(name::Symbol, dfttype::Symbol, output_type::Symbol,
                              outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
   
    allargs = dfttype == :lda ? (:ρ, outputs...): (:ρ, :σ, outputs...)

    # function arguments
    ddargs = (arg -> :($arg::DD.AxisArrays.$arg)).(allargs)
    dhargs = (arg -> :($arg::DH.AxisArrays.$arg)).(allargs)

    # conversion to axisarrays compatible with LibXC
    conversions = (arg -> :(_conversion(DH.Scalars.$arg{Float64}, $arg))).(allargs)
    # converts backt to original axisarrays
    back = if length(outputs) == 1
        arg = outputs[1]
        [:(AxisArray(convert(typeof($(arg).data), result.data), axes($(arg))))]
    else
        (x -> :(AxisArray(convert(typeof($x.data), result.$x.data), axes($x)))).(outputs)
    end
    # mutating function names
    name! = Symbol(name, :!)
    # input sanity checks
    arg_checks = x -> begin 
        msg = "Axes of ρ and $x do not match"
        spinmsg = "First axis of $x is not spin"
        spinerror = "Axis values of spin axis of $x are incorrect"
        comps = gensym("comps")
        quote
            if is_spin_polarized(ρ)
                let comps = components(eltype($x), ColinearSpinFirst())
                    if length(comps) > 1
                        axisnames($x)[1] ≠ :spin && throw(ArgumentError($spinmsg))
                        axisvalues($x)[1] ≠ comps && throw(ArgumentError($spinerror))
                        axes(ρ)[2:end] == axes($x)[2:end] || throw(ArgumentError($msg))
                    end
                end
            elseif axes(ρ) ≠ axes($x)
                throw(ArgumentError($msg))
            end
        end
    end
    # converts to type understood by C LibXC
    ctypes = (x -> :(reinterpret(typeof(one(eltype($x))), convert(Array, $x)))).(allargs)

    quote
        $(esc(name!))(func::AbstractLibXCFunctional, $(ddargs...)) = begin
            result = $(esc(name!))(func, $(conversions...))
            $(esc(output_type))($(back...))
        end
        $(esc(name!))(func::AbstractLibXCFunctional, $(dhargs...)) = begin
            $(arg_checks.(allargs)...)
            if family(func) ≠ getfield(Constants, $(QuoteNode(dfttype)))
                throw(ArgumentError("Incorrect functional type $($(QuoteNode(dfttype)))"))
            end
            reqout = _requested_outputs($(outputs...))
            if !(reqout ⊆ flags(func))
                throw(ArgumentError("This functional does not implement all of $reqout"))
            end

            $(esc(name!))(func, $(ctypes...))

            $(esc(output_type))($(outputs...))
        end
    end
end

""" Constructs non-mutating functions taking AxisArray arguments """
_nonmutating_wrapper_functionals(name::Symbol, dfttype::Symbol,
                                 outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
    @lintpragma("Ignore unused outputs")
    @lintpragma("Ignore unused x")
   
    name! = Symbol(name, :!)
    inputs = dfttype == :lda ? (:ρ, ) : (:ρ, :σ)
    args = (x -> :($x::DD.AxisArrays.$x)).(inputs)
    similars = (x -> :(similar(DH.Scalars.$x{Float64}, ρ))).(outputs)

    quote
        $(esc(name))(func::AbstractLibXCFunctional, $(args...)) =
            $(esc(name!))(func, $(inputs...), $(similars...))

        $(esc(name))(name::Symbol, $(args...)) =
            $(esc(name))(XCFunctional(name, is_spin_polarized(ρ)), $(inputs...))
    end
end

""" Constructs functions which actually call the C implementations """
_c_wrapper_functionals(name::Symbol, dfttype::Symbol, cname::Symbol,
                       outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
    allargs = dfttype == :lda ? (:ρ, outputs...): (:ρ, :σ, outputs...)
    name! = Symbol(name, :!)
    nonnull = collect(Iterators.filter(x -> x ≠ :C_NULL, allargs))
    dargs = (arg -> :($arg::DenseArray{Float64})).(nonnull)
    quote
        $(esc(name!))(func::AbstractLibXCFunctional, $(dargs...)) = begin
            @lintpragma("Ignore use of undeclared variable ρ")
            @lintpragma("Ignore unused u")
            @lintpragma("Ignore use of undeclared variable ccall")
            ccall(($(QuoteNode(cname)), $libxc), Void,
                  (Ptr{CFuncType}, Cint, $([:(Ptr{Float64}) for u in allargs]...)),
                  func.c_ptr, length(ρ) / (spin(func) == Constants.polarized ? 2: 1),
                  $(allargs...))
            tuple($(nonnull...))
        end
    end
end

macro _wrapper_functionals(func, ftype, cname, output_type, outputs...)
    quote
        $(_mutating_wrapper_functionals(func, ftype, output_type, outputs))
        $(_nonmutating_wrapper_functionals(func, ftype, outputs))
        $(_c_wrapper_functionals(func, ftype, cname, outputs))
    end
end

macro _all_wrapper_functionals(n::Int, func, cname, output_type, outputs...)
    @lintpragma("Ignore unused i")
    o = (outputs..., (:C_NULL for i in (length(outputs) + 1):n)...)
    quote
        $(_mutating_wrapper_functionals(func, func, output_type, outputs))
        $(_c_wrapper_functionals(func, func, cname, o))
    end
end
end
