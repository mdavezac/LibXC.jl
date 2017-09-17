module FunctionalMacros
import LibXC: energy, energy!, potential, potential!, energy_and_potential,
              energy_and_potential!, second_energy_derivative, second_energy_derivative!,
              third_energy_derivative, third_energy_derivative!, lda, lda!
using ..Constants
using ..Internals
using ..Internals: AbstractLibXCFunctional, CFuncType, libxc
using ..ArrayManips: from_libxc_array!, to_libxc_array, valid_array
using ..OutputTuples 

using AxisArrays
using ArgCheck

using DFTShims: ColinearSpinFirst, SpinDegenerate, Dispatch, is_spin_polarized, components, 
                SpinCategory, SpinAware
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end

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
    axis_args = (arg -> :($arg::DD.AxisArrays.$arg)).(allargs)
    array_args = (arg -> :($arg::DD.DenseArrays.$arg{Float64})).(allargs)

    # converts backt to original axisarrays
    name! = Symbol(name, :!)
    # input sanity checks
    arg_checks = x -> :(valid_array(ρ, $x))
    arg_checks_array = x -> :(valid_array(ρ, $x, Spin))
    # converts to type understood by C LibXC
    to_libxc_types = x -> :($(Symbol(:c, x)) = to_libxc_array(ρ, $x))
    to_libxc_array_types = x -> :($(Symbol(:c, x)) = to_libxc_array($x))
    to_dense_types = x -> :(reinterpret(Float64, $(Symbol(:c, x)).data))
    to_c_types = x -> :(reinterpret(Float64, $x))
    from_libxc_types = x -> :(from_libxc_array!($x, $(Symbol(:c, x))))

    quote
        $(esc(name!))(func::AbstractLibXCFunctional, $(axis_args...)) = begin
            $(arg_checks.(allargs)...)
            if family(func) ≠ getfield(Constants, $(QuoteNode(dfttype)))
                throw(ArgumentError("Incorrect functional type $($(QuoteNode(dfttype)))"))
            end
            reqout = _requested_outputs($(outputs...))
            if !(reqout ⊆ flags(func))
                throw(ArgumentError("This functional does not implement all of $reqout"))
            end

            $(to_libxc_types.(allargs)...)
            $(esc(name!))(func, $(to_dense_types.(allargs)...))
            $(from_libxc_types.(outputs)...)
            $(esc(output_type))($(outputs...))
        end

        $(esc(name!))(func::AbstractLibXCFunctional, $(array_args...)) = begin
            const Spin = SpinCategory(func)
            $(arg_checks_array.(allargs)...)
            $(to_libxc_array_types.(allargs)...)
            $(esc(name!))(func, $(to_c_types.(array_args)...))
            $(esc(output_type))($(from_libxc_types.(outputs)...))
        end
    end
end

""" Constructs non-mutating functions taking AxisArray arguments """
_nonmutating_wrapper_functionals(name::Symbol, dfttype::Symbol,
                                 outputs::Tuple{Symbol, Vararg{Symbol}}) = begin
    @lintpragma("Ignore unused outputs")

    name! = Symbol(name, :!)
    inputs = dfttype == :lda ? (:ρ, ) : (:ρ, :σ)
    args_axis = x -> :($x::DD.AxisArrays.$x)
    args_array = x -> :($x::DD.Arrays.$x)
    similars = x -> begin
        if x == :ϵ
            :(similar(ρ, DH.Scalars.$x{Float64}, SpinDegenerate()))
        else
            :(similar(ρ, DH.Scalars.$x{Float64}, SpinAware()))
        end
    end
    similars_nospin = x -> begin
        @lintpragma("Ignore unused x")
        :(similar(ρ, DH.Scalars.$x{Float64}))
    end
    similars_spin = x -> begin
        if x == :ϵ
            :(similar(ρ, DH.Scalars.$x{Float64}, Base.tail(size(ρ))))
        else
            :(similar(ρ, DH.Scalars.$x{Float64},
                      length(components(DH.Scalars.$x, ColinearSpinFirst())),
                      Base.tail(size(ρ))...))
        end
    end

    quote
        $(esc(name))(func::AbstractLibXCFunctional, $(args_axis.(inputs)...)) =
            $(esc(name!))(func, $(inputs...), $(similars.(outputs)...))
        $(esc(name))(name::Symbol, $(args_axis.(inputs)...)) =
            $(esc(name))(XCFunctional(name, is_spin_polarized(ρ)), $(inputs...))
        $(esc(name))(func::AbstractLibXCFunctional, $(args_array.(inputs)...)) = begin
            if spin(func) == Constants.polarized
                $(esc(name!))(func, $(inputs...), $(similars_nospin.(outputs)...))
            else
                $(esc(name!))(func, $(inputs...), $(similars_nospin.(outputs)...))
            end
        end
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
