for (funcname, name, factor) ∈ [(:xc_lda_exc, :energy, 1),
                                (:xc_lda_vxc, :potential, 2),
                                (:xc_lda_fxc, :second_energy_derivative, 3),
                                (:xc_lda_kxc, :third_energy_derivative, 4)]
    local name! = Symbol("$(name)!")
    local has_it = Dict(:energy => Constants.exc, :potential => Constants.vxc,
                        :second_energy_derivative => Constants.fxc,
                        :third_energy_derivative => Constants.kxc)[name]
    @eval begin
        function $name!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                        output::DenseArray{Cdouble})
            if family(func) ≠ Constants.lda
                msg = "Incorrect number of arguments: input is not an LDA functional"
                throw(ArgumentError(msg))
            end
            if $has_it ∉ flags(func)
                error("Functional does not implement the $(replace(name, "_", " ")).")
            end
            if size(output) ≠ output_size(func, ρ, $factor)
                throw(ArgumentError("sizes of ρ and input are incompatible"))
            end

            ccall(($(parse(":$funcname")), libxc), Void,
                  (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
                  func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, output)
            output
        end

        function $name(func::AbstractLibXCFunctional, ρ::DenseArray)
            $name!(func, ρ, similar(ρ, eltype(ρ), output_size(func, ρ, $factor)))
        end
    end
end

# Adds simplifying overloads
for name ∈ [:energy, :potential,:second_energy_derivative, :third_energy_derivative]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(name::Symbol, ρ::DenseArray{Cdouble}, args...)
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name!(name::Symbol, s::Union{Constants.SPIN, Bool},
                        ρ::DenseArray{Cdouble}, args...)
            $name!(XCFunctional(name, s), ρ, args...)
        end
        function $name(name::Symbol, ρ::DenseArray, args...)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name(name::Symbol, s::Union{Bool, Constants.SPIN}, ρ::DenseArray, args...)
            $name(XCFunctional(name, s), ρ, args...)
        end
    end
end

""" All outputs from LDA """
typealias AllLDA @NT(energy, potential, second_derivative, third_derivative)
""" LDA energy, first, second, and third derivatives """
function lda!{T <: DenseArray{Cdouble}}(func::AbstractLibXCFunctional{Cdouble},
                                        ρ::DenseArray{Cdouble}, εxc::DenseArray{Cdouble},
                                        outputs::Vararg{T})
    if family(func) ≠ Constants.lda
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end

    if length(outputs) == 0
        return AllLDA(energy!(func, ρ, εxc), [], [], [])
    elseif length(outputs) == 1
        result = energy_and_potential!(func, ρ, εxc, outputs...)
        return AllLDA(result[1], [], [])
    end
    if length(outputs) ∈ (2, 3) && Constants.fxc ∉ flags(func)
        throw(ArgumentError("Functional does not implement second energy derivative"))
    elseif length(outputs) == 3 && Constants.kxc ∉ flags(func)
        throw(ArgumentError("Functional does not implement third energy derivative"))
    elseif length(outputs) ∉ (2, 3)
        throw(ArgumentError("Incorrect number of outputs, expected between 1 and 4"))
    end
    if size(εxc) ≠ output_size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(outputs[1]) ≠ output_size(func, ρ, 2)
        throw(ArgumentError("sizes of ρ and potential are incompatible"))
    end
    if size(outputs[2]) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and second derivative are incompatible"))
    end
    if size(outputs[3]) ≠ output_size(func, ρ, 4)
        throw(ArgumentError("sizes of ρ and third derivative are incompatible"))
    end

    args = tuple(outputs..., (C_NULL for i in 1:(3 - length(outputs)))...)

    ccall((:xc_lda, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, εxc, args[1], args[2], args[3])

    if length(outputs) == 2
        AllLDA(εxc, outputs..., [], [])
    else
        AllLDA(εxc, outputs...)
    end
end

""" Energy and potential from LDA """
typealias LDAEnergyPotential @NT(energy, potential)
""" LDA energy and first derivative """
function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               εxc::DenseArray{Cdouble}, potential::DenseArray{Cdouble})
    if family(func) ≠ Constants.lda
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end
    if size(εxc) ≠ output_size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(potential) ≠ output_size(func, ρ, 2)
        throw(ArgumentError("sizes of ρ and potential are incompatible"))
    end

    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, εxc, potential)
    LDAEnergyPotential(εxc, potential)
end

for name ∈ [:energy_and_potential, :lda, :gga]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(name::Symbol, ρ::DenseArray{Cdouble}, args...)
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name!(name::Symbol, s::Union{Bool, Constants.SPIN},
                        ρ::DenseArray{Cdouble}, args...)
            $name!(XCFunctional(name, s), ρ, args...)
        end
        function $name(name::Symbol, ρ::DenseArray{Cdouble}, args...)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name(name::Symbol, s::Union{Bool, Constants.SPIN},
                       ρ::DenseArray{Cdouble}, args...)
            $name(XCFunctional(name, s), ρ, args...)
        end
    end
end
function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Cdouble})
    energy_and_potential!(func, ρ,
                          similar(ρ, eltype(ρ), output_size(func, ρ, 1)), similar(ρ))
end
function lda(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble})
    if family(func) ≠ Constants.lda
        throw(ArgumentError("Functional is not a GGA functional"))
    end

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        lda!(func, ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 2)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 4)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        lda!(func, ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 2)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        lda!(func, ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 2)))
    elseif Constants.exc ∈ f
        lda!(func, ρ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
