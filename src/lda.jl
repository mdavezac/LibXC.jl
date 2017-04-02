""" Adds check for functional type """
macro check_functional funcname functype
    msg = "Incorrect number of arguments: input is not an $functype functional"
    quote
        family($funcname) ≠ Constants.$functype && throw(ArgumentError($msg))
    end
end

""" Adds argument check for energy derivative availability """
macro check_availability funcname functype
    msg = "Functional does not implement energy."
    quote
        Constants.$functype ∉ flags($funcname) && error($msg)
    end
end
""" Adds argument check for size compatibility """
macro check_size funcname rhoname outname factor
    msg = "sizes of $rhoname and $outname are incompatible"
    quote
        if size($outname) ≠ output_size($funcname, $rhoname, $factor)
            throw(ArgumentError($msg))
        end
    end
end

"""
    $(SIGNATURES)

Computes the energy in-place for a given LDA functional. For spin-unpolarized functionals,
the output array has the dimensions of `ρ`. For spin-polarized functionals, assuming
`ndims(ρ) > 1 && size(ρ, 1) == 2`, it is `size(ρ)[2:end]`.
"""
function energy!(func::AbstractLibXCFunctional{Cdouble},
                 ρ::DenseArray{Units.ρ{Cdouble}},
                 ϵ::DenseArray{Units.ϵ{Cdouble}})
    energy!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ϵ))
    ϵ
end
function energy!(func::AbstractLibXCFunctional{Cdouble},
                 ρ::DenseArray{Cdouble},
                 ϵ::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func exc
    @check_size func ρ ϵ 1

    ccall((:xc_lda_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ϵ)
    ϵ
end


"""
    $(SIGNATURES)

Computes the potential in-place for a given LDA functional. For spin-unpolarized
functionals, the output array has the dimensions of `ρ`. For spin-polarized functionals,
assuming `ndims(ρ) > 1 && size(ρ, 1) == 2`, it is `size(ρ)`.
"""
function potential!(func::AbstractLibXCFunctional{Cdouble},
                    ρ::DenseArray{Units.ρ{Cdouble}},
                    ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}})
    potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂ϵ_∂ρ))
    ∂ϵ_∂ρ
end
function potential!(func::AbstractLibXCFunctional{Cdouble},
                    ρ::DenseArray{Cdouble},
                    ∂ϵ_∂ρ::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func vxc
    @check_size func ρ ∂ϵ_∂ρ 2

    ccall((:xc_lda_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ∂ϵ_∂ρ)
    ∂ϵ_∂ρ
end

"""
    $(SIGNATURES)

Computes the second energy derivative in-place for a given LDA functional. For
spin-unpolarized functionals, the output array has the dimensions of `ρ`. For spin-polarized
functionals, assuming `ndims(ρ) > 1 && size(ρ, 1) == 2`, it is `(3, size(ρ)[2:end]...)`.
"""
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Units.ρ{Cdouble}},
                                   ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}})
    second_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂²ϵ_∂²ρ))
    ∂²ϵ_∂²ρ
end
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble},
                                   ∂²ϵ_∂ρ²::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func fxc
    @check_size func ρ ∂²ϵ_∂ρ² 3

    ccall((:xc_lda_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ∂²ϵ_∂ρ²)

    ∂²ϵ_∂ρ²
end

"""
    $(SIGNATURES)

Computes the third energy derivative in-place for a given LDA functional. For
spin-unpolarized functionals, the output array has the dimensions of `ρ`. For spin-polarized
functionals, assuming `ndims(ρ) > 1 && size(ρ, 1) == 2`, it is `(4, size(ρ)[2:end]...)`.
"""
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Units.ρ{Cdouble}},
                                   ∂³ϵ_∂ρ³::DenseArray{Units.∂⁴ϵ_∂ρ⁴{Cdouble}})
    second_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂³ϵ_∂³ρ))
    ∂³ϵ_∂³ρ
end
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble},
                                   ∂³ϵ_∂ρ³::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func kxc
    @check_size func ρ ∂⁴ϵ_∂ρ⁴ 4

    ccall((:xc_lda_kxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ∂³ϵ_∂ρ³)

    ∂³ϵ_∂ρ³
end

""" Stateless function for computing LDA energies """
function energy(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}})
    energy!(func, ρ, similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
end
""" Stateless function for computing LDA first derivatives """
function potential(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}})
    potential!(func, ρ, similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)))
end
""" Stateless function for computing LDA second derivatives """
function second_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                                  ρ::DenseArray{Units.ρ{Cdouble}})
    second_energy_derivative!(func, ρ,
                              similar(ρ, Units.∂²ϵ_∂²ρ{Cdouble}, output_size(func, ρ, 3)))
end
""" Stateless function for computing LDA third derivatives """
function third_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                                 ρ::DenseArray{Units.ρ{Cdouble}})
    third_energy_derivative!(func, ρ,
                             similar(ρ, Units.∂³ϵ_∂³ρ{Cdouble}, output_size(func, ρ, 4)))
end

# Adds simplifying overloads
for name ∈ [:energy, :potential,:second_energy_derivative, :third_energy_derivative]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name{T <: Number}(func::AbstractLibXCFunctional{Cdouble},
                                    ρ::DenseArray{T})
            $name!(func, Units.conversion(Units.ρ, ρ))
        end

        @doc """
            $(SIGNATURES)

        A simple way to call a functional using it's name. Spin polarization is determined
        from the dimensionality of ρ: `ndims(ρ) > 1 && size(ρ, 1) == 2`.
        """ ->
        function $name!(name::Symbol, ρ::DenseArray, args...)
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end

        @doc """
            $(SIGNATURES)

        A simple way to call a functional using it's name. Spin polarization is explicitly
        requested.
        """ ->
        function $name!(name::Symbol, spin::Union{Constants.SPIN, Bool},
                        ρ::DenseArray, args...)
            $name!(XCFunctional(name, spin), ρ, args...)
        end

        @doc """
            $(SIGNATURES)

        A simple way to call a functional using it's name. Spin polarization is determined
        from the dimensionality of ρ: `ndims(ρ) > 1 && size(ρ, 1) == 2`.
        """ ->
        function $name(name::Symbol, ρ::DenseArray, args...)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end

        @doc """
        $(SIGNATURES)

        A simple way to call a functional using it's name. Spin polarization is explicitly
        specified.
        """ ->
        function $name(name::Symbol, spin::Union{Bool, Constants.SPIN},
                       ρ::DenseArray, args...)
            $name(XCFunctional(name, spin), ρ, args...)
        end
    end
end

""" All outputs from LDA """
typealias AllLDA @NT(energy, potential, second_derivative, third_derivative)
"""
    $(SIGNATURES)

LDA energy, first, second, and third derivatives.
The first, second, and third derivatives are optional. It is an error to request a
derivative of the functional that is not implemented in the underlying C library.
"""
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
        AllLDA(reinterpret(Units.ϵ{Cdouble}, εxc),
               map(i -> reinterpret(EnergyDensity{Cdouble}, outputs[i]), outputs)..., [], [])
    else
        AllLDA(reinterpret(EnergyDensity{Cdouble}, εxc),
               map(x -> reinterpret(EnergyDensity{Cdouble}, x), outputs)...)
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
    LDAEnergyPotential(reinterpret(EnergyDensity{Cdouble}, εxc),
                       reinterpret(EnergyDensity{Cdouble}, potential))
end

function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               εxc::EnergyDensityArray{Cdouble},
                               potential::EnergyDensityArray{Cdouble})
    energy_and_potential!(func, ρ, reinterpret(Cdouble, ϵxc),
                          reinterpret(Cdouble, potential))
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

""" Computes the energy and all available derivatives for the given functional """
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
