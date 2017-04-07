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
    second_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂²ϵ_∂ρ²))
    ∂²ϵ_∂ρ²
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
function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                  ρ::DenseArray{Units.ρ{Cdouble}},
                                  ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}})
    third_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂³ϵ_∂ρ³))
    ∂³ϵ_∂ρ³
end
function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                  ρ::DenseArray{Cdouble},
                                  ∂³ϵ_∂ρ³::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func kxc
    @check_size func ρ ∂³ϵ_∂ρ³ 4

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
                              similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)))
end
""" Stateless function for computing LDA third derivatives """
function third_energy_derivative(func::AbstractLibXCFunctional{Cdouble},
                                 ρ::DenseArray{Units.ρ{Cdouble}})
    third_energy_derivative!(func, ρ,
                             similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)))
end

""" All outputs from LDA """
typealias AllLDA @NT(ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³)
"""
    $(SIGNATURES)

LDA energy, first, second, and third derivatives.
The first, second, and third derivatives are optional. It is an error to request a
derivative of the functional that is not implemented in the underlying C library.
"""
function lda!{T <: DenseArray{Cdouble}}(func::AbstractLibXCFunctional{Cdouble},
                                        ρ::DenseArray{Cdouble}, ϵ::DenseArray{Cdouble},
                                        outputs::Vararg{T})
    if length(outputs) == 0
        return AllLDA(energy!(func, ρ, ϵ), [], [], [])
    elseif length(outputs) == 1
        return AllLDA(energy_and_potential!(func, ρ, ϵ, outputs[1])..., [], [])
    elseif length(outputs) > 4
        throw(ArgumentError("Too many arguments"))
    end

    @check_functional func lda
    @check_availability func exc
    @check_availability func vxc
    length(outputs) >= 2 && @check_availability func fxc
    length(outputs) == 3 && @check_availability func kxc
    @check_size func ρ ϵ 1
    @check_size func ρ outputs[1] 2
    length(outputs) >= 2 && @check_size func ρ outputs[2] 3
    length(outputs) == 3 && @check_size func ρ outputs[3] 4
    length(outputs) > 3 && throw(ArgumentError("Too many arguments"))

    args = tuple(outputs..., (C_NULL for i in 1:(3 - length(outputs)))...)

    ccall((:xc_lda, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ϵ, args[1], args[2], args[3])

    if length(outputs) == 2
        AllLDA(ϵ, outputs..., [], [])
    else
        AllLDA(ϵ, outputs...)
    end
end
function lda!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
              ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
              ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}})
    result = lda!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ϵ),
                  reinterpret(Cdouble, ∂ϵ_∂ρ), reinterpret(Cdouble, ∂²ϵ_∂ρ²),
                  reinterpret(Cdouble, ∂³ϵ_∂ρ³))
    AllLDA(ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², ∂³ϵ_∂ρ³)
end
function lda!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
              ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}})
    result = lda!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ϵ),
                  reinterpret(Cdouble, ∂ϵ_∂ρ), reinterpret(Cdouble, ∂²ϵ_∂ρ²))
    AllLDA(ϵ, ∂ϵ_∂ρ, ∂²ϵ_∂ρ², Units.∂³ϵ_∂ρ³{Cdouble}[])
end
function lda!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}})
    AllLDA(energy_and_potential!(ρ, ϵ, ∂ϵ_∂ρ)...,
           Units.∂²ϵ_∂ρ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[])
end
function lda!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}})
    AllLDA(energy!(ρ, ϵ), Units.∂ϵ_∂ρ{Cdouble}[],
           Units.∂²ϵ_∂ρ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[])
end

""" Energy and potential from LDA """
typealias LDAEnergyPotential @NT(ϵ, ∂ϵ_∂ρ)
""" LDA energy and first derivative """
function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               ϵ::DenseArray{Cdouble}, ∂ϵ_∂ρ::DenseArray{Cdouble})
    @check_functional func lda
    @check_availability func exc
    @check_availability func vxc
    @check_size func ρ ϵ 1
    @check_size func ρ ∂ϵ_∂ρ 2

    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, ϵ, ∂ϵ_∂ρ)

    LDAEnergyPotential(ϵ, ∂ϵ_∂ρ)
end

function energy_and_potential!(func::AbstractLibXCFunctional,
                               ρ::DenseArray{Units.ρ{Cdouble}},
                               ϵ::DenseArray{Units.ϵ{Cdouble}},
                               ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}})
    energy_and_potential!(func, reinterpret(Cdouble, ρ),
                          reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ))
    LDAEnergyPotential(ϵ, ∂ϵ_∂ρ)
end

for name ∈ [:energy_and_potential, :lda, :gga]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(name::Symbol,
                        ρ::Union{DenseArray{Units.ρ{Cdouble}}, DenseArray{Cdouble}},
                        args...)
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name!(name::Symbol, s::Union{Bool, Constants.SPIN},
                        ρ::Union{DenseArray{Units.ρ{Cdouble}}, DenseArray{Cdouble}},
                        args...)
            $name!(XCFunctional(name, s), ρ, args...)
        end
        function $name(name::Symbol,
                       ρ::Union{DenseArray{Units.ρ{Cdouble}}, DenseArray{Cdouble}}, args...)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name(name::Symbol, s::Union{Bool, Constants.SPIN},
                       ρ::Union{DenseArray{Units.ρ{Cdouble}}, DenseArray{Cdouble}}, args...)
            $name(XCFunctional(name, s), ρ, args...)
        end
    end
end
function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Cdouble})
    energy_and_potential!(func, ρ,
                          similar(ρ, eltype(ρ), output_size(func, ρ, 1)), similar(ρ))
end
function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Units.ρ{Cdouble}})
    energy_and_potential!(func, ρ,
                          similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
                          similar(ρ, Units.∂ϵ_∂ρ{Cdouble}))
end

""" Computes the energy and all available derivatives for the given functional """
function lda(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble})
    @check_functional func lda

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

""" Computes the energy and all available derivatives for the given functional """
function lda(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}})
    @check_functional func lda

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        lda!(func, ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
             similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
             similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        lda!(func, ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
             similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        lda!(func, ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)))
    elseif Constants.exc ∈ f
        lda!(func, ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
