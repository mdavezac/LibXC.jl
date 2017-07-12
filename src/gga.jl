# Mostly so linting doesn't throw false positives when processing this file alone
using LibXC: Constants, CFuncType, libxc, spin, output_size, flags
using LibXC: GGAPotential, GGASecondDerivative, GGAThirdDerivative
using LibXC: GGAEnergyAndPotential, AllGGA
using LibXC.Checks: @check_functional, @check_availability, @check_size
@lintpragma("Ignore use of undeclared variable ccall")

"""
    $(SIGNATURES)

GGA energy as a function of ρ and ∇ρ=|∇ρ|. The dimensionality is as follows:

|GGA | unpolarized | polarized                |
|----|-------------|--------------------------|
|ρ   | any         | `(2, ...)`               |
|∇ρ  | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
|ϵ   | `size(ρ)`   | `size(ρ)[2:end]`         |
"""
function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}},
                 ∇ρ::DenseArray{Units.∇ρ{Cdouble}}, ϵ::DenseArray{Units.ϵ{Cdouble}})
    energy!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ), reinterpret(Cdouble, ϵ))
    ϵ
end
function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                 ∇ρ::DenseArray{Cdouble}, ϵ::DenseArray{Cdouble})
    @check_functional func gga
    @check_availability func exc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ϵ 1

    ccall((:xc_gga_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) /convert(Int64, spin(func)), ρ, ∇ρ, ϵ)
    ϵ
end
function energy{Ρ <: Quantity, ∇Σ <: Quantity}(func::AbstractLibXCFunctional,
                                               ρ::DenseArray{Ρ}, ∇ρ::DenseArray{∇Σ})
    energy!(func, Units.conversion(Units.ρ, ρ), Units.conversion(Units.∇ρ, ∇ρ),
            similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
end
function energy(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                ∇ρ::DenseArray{Cdouble})
    energy!(func, ρ, ∇ρ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
end


"""
    $(SIGNATURES)

GGA potential computed in place. The dimensionality of the different arrays are as follows:

|GGA     | unpolarized | polarized                |
|--------|-------------|--------------------------|
|ρ       | any         | `(2, ...)`               |
|∇ρ      | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
|∂ϵ/∂ρ   | `size(ρ)`   | `size(ρ)`                |
|∂ϵ/∂∇ρ  | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
"""
function potential!(func::AbstractLibXCFunctional{Cdouble},
                    ρ::DenseArray{Units.ρ{Cdouble}},
                    ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
                    ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
                    ∂ϵ_∂∇ρ::DenseArray{Units.∂ϵ_∂∇ρ{Cdouble}})
    potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
               reinterpret(Cdouble, ∂ϵ_∂ρ), reinterpret(Cdouble, ∂ϵ_∂∇ρ))
    GGAPotential(∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
end
function potential!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                    ∇ρ::DenseArray{Cdouble}, ∂ϵ_∂ρ::DenseArray{Cdouble},
                    ∂ϵ_∂∇ρ::DenseArray{Cdouble})
    @check_functional func gga
    @check_availability func vxc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ∂ϵ_∂ρ 2
    @check_size func ρ ∂ϵ_∂∇ρ 3

    ccall((:xc_gga_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, ∇ρ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
    GGAPotential(∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
end
function potential{Ρ <: Quantity, ∇Σ <: Quantity}(func::AbstractLibXCFunctional,
                                                  ρ::DenseArray{Ρ},
                                                  ∇ρ::DenseArray{∇Σ})
    potential!(func, Units.conversion(Units.ρ, ρ), Units.conversion(Units.∇ρ, ∇ρ),
               similar(ρ, Units.∂ϵ_∂ρ{Cdouble}),
               similar(ρ, Units.∂ϵ_∂∇ρ{Cdouble}, output_size(func, ρ, 3)))
end
function potential(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                   ∇ρ::DenseArray{Cdouble})
    potential!(func, ρ, ∇ρ, similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
end

"""
    $(SIGNATURES)

Second derivatives of GGA energy w.r.t. ρ and ∇ρ=|∇ρ|. The dimensionality of the arrays is as
follows:

|GGA       | unpolarized | polarized                |
|----------|-------------|--------------------------|
|ρ         | any         | `(2, ...)`               |
|∇ρ        | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
|∂²ϵ/∂ρ²   | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
|∂²ϵ/∂ρ∂∇ρ | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
|∂²ϵ/∂∇ρ²  | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
"""
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Units.ρ{Cdouble}},
                                   ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
                                   ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
                                   ∂²ϵ_∂ρ∂∇ρ::DenseArray{Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}},
                                   ∂²ϵ_∂∇ρ²::DenseArray{Units.∂²ϵ_∂∇ρ²{Cdouble}})
    second_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                              reinterpret(Cdouble, ∂²ϵ_∂ρ²), reinterpret(Cdouble, ∂²ϵ_∂ρ∂∇ρ),
                              reinterpret(Cdouble, ∂²ϵ_∂∇ρ²))
    GGASecondDerivative(∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ²)
end
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble}, ∇ρ::DenseArray{Cdouble},
                                   ∂²ϵ_∂ρ²::DenseArray{Cdouble},
                                   ∂²ϵ_∂ρ∂∇ρ::DenseArray{Cdouble},
                                   ∂²ϵ_∂∇ρ²::DenseArray{Cdouble})
    @check_functional func gga
    @check_availability func fxc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ∂²ϵ_∂ρ² 3
    @check_size func ρ ∂²ϵ_∂ρ∂∇ρ 6
    @check_size func ρ ∂²ϵ_∂∇ρ² 6

    ccall((:xc_gga_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, ∇ρ, ∂²ϵ_∂ρ²,
          ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ² )
    GGASecondDerivative(∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ²)
end
function second_energy_derivative{Ρ <: Quantity, ∇Σ <: Quantity}(
                    func::AbstractLibXCFunctional, ρ::DenseArray{Ρ}, ∇ρ::DenseArray{∇Σ})
    second_energy_derivative!(func,
                              Units.conversion(Units.ρ, ρ),
                              Units.conversion(Units.∇ρ, ∇ρ),
                              similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
                              similar(ρ, Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}, output_size(func, ρ, 6)),
                              similar(ρ, Units.∂²ϵ_∂∇ρ²{Cdouble}, output_size(func, ρ, 6)))
end
function second_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                                  ∇ρ::DenseArray{Cdouble})
    second_energy_derivative!(func, ρ, ∇ρ,
                              similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
end

"""
    $(SIGNATURES)

Third derivatives of GGA energy w.r.t. ρ and ∇ρ=|∇ρ|. The dimensionality of the arrays is as
follows:

|GGA        | unpolarized | polarized                 |
|-----------|-------------|---------------------------|
|ρ          | any         | `(2, ...)`                |
|∇ρ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
|∂³ϵ/∂ρ³    | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
|∂³ϵ/∂ρ²∂∇ρ | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
|∂³ϵ/∂ρ∂∇ρ² | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
|∂³ϵ/∂∇ρ³   | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |
"""
function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                  ρ::DenseArray{Units.ρ{Cdouble}},
                                  ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
                                  ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}},
                                  ∂³ϵ_∂ρ²∂∇ρ::DenseArray{Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}},
                                  ∂³ϵ_∂ρ∂∇ρ²::DenseArray{Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}},
                                  ∂³ϵ_∂∇ρ³::DenseArray{Units.∂³ϵ_∂∇ρ³{Cdouble}})
    third_energy_derivative!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∂³ϵ_∂ρ³),
                             reinterpret(Cdouble, ∂³ϵ_∂ρ²∂∇ρ),
                             reinterpret(Cdouble, ∂³ϵ_∂ρ∂∇ρ²),
                             reinterpret(Cdouble, ∂³ϵ_∂∇ρ³))
    GGAThirdDerivative(∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂∇ρ, ∂³ϵ_∂ρ∂∇ρ², ∂³ϵ_∂∇ρ³)
end
function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                  ρ::DenseArray{Cdouble}, ∇ρ::DenseArray{Cdouble},
                                  ∂³ϵ_∂ρ³::DenseArray{Cdouble},
                                  ∂³ϵ_∂ρ²∂∇ρ::DenseArray{Cdouble},
                                  ∂³ϵ_∂ρ∂∇ρ²::DenseArray{Cdouble},
                                  ∂³ϵ_∂∇ρ³::DenseArray{Cdouble})
    @check_functional func gga
    @check_availability func kxc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ∂³ϵ_∂ρ³ 4
    @check_size func ρ ∂³ϵ_∂ρ²∂∇ρ 9
    @check_size func ρ ∂³ϵ_∂ρ∂∇ρ² 12
    @check_size func ρ ∂³ϵ_∂∇ρ³ 10

    ccall((:xc_gga_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, ∇ρ, ∂³ϵ_∂ρ³,
          ∂³ϵ_∂ρ²∂∇ρ, ∂³ϵ_∂ρ∂∇ρ², ∂³ϵ_∂∇ρ³)
    GGAThirdDerivative(∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂∇ρ, ∂³ϵ_∂ρ∂∇ρ², ∂³ϵ_∂∇ρ³)
end
function third_energy_derivative{Ρ <: Quantity, ∇Σ <: Quantity}(
                func::AbstractLibXCFunctional, ρ::DenseArray{Ρ}, ∇ρ::DenseArray{∇Σ})
    third_energy_derivative!(func,
                             Units.conversion(Units.ρ, ρ), Units.conversion(Units.∇ρ, ∇ρ),
                             similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)),
                             similar(ρ, Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}, output_size(func, ρ, 9)),
                             similar(ρ, Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}, output_size(func, ρ, 12)),
                             similar(ρ, Units.∂³ϵ_∂∇ρ³{Cdouble}, output_size(func, ρ, 10)))
end
function third_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                                 ∇ρ::DenseArray{Cdouble})
    second_energy_derivative!(func, ρ, ∇ρ,
                              similar(ρ, eltype(ρ), output_size(func, ρ, 4)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 9)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 12)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 10)))
end

"""
    $(SIGNATURES)

GGA energy and potential
"""
function energy_and_potential!(func::AbstractLibXCFunctional{Cdouble},
                               ρ::DenseArray{Units.ρ{Cdouble}},
                               ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
                               ϵ::DenseArray{Units.ϵ{Cdouble}},
                               ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
                               ∂ϵ_∂∇ρ::DenseArray{Units.∂ϵ_∂∇ρ{Cdouble}})
    energy_and_potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                          reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
                          reinterpret(Cdouble, ∂ϵ_∂∇ρ))
    GGAEnergyAndPotential(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
end
function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               ∇ρ::DenseArray{Cdouble}, ϵ::DenseArray{Cdouble},
                               ∂ϵ_∂ρ::DenseArray{Cdouble}, ∂ϵ_∂∇ρ::DenseArray{Cdouble})
    @check_functional func gga
    @check_availability func vxc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ϵ 1
    @check_size func ρ ∂ϵ_∂ρ 2
    @check_size func ρ ∂ϵ_∂∇ρ 3

    ccall((:xc_gga_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, ∇ρ, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
    GGAEnergyAndPotential(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)
end
function energy_and_potential(func::AbstractLibXCFunctional, ρ::DenseArray{Units.ρ{Cdouble}},
                ∇ρ::DenseArray{Units.∇ρ{Cdouble}})
    energy_and_potential!(func, ρ, ∇ρ, similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
                          similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
                          similar(ρ, Units.∂ϵ_∂∇ρ{Cdouble}, output_size(func, ρ, 3)))
end
function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Cdouble}, ∇ρ::DenseArray{Cdouble})
    energy_and_potential!(func, ρ, ∇ρ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
                          similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
end


"""
    $(SIGNATURES)

GGA energy, first, second, and third derivatives, in place. Arrays for the first, second,
and third derivatives are optional. They should be given only if available for that
particular functional. When requesting higher derivatives, arrays to store the lower
derivatives should also be given.
"""
function gga!{T <: DenseArray{Cdouble}}(func::AbstractLibXCFunctional{Cdouble},
                                        ρ::DenseArray{Cdouble}, ∇ρ::DenseArray{Cdouble},
                                        ϵ::DenseArray{Cdouble}, outputs::Vararg{T})
    if length(outputs) == 0
        return AllGGA(energy!(func, ρ, ∇ρ, ϵ), [], [], [], [], [], [], [], [])
    elseif length(outputs) == 2
        result = energy_and_potential!(func, ρ, ∇ρ, ϵ, outputs...)
        return AllGGA(result[1], result[2], [], [], [], [], [], [], [])
    end

    @check_functional func gga
    @check_availability func exc
    @check_availability func vxc
    length(outputs) ∉ (5, 9) && throw(ArgumentError("Incorrect number of arguments"))
    length(outputs) ∈ (5, 9) && @check_availability func fxc
    length(outputs) == 9 && @check_availability func kxc
    @check_size func ρ ∇ρ 3
    @check_size func ρ ϵ 1
    @check_size func ρ outputs[1] 2
    @check_size func ρ outputs[2] 3
    @check_size func ρ outputs[3] 3
    @check_size func ρ outputs[4] 6
    @check_size func ρ outputs[5] 6
    if length(outputs) == 9
        @check_size func ρ outputs[6] 4
        @check_size func ρ outputs[7] 9
        @check_size func ρ outputs[8] 12
        @check_size func ρ outputs[9] 10
    end

    args = tuple(outputs..., (C_NULL for i in 1:(9 - length(outputs)))...)

    ccall((:xc_gga, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, ∇ρ, ϵ, args[1], args[2], args[3], args[4], args[5], args[6], args[7],
          args[8], args[9])

    if length(outputs) == 5
        AllGGA(ϵ, outputs..., [], [], [], [])
    else
        AllGGA(ϵ, outputs...)
    end
end
function gga!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}},
              ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
              ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
              ∂ϵ_∂∇ρ::DenseArray{Units.∂ϵ_∂∇ρ{Cdouble}},
              ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
              ∂²ϵ_∂ρ∂∇ρ::DenseArray{Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}},
              ∂²ϵ_∂∇ρ²::DenseArray{Units.∂²ϵ_∂∇ρ²{Cdouble}},
              ∂³ϵ_∂ρ³::DenseArray{Units.∂³ϵ_∂ρ³{Cdouble}},
              ∂³ϵ_∂ρ²∂∇ρ::DenseArray{Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}},
              ∂³ϵ_∂ρ∂∇ρ²::DenseArray{Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}},
              ∂³ϵ_∂∇ρ³::DenseArray{Units.∂³ϵ_∂∇ρ³{Cdouble}})
    result = gga!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                  reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
                  reinterpret(Cdouble, ∂ϵ_∂∇ρ), reinterpret(Cdouble, ∂²ϵ_∂ρ²),
                  reinterpret(Cdouble, ∂²ϵ_∂ρ∂∇ρ), reinterpret(Cdouble, ∂²ϵ_∂∇ρ²),
                  reinterpret(Cdouble, ∂³ϵ_∂ρ³), reinterpret(Cdouble, ∂³ϵ_∂ρ²∂∇ρ),
                  reinterpret(Cdouble, ∂³ϵ_∂ρ∂∇ρ²), reinterpret(Cdouble, ∂³ϵ_∂∇ρ³))
    AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ², ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂∇ρ, ∂³ϵ_∂ρ∂∇ρ²,
           ∂³ϵ_∂∇ρ³)
end
function gga!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
              ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
              ∂ϵ_∂∇ρ::DenseArray{Units.∂ϵ_∂∇ρ{Cdouble}},
              ∂²ϵ_∂ρ²::DenseArray{Units.∂²ϵ_∂ρ²{Cdouble}},
              ∂²ϵ_∂ρ∂∇ρ::DenseArray{Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}},
              ∂²ϵ_∂∇ρ²::DenseArray{Units.∂²ϵ_∂∇ρ²{Cdouble}})
    result = gga!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                  reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
                  reinterpret(Cdouble, ∂ϵ_∂∇ρ), reinterpret(Cdouble, ∂²ϵ_∂ρ²),
                  reinterpret(Cdouble, ∂²ϵ_∂ρ∂∇ρ), reinterpret(Cdouble, ∂²ϵ_∂∇ρ²))
    AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ², Units.∂³ϵ_∂ρ³{Cdouble}[],
           Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}[], Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}[],
           Units.∂³ϵ_∂∇ρ³{Cdouble}[])
end
function gga!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
              ϵ::DenseArray{Units.ϵ{Cdouble}},
              ∂ϵ_∂ρ::DenseArray{Units.∂ϵ_∂ρ{Cdouble}},
              ∂ϵ_∂∇ρ::DenseArray{Units.∂ϵ_∂∇ρ{Cdouble}})
    result = energy_and_potential!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                                   reinterpret(Cdouble, ϵ), reinterpret(Cdouble, ∂ϵ_∂ρ),
                                   reinterpret(Cdouble, ∂ϵ_∂∇ρ))
    AllGGA(ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ, Units.∂²ϵ_∂ρ²{Cdouble}[], Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}[],
           Units.∂²ϵ_∂∇ρ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[], Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}[],
           Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}[], Units.∂³ϵ_∂∇ρ³{Cdouble}[])
end
function gga!(func::AbstractLibXCFunctional{Cdouble},
              ρ::DenseArray{Units.ρ{Cdouble}}, ∇ρ::DenseArray{Units.∇ρ{Cdouble}},
              ϵ::DenseArray{Units.ϵ{Cdouble}})
    result = energy!(func, reinterpret(Cdouble, ρ), reinterpret(Cdouble, ∇ρ),
                     reinterpret(Cdouble, ϵ))
    AllGGA(ϵ, Units.∂ϵ_∂ρ{Cdouble}[], Units.∂ϵ_∂∇ρ{Cdouble}[], Units.∂²ϵ_∂ρ²{Cdouble}[],
           Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}[], Units.∂²ϵ_∂∇ρ²{Cdouble}[], Units.∂³ϵ_∂ρ³{Cdouble}[],
           Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}[], Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}[], Units.∂³ϵ_∂∇ρ³{Cdouble}[])
end


"""
    $(SIGNATURES)

Computes the energy and all available derivatives for the given functional
"""
function gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Units.ρ{Cdouble}},
             ∇ρ::DenseArray{Units.∇ρ{Cdouble}})
    @check_functional func gga
    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
             similar(ρ, Units.∂ϵ_∂∇ρ{Cdouble}, output_size(func, ρ, 3)),
             similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
             similar(ρ, Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}, output_size(func, ρ, 6)),
             similar(ρ, Units.∂²ϵ_∂∇ρ²{Cdouble}, output_size(func, ρ, 6)),
             similar(ρ, Units.∂³ϵ_∂ρ³{Cdouble}, output_size(func, ρ, 4)),
             similar(ρ, Units.∂³ϵ_∂ρ²∂∇ρ{Cdouble}, output_size(func, ρ, 9)),
             similar(ρ, Units.∂³ϵ_∂ρ∂∇ρ²{Cdouble}, output_size(func, ρ, 10)),
             similar(ρ, Units.∂³ϵ_∂∇ρ³{Cdouble}, output_size(func, ρ, 12)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
             similar(ρ, Units.∂ϵ_∂∇ρ{Cdouble}, output_size(func, ρ, 3)),
             similar(ρ, Units.∂²ϵ_∂ρ²{Cdouble}, output_size(func, ρ, 3)),
             similar(ρ, Units.∂²ϵ_∂ρ∂∇ρ{Cdouble}, output_size(func, ρ, 6)),
             similar(ρ, Units.∂²ϵ_∂∇ρ²{Cdouble}, output_size(func, ρ, 6)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)),
             similar(ρ, Units.∂ϵ_∂ρ{Cdouble}, output_size(func, ρ, 2)),
             similar(ρ, Units.∂ϵ_∂∇ρ{Cdouble}, output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, Units.ϵ{Cdouble}, output_size(func, ρ, 1)))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
function gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
             ∇ρ::DenseArray{Cdouble})
    @check_functional func gga

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 4)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 9)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 12)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 10)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, ∇ρ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f
        gga!(func, ρ, ∇ρ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
