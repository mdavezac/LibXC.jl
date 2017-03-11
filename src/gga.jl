"""
    $(SIGNATURES)

GGA energy as a function of ρ and σ=|∇ρ|. The dimensionality is as follows:

    |GGA       | unpolarized | polarized                |
    |----------+-------------+--------------------------|
    |ρ         | any         | `(2, ...)`               |
    |σ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
    |εxc       | `size(ρ)`   | `size(ρ)[2:end]`         |
"""
function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                 σ::DenseArray{Cdouble}, output::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if Constants.exc ∉ flags(func)
        error("Functional does not implement the energy")
    end
    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(output) ≠ output_size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and input are incompatible"))
    end

    ccall((:xc_gga_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, output)
    output
end

function energy(func::AbstractLibXCFunctional, ρ::DenseArray, σ::DenseArray)
    energy!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
end

""" Potential from GGA """
typealias GGAPotential @NT(rho, sigma)

"""
    $(SIGNATURES)

GGA potential computed in place. The dimensionality of the different arrays are as follows:

    |GGA       | unpolarized | polarized                |
    |----------+-------------+--------------------------|
    |ρ         | any         | `(2, ...)`               |
    |σ         | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
    |∂εxc/∂ρ   | `size(ρ)`   | `size(ρ)`                |
    |∂εxc/∂σ   | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
"""
function potential!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                    σ::DenseArray{Cdouble}, pot_rho::DenseArray{Cdouble},
                    pot_sigma::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if Constants.vxc ∉ flags(func)
        error("Functional does not implement first derivatives of the energy")
    end
    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(pot_rho) ≠ size(ρ)
        throw(ArgumentError("sizes of ρ and output pot_rho are incompatible"))
    end
    if size(pot_sigma) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and output pot_sigma are incompatible"))
    end

    ccall((:xc_gga_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, pot_rho, pot_sigma)
    GGAPotential(pot_rho, pot_sigma)
end
function potential(func::AbstractLibXCFunctional, ρ::DenseArray, σ::DenseArray)
    potential!(func, ρ, σ, similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
end

""" Second derivative from GGA

Include the second derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
typealias GGASecondDerivative @NT(rho2, rho_sigma, sigma2)

"""
    $(SIGNATURES)

Second derivatives of GGA energy w.r.t. ρ and σ=|∇ρ|. The dimensionality of the arrays is as
follows:

    |GGA        | unpolarized | polarized                |
    |-----------+-------------+--------------------------|
    |ρ          | any         | `(2, ...)`               |
    |σ          | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
    |∂²εxc/∂ρ²  | `size(ρ)`   | `(3, size(ρ)[2:end]...)` |
    |∂²εxc/∂ρ∂σ | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
    |∂²εxc/∂σ²  | `size(ρ)`   | `(6, size(ρ)[2:end]...)` |
"""
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
                                   ∂ρ²::DenseArray{Cdouble}, ∂ρ∂σ::DenseArray{Cdouble},
                                   ∂σ²::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if Constants.fxc ∉ flags(func)
        error("Functional does not implement second derivatives of the energy")
    end
    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(∂ρ²) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and output ∂²εxc/∂ρ² are incompatible"))
    end
    if size(∂ρ∂σ) ≠ output_size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and output ∂²εxc/∂ρ∂σ are incompatible"))
    end
    if size(∂σ² ) ≠ output_size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and output ∂²εxc/∂σ² are incompatible"))
    end

    ccall((:xc_gga_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ∂ρ²,
          ∂ρ∂σ, ∂σ² )
    GGASecondDerivative(∂ρ², ∂ρ∂σ, ∂σ²)
end
function second_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray,
                                  σ::DenseArray)
    second_energy_derivative!(func, ρ, σ,
                              similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
end

""" Third derivative from GGA

Include the third derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
typealias GGAThirdDerivative @NT(rho3, rho2_sigma, rho_sigma2, sigma3)

"""
    $(SIGNATURES)

Third derivatives of GGA energy w.r.t. ρ and σ=|∇ρ|. The dimensionality of the arrays is as
follows:

    |GGA         | unpolarized | polarized                 |
    |------------+-------------+---------------------------|
    |ρ           | any         | `(2, ...)`                |
    |σ           | `size(ρ)`   | `(3, size(ρ)[2:end]...)`  |
    |∂³εxc/∂ρ³   | `size(ρ)`   | `(4, size(ρ)[2:end]...)`  |
    |∂³εxc/∂ρ²∂σ | `size(ρ)`   | `(9, size(ρ)[2:end]...)`  |
    |∂³εxc/∂ρ∂σ² | `size(ρ)`   | `(10, size(ρ)[2:end]...)` |
    |∂³εxc/∂σ³   | `size(ρ)`   | `(12, size(ρ)[2:end]...)` |
"""
function third_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
                                   ∂ρ³::DenseArray{Cdouble},
                                   ∂ρ²∂σ::DenseArray{Cdouble},
                                   ∂ρ∂σ²::DenseArray{Cdouble},
                                   ∂σ³::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if Constants.kxc ∉ flags(func)
        error("Functional does not implement third derivatives of the energy")
    end
    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(∂ρ³) ≠ output_size(func, ρ, 4)
        throw(ArgumentError("sizes of ρ and output ∂³εxc/∂ρ³ are incompatible"))
    end
    if size(∂ρ²∂σ) ≠ output_size(func, ρ, 9)
        throw(ArgumentError("sizes of ρ and output ∂³εxc/∂ρ²∂σ are incompatible"))
    end
    if size(∂ρ²∂σ) ≠ output_size(func, ρ, 12)
        throw(ArgumentError("sizes of ρ and output ∂³εxc/∂ρ∂σ² are incompatible"))
    end
    if size(deriv_sigma) ≠ output_size(func, ρ, 10)
        throw(ArgumentError("sizes of ρ and output ∂³εxc/∂σ³ are incompatible"))
    end

    ccall((:xc_gga_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, ∂ρ³,
          ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
    GGAThirdDerivative(∂ρ³, ∂ρ²∂σ, ∂ρ∂σ², ∂σ³)
end
function third_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray,
                                  σ::DenseArray)
    second_energy_derivative!(func, ρ, σ,
                              similar(ρ, eltype(ρ), output_size(func, ρ, 4)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 9)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 12)),
                              similar(ρ, eltype(ρ), output_size(func, ρ, 10)))
end

""" Holds GGA energy and first derivatives """
typealias GGAEnergyPotential @NT(energy, pot_rho, pot_sigma)

""" GGA energy and potential """
function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               σ::DenseArray{Cdouble}, εxc::DenseArray{Cdouble},
                               pot_rho::DenseArray{Cdouble}, pot_sigma::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end
    if Constants.exc ∉ flags(func) && Constants.vxc ∉ flags(func)
        error("Functional does not implement energy or potential")
    end
    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(εxc) ≠ output_size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(pot_rho) ≠ size(ρ)
        throw(ArgumentError("sizes of ρ and output pot_rho are incompatible"))
    end
    if size(pot_sigma) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and output pot_sigma are incompatible"))
    end

    ccall((:xc_gga_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, εxc, pot_rho, pot_sigma)
    GGAEnergyPotential(εxc, pot_rho, pot_sigma)
end

function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble})
    energy_and_potential!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
                          similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
end


""" All outputs from LDA """
typealias AllGGA @NT(energy, first_rho, first_sigma, second_rho2, second_rho_sigma,
                     second_sigma2, third_rho3, third_rho2_sigma, third_rho_sigma3,
                     third_sigma3)

"""
    $(SIGNATURES)

GGA energy, first, second, and third derivatives, in place. Arrays for the first, second,
and third derivatives are optional. They should be given only if available for that
particular functional. When requesting higher derivatives, arrays to store the lower
derivatives should also be given.
"""
function gga!{T <: DenseArray{Cdouble}}(func::AbstractLibXCFunctional{Cdouble},
                                        ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
                                        εxc::DenseArray{Cdouble}, outputs::Vararg{T})
    if length(outputs) == 0
        return AllGGA(energy!(func, ρ, σ, εxc), [], [], [], [], [], [], [], [])
    elseif length(outputs) == 2
        result = energy_and_potential!(func, ρ, σ, εxc, outputs...)
        return AllGGA(result[1], result[2], [], [], [], [], [], [], [])
    end
    if length(outputs) ∈ (5, 9) && Constants.fxc ∉ flags(func)
        throw(ArgumentError("Functional does not implement second energy derivative"))
    elseif length(outputs) == 9 && Constants.kxc ∉ flags(func)
        throw(ArgumentError("Functional does not implement third energy derivative"))
    elseif length(outputs) ∉ (5, 9)
        throw(ArgumentError("Incorrect number of outputs, expected 1, 3, 6, or 10"))
    end

    if family(func) ≠ Constants.gga
        throw(ArgumentError("Functional is not a GGA functional"))
    end

    if size(σ) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(εxc) ≠ output_size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end

    if size(outputs[1]) ≠ size(ρ)
        throw(ArgumentError("sizes of ρ and first_rho are incompatible"))
    end
    if size(outputs[2]) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and first_sigma are incompatible"))
    end
    if size(outputs[3]) ≠ output_size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and second_rho2 are incompatible"))
    end
    if size(outputs[4]) ≠ output_size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and secondrho_sigma are incompatible"))
    end
    if size(outputs[5]) ≠ output_size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and second_sigma2 are incompatible"))
    end
    if length(outputs) == 9
        if  size(outputs[6]) ≠ output_size(func, ρ, 4)
            throw(ArgumentError("sizes of ρ and third_rho3 are incompatible"))
        end
        if size(outputs[7]) ≠ output_size(func, ρ, 9)
            throw(ArgumentError("sizes of ρ and third_rho2_sigma are incompatible"))
        end
        if size(outputs[8]) ≠ output_size(func, ρ, 12)
            throw(ArgumentError("sizes of ρ and third_rho_sigma2 are incompatible"))
        end
        if size(outputs[9]) ≠ output_size(func, ρ, 10)
            throw(ArgumentError("sizes of ρ and third_sigma3 are incompatible"))
        end
    end

    args = tuple(outputs..., (C_NULL for i in 1:(9 - length(outputs)))...)

    ccall((:xc_gga, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, σ, εxc, args[1], args[2], args[3], args[4], args[5], args[6], args[7],
          args[8], args[9])

    if length(outputs) == 5
        AllGGA(εxc, outputs..., [], [], [], [])
    else
        AllGGA(εxc, outputs...)
    end
end

""" Computes the energy and all available derivatives for the given functional """
function gga(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
             σ::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        throw(ArgumentError("Functional is not a GGA functional"))
    end

    const f = flags(func)
    if Constants.exc ∈ f && Constants.vxc ∈ f && Constants.fxc ∈ f && Constants.kxc ∈ f
        gga!(func, ρ, σ,
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
        gga!(func, ρ, σ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 3)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)),
             similar(ρ, eltype(ρ), output_size(func, ρ, 6)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f && Constants.vxc ∈ f
        gga!(func, ρ, σ,
             similar(ρ, eltype(ρ), output_size(func, ρ, 1)),
             similar(ρ), similar(ρ, eltype(ρ), output_size(func, ρ, 3)))
    elseif Constants.exc ∈ f
        gga!(func, ρ, σ, similar(ρ, eltype(ρ), output_size(func, ρ, 1)))
    else
        throw(ArgumentError("Not sure what this functional can do"))
    end
end
