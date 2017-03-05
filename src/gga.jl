""" GGA energy as a function of ρ and σ=|∇ρ| """
function energy!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                 σ::DenseArray{Cdouble}, output::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if size(output) ≠ size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and input are incompatible"))
    end
    if size(σ) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end

    ccall((:xc_gga_exc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, output)
    output
end

function energy(func::AbstractLibXCFunctional, ρ::DenseArray, σ::DenseArray)
    energy!(func, ρ, σ, similar(ρ, eltype(ρ), esize(func, ρ)))
end

""" Potential from GGA """
typealias GGAPotential @NT(rho, sigma)

""" First derivatives of GGA energy w.r.t. ρ and σ=|∇ρ| """
function potential!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                    σ::DenseArray{Cdouble}, pot_rho::DenseArray{Cdouble},
                    pot_sigma::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if size(σ) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(pot_rho) ≠ size(ρ)
        throw(ArgumentError("sizes of ρ and output pot_rho are incompatible"))
    end
    if size(pot_sigma) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and output pot_sigma are incompatible"))
    end

    ccall((:xc_gga_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, pot_rho, pot_sigma)
    GGAPotential(pot_rho, pot_sigma)
end
function potential(func::AbstractLibXCFunctional, ρ::DenseArray, σ::DenseArray)
    potential!(func, ρ, σ, similar(ρ), similar(ρ, eltype(ρ), size(func, ρ, 3)))
end

""" Second derivative from GGA

Include the second derivative of the energy with respect to ρ, σ, and both ρ and σ.
"""
typealias GGASecondDerivative @NT(rho2, rho_sigma, sigma2)

""" Second derivatives of GGA energy w.r.t. ρ and σ=|∇ρ| """
function second_energy_derivative!(func::AbstractLibXCFunctional{Cdouble},
                                   ρ::DenseArray{Cdouble}, σ::DenseArray{Cdouble},
                                   deriv_rho::DenseArray{Cdouble},
                                   deriv_rho_sigma::DenseArray{Cdouble},
                                   deriv_sigma::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an GGA functional"
        throw(ArgumentError(msg))
    end
    if size(σ) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(deriv_rho) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and output deriv_rho are incompatible"))
    end
    if size(deriv_rho_sigma) ≠ size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and output deriv_rho_sigma are incompatible"))
    end
    if size(deriv_sigma) ≠ size(func, ρ, 6)
        throw(ArgumentError("sizes of ρ and output deriv_sigma are incompatible"))
    end

    ccall((:xc_gga_fxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, σ, deriv_rho,
          deriv_rho_sigma, deriv_sigma)
    GGASecondDerivative(deriv_rho, deriv_rho_sigma, deriv_sigma)
end
function second_energy_derivative(func::AbstractLibXCFunctional, ρ::DenseArray,
                                  σ::DenseArray)
    second_energy_derivative!(func, ρ, σ,
                              similar(ρ, eltype(ρ), size(func, ρ, 3)),
                              similar(ρ, eltype(ρ), size(func, ρ, 6)),
                              similar(ρ, eltype(ρ), size(func, ρ, 6)))
end

""" Holds GGA energy and first derivatives """
typealias GGAEnergyPotential @NT(energy, pot_rho, pot_sigma)

""" GGA energy and first derivatives """
function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               σ::DenseArray{Cdouble}, εxc::DenseArray{Cdouble},
                               pot_rho::DenseArray{Cdouble}, pot_sigma::DenseArray{Cdouble})
    if family(func) ≠ Constants.gga
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end
    if size(σ) ≠ size(func, ρ, 3)
        throw(ArgumentError("sizes of ρ and σ are incompatible"))
    end
    if size(εxc) ≠ size(func, ρ, 1)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(pot_rho) ≠ size(ρ)
        throw(ArgumentError("sizes of ρ and output pot_rho are incompatible"))
    end
    if size(pot_sigma) ≠ size(func, ρ, 3)
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
    energy_and_potential!(func, ρ, σ, similar(ρ, eltype(ρ), size(func, ρ, 1)),
                          similar(ρ), similar(ρ, eltype(ρ), size(func, ρ, 3)))
end
