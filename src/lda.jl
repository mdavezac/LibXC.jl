for (funcname, name, sizer) ∈ [
                                 (:xc_lda_exc, :energy, :esize),
                                 (:xc_lda_vxc, :potential, :vsize),
                                 (:xc_lda_fxc, :second_energy_derivative, :fsize),
                                 (:xc_lda_kxc, :third_energy_derivative, :ksize)]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                        output::DenseArray{Cdouble})
            if family(func) ≠ Constants.lda
                msg = "Incorrect number of arguments: input is not an LDA functional"
                throw(ArgumentError(msg))
            end
            if size(output) ≠ $sizer(func, ρ)
                throw(ArgumentError("sizes of ρ and input are incompatible"))
            end

            ccall(($(parse(":$funcname")), libxc), Void,
                  (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
                  func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, output)
            output
        end

        function $name!(name::Symbol, ρ::DenseArray{Cdouble}, $name::DenseArray{Cdouble})
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, $name)
        end
        function $name!(name::Symbol, s::Union{Constants.SPIN, Bool},
                        ρ::DenseArray{Cdouble}, $name::DenseArray{Cdouble})
            $name!(XCFunctional(name, s), ρ, $name)
        end
        $name(name::Symbol, ρ::DenseArray) = $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ)
        function $name(name::Symbol, s::Union{Bool, Constants.SPIN}, ρ::DenseArray)
            $name(XCFunctional(name, s), ρ)
        end
        function $name(func::AbstractLibXCFunctional, ρ::DenseArray)
            $name!(func, ρ, similar(ρ, eltype(ρ), $sizer(func, ρ)))
        end
    end
end

function lda!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
              εxc::DenseArray{Cdouble}, potential::DenseArray{Cdouble},
              second_deriv::DenseArray{Cdouble}, third_deriv::DenseArray{Cdouble})
    if family(func) ≠ Constants.lda
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end
    if size(εxc) ≠ esize(func, ρ)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(potential) ≠ vsize(func, ρ)
        throw(ArgumentError("sizes of ρ and potential are incompatible"))
    end
    if size(second_deriv) ≠ fsize(func, ρ)
        throw(ArgumentError("sizes of ρ and second derivative are incompatible"))
    end
    if size(third_deriv) ≠ ksize(func, ρ)
        throw(ArgumentError("sizes of ρ and third derivative are incompatible"))
    end

    ccall((:xc_lda, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble},
           Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)),
          ρ, εxc, potential, second_deriv, third_deriv)

    εxc, potential, second_deriv, third_deriv
end

function energy_and_potential!(func::AbstractLibXCFunctional, ρ::DenseArray{Cdouble},
                               εxc::DenseArray{Cdouble}, potential::DenseArray{Cdouble})
    if family(func) ≠ Constants.lda
        msg = "Incorrect number of arguments: input is not an LDA functional"
        throw(ArgumentError(msg))
    end
    if size(εxc) ≠ esize(func, ρ)
        throw(ArgumentError("sizes of ρ and εxc are incompatible"))
    end
    if size(potential) ≠ vsize(func, ρ)
        throw(ArgumentError("sizes of ρ and potential are incompatible"))
    end

    ccall((:xc_lda_exc_vxc, libxc), Void,
          (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          func.c_ptr, length(ρ) / convert(Int64, spin(func)), ρ, εxc, potential)
    εxc, potential
end

for name ∈ [:energy_and_potential, :lda]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(name::Symbol, ρ::DenseArray{Cdouble}, args...)
            $name!(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, args...)
        end
        function $name!(name::Symbol, s::Union{Bool, Constants.SPIN},
                        ρ::DenseArray{Cdouble}, args...)
            $name!(XCFunctional(name, s), ρ, args...)
        end
        function $name(name::Symbol, ρ::DenseArray{Cdouble})
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ)
        end
        function $name(name::Symbol, s::Union{Bool, Constants.SPIN},
                       ρ::DenseArray{Cdouble})
            $name(XCFunctional(name, s), ρ)
        end
    end
end
function energy_and_potential(func::AbstractLibXCFunctional{Cdouble},
                              ρ::DenseArray{Cdouble})
    energy_and_potential!(func, ρ, similar(ρ, eltype(ρ), esize(func, ρ)), similar(ρ))
end
function lda(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble})
    lda!(func, ρ,
         similar(ρ, eltype(ρ), esize(func, ρ)),
         similar(ρ, eltype(ρ), vsize(func, ρ)),
         similar(ρ, eltype(ρ), fsize(func, ρ)),
         similar(ρ, eltype(ρ), ksize(func, ρ)))
end
