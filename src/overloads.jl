for name ∈ [:energy, :potential,:second_energy_derivative, :third_energy_derivative, :lda,
            :energy_and_potential]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name{T <: Number}(func::AbstractLibXCFunctional{Cdouble},
                                    ρ::DenseArray{T})
            $name(func, Units.conversion(Units.ρ, ρ))
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
        function $name(name::Symbol, ρ::DenseArray)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ)
        end

        @doc """
            $(SIGNATURES)

        A simple way to call a functional using it's name. Spin polarization is explicitly
        specified.
        """ ->
        function $name(name::Symbol, spin::Union{Bool, Constants.SPIN}, ρ::DenseArray)
            $name(XCFunctional(name, spin), ρ)
        end
    end
end

for name ∈ [:energy, :potential, :second_energy_derivative, :third_energy_derivative, :gga,
            :energy_and_potential]
    @eval begin
        function $name{T <: Real, TT <: Real}(func::AbstractLibXCFunctional{Cdouble},
                                              ρ::DenseArray{T}, ∇ρ::DenseArray{TT})
            $name(func, Units.conversion(Units.ρ, ρ), Units.conversion(Units.∇ρ, ∇ρ))
        end
        function $name(name::Symbol, ρ::DenseArray, ∇ρ::DenseArray)
            $name(name, ndims(ρ) > 1 && size(ρ, 1) == 2, ρ, ∇ρ)
        end
        function $name(name::Symbol, spin::Union{Bool, Constants.SPIN}, ρ::DenseArray,
                       ∇ρ::DenseArray)
            $name(XCFunctional(name, spin), ρ, ∇ρ)
        end
    end
end
