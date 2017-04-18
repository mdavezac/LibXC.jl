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

for name ∈ [:energy, :potential, :second_energy_derivative, :third_energy_derivative, :gga,
            :energy_and_potential]
    @eval function $name{T <: Real, TT <: Real}(func::AbstractLibXCFunctional{Cdouble},
                                                ρ::DenseArray{T}, σ::DenseArray{TT})
        $name(func, convert(Array{Cdouble}, ρ), convert(Array{Cdouble}, σ))
    end
end
