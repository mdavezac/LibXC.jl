using AxisArrays, Unitful, DFTShims
using LibXC
const DH = Dispatch.Hartree
const DD = Dispatch.Dimensions

macro test_nothrow(expr)
    quote
        Base.Test.do_test(
            Base.Test.Returned(
                try $(esc(expr)); true; catch false end,
                Base.Test.nothing
            ), $(QuoteNode(expr))
        )
    end
end

@testset "Validating axis arrays" begin
    const valid_array = LibXC.ArrayManips.valid_array
    AXES = Axis{:x}(Base.OneTo(4)), Axis{:y}(Base.OneTo(5)), Axis{:z}(Base.OneTo(6))
    AXES′ = Axis{:xx}(Base.OneTo(4)), AXES[2:end]...
    ρ = rand(typeof(1u"nm^-3"), ColinearSpinLast(), AXES)
    @test_nothrow valid_array(ρ, ρ)
    @test_nothrow valid_array(ρ, zeros(ρ, DH.Scalars.ϵ{Float64}, SpinDegenerate()))
    @test_throws ArgumentError valid_array(ρ, zeros(ρ, SpinDegenerate()))
    @test_throws(ArgumentError, 
                 valid_array(ρ, zeros(eltype(ρ), ColinearSpinLast(), AXES[1:(end - 1)])))
    @test_throws(ArgumentError, 
                 valid_array(ρ, rand(DH.Scalars.∂ϵ_∂ρ{Int64}, ColinearSpinFirst(), AXES′)))
end

@testset "Validating dense arrays" begin
    const valid_array = LibXC.ArrayManips.valid_array
    ρ = zeros(typeof(1u"nm^-3"), 2, 3, 4)
    @test_nothrow valid_array(ρ, ρ, ColinearSpinFirst())
    @test_nothrow valid_array(ρ, ρ, SpinDegenerate())
    @test_nothrow valid_array(ρ, zeros(typeof(1.0u"∂ϵ_∂σ"), size(ρ)), SpinDegenerate())
    n = length(components(typeof(1.0u"∂ϵ_∂σ"), ColinearSpinFirst()))
    @test_nothrow valid_array(ρ, zeros(typeof(1.0u"∂ϵ_∂σ"), n, 3, 4), ColinearSpinFirst())
    @test_throws(ArgumentError,
                 valid_array(ρ, zeros(typeof(1.0u"∂ϵ_∂σ"), size(ρ)), ColinearSpinFirst()))

    @test_nothrow valid_array(ρ, zeros(typeof(1.0u"ϵ"), size(ρ)), SpinDegenerate())
    @test_nothrow valid_array(ρ, zeros(typeof(1.0u"ϵ"), Base.tail(size(ρ))),
                              ColinearSpinFirst())
    @test_throws(ArgumentError,
                 valid_array(ρ, zeros(typeof(1.0u"ϵ"), size(ρ)), ColinearSpinFirst()))
end

@testset "AxisArray to LibXC style" begin
    const to_libxc_array = LibXC.ArrayManips.to_libxc_array
    AXES = Axis{:x}(Base.OneTo(4)), Axis{:y}(Base.OneTo(5)), Axis{:z}(Base.OneTo(6))
    AXES′ = Axis{:xx}(Base.OneTo(4)), AXES[2:end]...

    ρ = rand(typeof(UInt16(1)u"nm^-3"), ColinearSpinLast(), AXES)
    ρ′= @inferred to_libxc_array(ρ, ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ permutedims(ρ, [4, 1, 2, 3]) atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′, ρ′)) === ρ′

    ρ = rand(typeof(Float64(1)u"nm^-3"), ColinearSpinLast(), AXES)
    ρ′= @inferred to_libxc_array(ρ, ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ permutedims(ρ, [4, 1, 2, 3]) atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′, ρ′)) === ρ′

    ρ = rand(typeof(Float64(1)u"ρ"), ColinearSpinLast(), AXES)
    ρ′= to_libxc_array(ρ, ρ)
    @test_broken @inferred to_libxc_array(ρ, ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ permutedims(ρ, [4, 1, 2, 3]) atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′, ρ′)) === ρ′

    ∂ϵ_∂ρ = rand(typeof(UInt16(1)u"eV*nm^3"), ColinearSpin{length(AXES)}(), AXES)
    ∂ϵ_∂ρ′= @inferred to_libxc_array(ρ, ∂ϵ_∂ρ)
    @test eltype(∂ϵ_∂ρ′) === DH.Scalars.∂ϵ_∂ρ{Float64}
    @test ∂ϵ_∂ρ′ ≈ permutedims(∂ϵ_∂ρ, [3, 1, 2, 4]) atol=1e-8u"∂ϵ_∂ρ"
    @test @inferred(to_libxc_array(ρ, ∂ϵ_∂ρ′)) === ∂ϵ_∂ρ′

    ϵ = rand(typeof(UInt16(1)u"eV"), AXES)
    ϵ′= @inferred to_libxc_array(ρ, ϵ)
    @test eltype(ϵ′) === DH.Scalars.ϵ{Float64}
    @test ϵ′ ≈ ϵ atol=1e-8u"ϵ"
    @test to_libxc_array(ρ, ϵ′) === ϵ′
end

@testset "DenseArray to LibXC style" begin
    const to_libxc_array = LibXC.ArrayManips.to_libxc_array
    SIZES = (2, 4, 5)

    ρ = rand(UInt16, SIZES...)u"nm^-3"
    ρ′= @inferred to_libxc_array(ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ ρ atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′)) === ρ′

    ρ = rand(UInt16, SIZES...)u"ρ"
    ρ′= @inferred to_libxc_array(ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ ρ atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′)) === ρ′

    ρ = rand(Float64, SIZES...)u"nm^-3"
    ρ′= @inferred to_libxc_array(ρ)
    @test eltype(ρ′) === DH.Scalars.ρ{Float64}
    @test ρ′ ≈ ρ atol=1e-8u"ρ"
    @test @inferred(to_libxc_array(ρ′)) === ρ′
end
