using LibXC: OutputTuples
using DFTShims: ColinearSpinFirst, Dispatch
const DH = Dispatch.Hartree
@testset "indexing and length" begin
    @test length(OutputTuples.AllLDA(1, 2, 3, 4)) == 4
    @test OutputTuples.AllLDA(1, 2, 3, 4)[1] == 1
    @test OutputTuples.AllLDA(1, 2, 3, 4)[2] == 2
    @test OutputTuples.AllLDA(1, 2, 3, 4)[3] == 3
    @test OutputTuples.AllLDA(1, 2, 3, 4)[4] == 4

    ∂ϵ_∂ρ = zeros(DH.Scalars.∂ϵ_∂ρ{Float64}, ColinearSpinFirst(), (2, 3))
    ∂³ϵ_∂ρ³ = zeros(DH.Scalars.∂³ϵ_∂ρ³{Float64}, ColinearSpinFirst(), (2, 3))

    @test OutputTuples.LDATuple(∂³ϵ_∂ρ³, ∂ϵ_∂ρ).ϵ === nothing
    @test OutputTuples.LDATuple(∂³ϵ_∂ρ³, ∂ϵ_∂ρ).∂ϵ_∂ρ === ∂ϵ_∂ρ
    @test OutputTuples.LDATuple(∂³ϵ_∂ρ³, ∂ϵ_∂ρ).∂²ϵ_∂ρ² === nothing
    @test OutputTuples.LDATuple(∂³ϵ_∂ρ³, ∂ϵ_∂ρ).∂³ϵ_∂ρ³ === ∂³ϵ_∂ρ³
end
