module Fixtures
export expected_data, input_data

using LibXC
using DFTShims
using DataFrames
using AxisArrays

const DHS = Dispatch.Hartree.Scalars

nospin(t, a) = reinterpret(t, SpinDegenerate(), a)
withspin(t, a) = reinterpret(t, ColinearSpinFirst(), a)

input_data(name::String) = begin
    data_file = joinpath(dirname(Base.source_path()), "data", name)

    input = (map(x -> parse(Float64, x), split(v))
             for v in readlines(data_file)[2:end])
    input = transpose(hcat(input...))
    input = DataFrame(Any[input[:, i] for i in 1:size(input, 2)],
                      [:ρ_a, :ρ_b, :σ_aa, :σ_ab, :σ_bb, :δ_a, :δ_b, :τ_a, :τ_b])
end

expected_data(name::String) = begin
    data_file = joinpath(dirname(Base.source_path()), "data", name)
    expected = (map(x -> parse(Float64, x), split(v))
                for v in readlines(data_file)[3:end])
    expected = transpose(hcat(expected...))
    if size(expected, 2) == 3
        DataFrame(Any[nospin(DHS.ϵ{Cdouble}, expected[:, 1]),
                      nospin(DHS.∂ϵ_∂ρ{Cdouble}, expected[:, 2]),
                      nospin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 3])],
                  [:ε, :v, :δv])
    elseif size(expected, 2) == 6
        DataFrame(Any[withspin(DHS.ϵ{Cdouble}, expected[:, 1]),
                      withspin(DHS.∂ϵ_∂ρ{Cdouble}, expected[:, 2]),
                      withspin(DHS.∂ϵ_∂ρ{Cdouble}, expected[:, 3]),
                      withspin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 4]),
                      withspin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 5]),
                      withspin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 6])],
                  [:ε, :v_a, :v_b, :δv_aa, :δv_ab, :δv_bb])
    elseif size(expected, 2) > 6
        DataFrame(Any[nospin(DHS.ϵ{Cdouble}, expected[:, 1]),
                      nospin(DHS.∂ϵ_∂ρ{Cdouble}, expected[:, 2]),
                      nospin(DHS.∂ϵ_∂ρ{Cdouble}, expected[:, 3]),
                      nospin(DHS.∂ϵ_∂σ{Cdouble}, expected[:, 4]),
                      nospin(DHS.∂ϵ_∂σ{Cdouble}, expected[:, 5]),
                      nospin(DHS.∂ϵ_∂σ{Cdouble}, expected[:, 6]),
                      nospin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 7]),
                      nospin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 8]),
                      nospin(DHS.∂²ϵ_∂ρ²{Cdouble}, expected[:, 9]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 10]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 11]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 12]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 13]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 14]),
                      nospin(DHS.∂²ϵ_∂σ²{Cdouble}, expected[:, 15]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 16]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 17]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 18]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 19]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 20]),
                      nospin(DHS.∂²ϵ_∂ρ∂σ{Cdouble}, expected[:, 21])],
                   [:ε, :vrho_a, :vrho_b, :vsigma_aa, :vsigma_ab, :vsigma_bb, :v2rho_aa,
                    :v2rho_ab, :v2rho_bb, :v2sigma2_aa_aa, :v2sigma2_aa_ab, :v2sigma2_aa_bb,
                    :v2sigma2_ab_ab, :v2sigma2_ab_bb, :v2sigma2_bb_bb, :v2rho_asigma_aa,
                    :v2rho_asigma_ab, :v2rho_asigma_bb, :v2rho_bsigma_aa, :v2rho_bsigma_ab,
                    :v2rho_bsigma_bb])
    end
end
end
