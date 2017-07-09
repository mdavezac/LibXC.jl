var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "CurrentModule = LibXC\nDocTestSetup = quote\n    using LibXC\n    using Unitful\n    using UnitfulHartree\n    func = XCFunctional(:lda_x, false)\nend"
},

{
    "location": "index.html#LibXC-1",
    "page": "Home",
    "title": "LibXC",
    "category": "section",
    "text": "This package brings julia bindings for libxc. Currently, only bindings for LDA and GGA functionals exist. Bindings for the routines returning information about the functional, e.g. a longer description or citations are also given."
},

{
    "location": "index.html#LibXC.libxc_functionals",
    "page": "Home",
    "title": "LibXC.libxc_functionals",
    "category": "Function",
    "text": "Names of the available libxc functionals \n\n\n\n"
},

{
    "location": "index.html#Creating-a-Functional-1",
    "page": "Home",
    "title": "Creating a Functional",
    "category": "section",
    "text": "The most efficient way to access the functionals is to create one:julia> functional = XCFunctional(:lda_x, true)\nLibXC.XCFunctional{Float64}(:lda_x, polarized)\n  - description: Slater exchange\n  - kind: exchange\n  - family: lda\n  - spin: polarized\n  - citations:\n      * P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)\n      * F. Bloch, Z. Phys. 57, 545 (1929)The first argument is the name of the functional, and the second (optional, defaults to true) is whether the functional is spin-polarized. The names of the available functionals can be obtained withlibxc_functionalsThe functionals can be queried for their kind, family, description, citations, and spin."
},

{
    "location": "index.html#A-word-about-input-dimensionality-1",
    "page": "Home",
    "title": "A word about input dimensionality",
    "category": "section",
    "text": "The functionals expect input arrays ρ and ∇ρ, and (optionally) a number of output arrays, for the energy ϵ, and the different derivatives, e.g. ∂ϵ/∂ρ. Because we are accessing a C library, some care must be taken when creating these arrays.All arrays must be dense (contiguous in memory) and the element type must match Cdouble (On most systems, Cdouble is an alias for Float64).\nnon-spin-polarized cases: ρ can have any dimensions. All other arrays must match ρ.\nspin-polarized cases: the first dimension of ρ must be 2: size(ρ) = (2, ....). Other arrays must match in very specific ways:\nLDA unpolarized polarized\nρ any (2, ...)\nϵ size(ρ) size(ρ)[2:end]\n∂ϵ/∂ρ size(ρ) size(ρ)\n∂²ϵ/∂ρ² size(ρ) (3, size(ρ)[2:end]...)\n∂³ϵ/∂ρ³ size(ρ) (4, size(ρ)[2:end]...)\nGGA unpolarized polarized\nρ any (2, ...)\n∇ρ size(ρ) (3, size(ρ)[2:end]...)\nϵ size(ρ) size(ρ)[2:end]\n∂ϵ/∂ρ size(ρ) size(ρ)\n∂ϵ/∂∇ρ size(ρ) (3, size(ρ)[2:end]...)\n∂²ϵ/∂ρ² size(ρ) (3, size(ρ)[2:end]...)\n∂²ϵ/∂ρ∂∇ρ size(ρ) (6, size(ρ)[2:end]...)\n∂²ϵ/∂∇ρ² size(ρ) (6, size(ρ)[2:end]...)\n∂³ϵ/∂ρ³ size(ρ) (4, size(ρ)[2:end]...)\n∂³ϵ/∂ρ²∂∇ρ size(ρ) (9, size(ρ)[2:end]...)\n∂³ϵ/∂ρ∂∇ρ² size(ρ) (10, size(ρ)[2:end]...)\n∂³ϵ/∂∇ρ³ size(ρ) (12, size(ρ)[2:end]...)For the exact meaning of each dimension in each array, please refer to libxc"
},

{
    "location": "index.html#A-word-about-physical-units-1",
    "page": "Home",
    "title": "A word about physical units",
    "category": "section",
    "text": "The underlying C library expects inputs in Hartree atomic units. It is possible (and recommended) to make units part of the type of the inputs and outputs. When using Hartree atomic units with Cdouble, as shown below, this will not incur any overhead. We use Unitful, UnitfulHartree, and to defined (within LibXC.DFTUnits) a set of units to represent the electronic density, its gradient, the exchange-correlation energy densities, and their derivatives. These units can be accessed in the standard way:julia> using LibXC, UnitfulHartree, Unitful\n\njulia> 1u\"ρ\" === 1u\"rho\" == 1u\"a₀^-3\"\ntrue\n\njulia> 1u\"∇ρ\" === 1u\"grho\" == 1u\"a₀^-4\"\ntrue\n\njulia> 1u\"ϵ\" === 1u\"Exc\" === 1u\"Eₕ\" ≈ 27.211386034310873u\"eV\"\ntrue\n\njulia> 1u\"∂ϵ_∂ρ\" == 1u\"Eₕ*a₀^3\"\ntrue\n\njulia> 1u\"∂²ϵ_∂ρ∂∇ρ\" == 1u\"Eₕ*a₀^7\"\ntrue\n\njulia> 1u\"∂³ϵ_∂ρ²∂∇ρ\" == 1u\"Eₕ*a₀^10\"\ntrueρ, ∇ρ (gradient of ρ) and ϵ have non-unicode aliases, for ease of access. The energy derivatives do not."
},

{
    "location": "index.html#Using-the-functionals-1",
    "page": "Home",
    "title": "Using the functionals",
    "category": "section",
    "text": "Once a functional is created, it can be called with a number of methods to compute the energy, the potential, as well as the second and third energy derivatives (when available for that functional).julia> func = XCFunctional(:lda_x, false);\n\njulia> energy(func, Cdouble[1, 2, 3]u\"rho\")\n3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:\n -0.738559 Eₕ\n -0.930526 Eₕ\n  -1.06519 EₕNote that we create an array of Cdouble (with the right units, as well). The underlying C library expects this type. Other types (and units, if not in Hartree atomic units) will incur the cost of creating of a new array with the right type.The following functions are available:energy\npotential\nenergy_and_potential\nsecond_energy_derivative\nthird_energy_derivative\nlda (all possible lda for the given functional outputs)\ngga (all possible gga outputs for the given functional)All these functions have overloads which hide the creation of a functional from the user:julia> energy(:lda_x, [1, 2, 3]u\"ρ\")\n3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:\n -0.738559 Eₕ\n -0.930526 Eₕ\n  -1.06519 Eₕ\n\njulia> energy(:lda_x, [1 2 3; 3 2 1]u\"ρ\")\n3-element Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},1}:\n -1.23917 Eₕ\n -1.17239 Eₕ\n -1.23917 Eₕ\n\njulia> energy(:lda_x, false, [1 2 3; 3 2 1]u\"ρ\")\n2×3 Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},2}:\n -0.738559 Eₕ  -0.930526 Eₕ   -1.06519 Eₕ\n  -1.06519 Eₕ  -0.930526 Eₕ  -0.738559 EₕIn most cases, the overhead of creating and destroying a C functional object at each call is likely too small to matter.The spin-polarization can be specified in the second argument (true for spin-polarized, false for spin-polarized). If this argument is not given, then a best-guess attempt is made: the functional is spin-polarized when ρ is at least two-dimensional and the first dimension of ρ is two (size(ρ, 1) == 2), and the functional is unpolarized in all other cases.Finally, it is possible to give inputs in different units. However, this will incur the cost of converting the array to the Hartree atomic units, both in terms of memory (an extra array is allocated) and in terms of compute (the actual conversion). The return is always in atomic units:julia> energy(:lda_x, false, [1 2 3; 3 2 1]u\"nm^-3\")\n2×3 Array{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},2}:\n -0.0390828 Eₕ  -0.0492413 Eₕ  -0.0563672 Eₕ\n -0.0563672 Eₕ  -0.0492413 Eₕ  -0.0390828 Eₕ"
},

{
    "location": "index.html#Using-pre-allocated-output-array-1",
    "page": "Home",
    "title": "Using pre-allocated output array",
    "category": "section",
    "text": "Similar functions exist that take pre-allocated output arrays. Following Julia conventions, these functions are named energy!, potential!, etc... Each function named above has an xxx! counterpart.DocTestSetup = quote\n    using LibXC\n    using Unitful\n    func = XCFunctional(:lda_x, false)\nendjulia> ρ = Cdouble[1, 2, 3]u\"rho\";\n\njulia> ϵ = similar(ρ, LibXC.Units.ϵ{Cdouble});\n\njulia> ∂ϵ_∂ρ = similar(ρ, LibXC.Units.∂ϵ_∂ρ{Cdouble});\n\njulia> result = energy_and_potential!(func, ρ, ϵ, ∂ϵ_∂ρ)\n(ϵ = [-0.738559,-0.930526,-1.06519]u\"Eₕ\", ∂ϵ_∂ρ = [-0.984745,-1.2407,-1.42025]u\"∂ϵ_∂ρ\")\n\njulia> result.ϵ === ϵ\ntrueFor convenience, some of the functions with more complex outputs return a named tuple. However, notice that the arrays in the tuple are aliases to the input arrays."
},

{
    "location": "index.html#LibXC.XCFunctional",
    "page": "Home",
    "title": "LibXC.XCFunctional",
    "category": "Type",
    "text": "Creates a functional from it's name and polarization\n\nThe name should be a one of the following symbols:\n\nLDA: lda_c_rpa, lda_c_vwn_rpa, lda_c_vwn_3, lda_x_2d, lda_x_1d, lda_c_pz, lda_c_gl, lda_c_gombas, lda_c_wigner, lda_c_1d_loos, lda_xc_ksdt, lda_c_xalpha, lda_c_vwn_4, lda_c_1d_csc, lda_x, lda_k_lp, lda_c_rc04, lda_c_ml2, lda_c_ob_pz, lda_c_pw_mod, lda_xc_teter93, lda_c_vwn_1, lda_c_pw_rpa, lda_c_vwn_2, lda_c_pz_mod, lda_c_2d_amgb, lda_c_pw, lda_xc_zlp, lda_c_vbh, lda_c_ob_pw, lda_c_vwn, lda_c_ml1, lda_k_tf, lda_c_hl, lda_c_2d_prm\nGGA: gga_xc_vv10, gga_x_lambda_ch_n, gga_x_lv_rpw86, gga_x_b86, gga_c_wl, gga_xc_hcth_p14, gga_x_mb88, gga_x_sogga, gga_k_pearson, gga_x_herman, gga_x_htbs, gga_c_am05, gga_x_c09x, gga_c_regtpss, gga_x_dk87_r2, gga_k_dk, gga_c_n12_sx, gga_c_wi, gga_x_n12, gga_k_fr_pw86, gga_x_hjs_b97x, gga_x_am05, gga_k_vjks, gga_c_pbe, gga_k_gp85, gga_c_pw91, gga_x_ev93, gga_k_perdew, gga_xc_th3, gga_xc_opbe_d, gga_k_llp, gga_xc_b97_d, gga_c_pbeloc, gga_x_pbe, gga_x_pbe_mol, gga_xc_hcth_407, gga_x_sogga11, gga_k_yt65, gga_k_gr, gga_c_pbefe, gga_x_ol2, gga_x_ssb, gga_x_2d_b86, gga_x_hjs_b88, gga_c_lm, gga_x_lbm, gga_k_tw2, gga_x_gam, gga_c_apbe, gga_k_ernzerhof, gga_x_b88, gga_x_2d_pbe, gga_c_pbe_sol, gga_k_fr_b88, gga_x_rpw86, gga_x_pbeint, gga_k_tfvw, gga_x_mpw91, gga_xc_mohlyp, gga_xc_mohlyp2, gga_xc_pbe1w, gga_c_op_pbe, gga_c_op_g96, gga_x_pbe_jsjr, gga_x_rge2, gga_c_pbe_jrgx, gga_k_lieb, gga_c_revtca, gga_x_ityh, gga_c_zpbesol, gga_k_tw1, gga_x_pbe_sol, gga_x_pw86, gga_k_vsk, gga_c_bgcp, gga_x_hjs_b88_v2, gga_x_dk87_r1, gga_x_hjs_pbe_sol, gga_x_ft97_b, gga_xc_hcth_407p, gga_x_lg93, gga_x_b86_r, gga_c_op_xalpha, gga_xc_kt2, gga_k_golden, gga_x_vmt_pbe, gga_x_wc, gga_c_rge2, gga_xc_th2, gga_x_optb88_vdw, gga_x_vmt84_pbe, gga_x_hjs_pbe, gga_x_g96, gga_x_rpbe, gga_x_2d_b88, gga_xc_edf1, gga_xc_th_fc, gga_k_ol2, gga_x_ssb_d, gga_c_xpbe, gga_k_apbe, gga_xc_th_fco, gga_c_op_pw91, gga_xc_b97_gga1, gga_x_vmt_ge, gga_c_sogga11, gga_x_2d_b86_mgc, gga_x_optpbe_vdw, gga_x_pbe_tca, gga_c_tca, gga_x_vmt84_ge, gga_x_pbea, gga_xc_hcth_p76, gga_k_vw, gga_x_ft97_a, gga_x_bayesian, gga_k_ludena, gga_xc_oblyp_d, gga_x_bpccac, gga_xc_opwlyp_d, gga_c_spbe, gga_x_mpbe, gga_x_pbek1_vdw, gga_c_sogga11_x, gga_c_optc, gga_xc_hcth_120, gga_c_zpbeint, gga_k_tw4, gga_c_lyp, gga_x_q2d, gga_x_wpbeh, gga_xc_th_fl, gga_xc_hcth_93, gga_xc_mpwlyp1w, gga_c_ft97, gga_x_apbe, gga_c_wi0, gga_x_kt1, gga_x_pw91, gga_x_airy, gga_k_lc94, gga_k_revapbe, gga_xc_xlyp, gga_k_revapbeint, gga_x_bgcp, gga_xc_pbelyp1w, gga_k_thakkar, gga_x_b86_mgc, gga_x_ssb_sw, gga_x_optx, gga_k_meyer, gga_k_baltin, gga_c_q2d, gga_x_pbefe, gga_c_p86, gga_c_op_b88, gga_k_apbeint, gga_x_pbe_r, gga_x_hcth_a, gga_c_pbeint, gga_xc_th4, gga_x_lambda_oc2_n, gga_c_n12, gga_x_cap, gga_xc_th1, gga_k_absp2, gga_k_ol1, gga_xc_hcth_147, gga_x_sfat, gga_xc_th_fcfo, gga_k_absp1, gga_x_xpbe, gga_c_gam, gga_c_hcth_a, gga_x_lag, gga_x_ak13, gga_k_tw3, gga_x_lambda_lo_n, gga_x_lb, gga_k_ge2\n\nThe second argument is true if the functional should be polarized.\n\n\n\n"
},

{
    "location": "index.html#LibXC.XCFunctional",
    "page": "Home",
    "title": "LibXC.XCFunctional",
    "category": "Type",
    "text": "Manages pointer to C libxc funtional data \n\n\n\n"
},

{
    "location": "index.html#LibXC.citations-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.citations",
    "category": "Method",
    "text": "List of journal references \n\n\n\n"
},

{
    "location": "index.html#LibXC.description-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.description",
    "category": "Method",
    "text": "Long name/description of a functional \n\n\n\n"
},

{
    "location": "index.html#LibXC.energy!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},N}}",
    "page": "Home",
    "title": "LibXC.energy!",
    "category": "Method",
    "text": "energy!(func, ρ, ∇ρ, ϵ)\n\n\nGGA energy as a function of ρ and ∇ρ=|∇ρ|. The dimensionality is as follows:\n\nGGA unpolarized polarized\nρ any (2, ...)\n∇ρ size(ρ) (3, size(ρ)[2:end]...)\nϵ size(ρ) size(ρ)[2:end]\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},N}}",
    "page": "Home",
    "title": "LibXC.energy!",
    "category": "Method",
    "text": "energy!(func, ρ, ϵ)\n\n\nComputes the energy in-place for a given LDA functional. For spin-unpolarized functionals, the output array has the dimensions of ρ. For spin-polarized functionals, assuming ndims(ρ) > 1 && size(ρ, 1) == 2, it is size(ρ)[2:end].\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.energy!",
    "category": "Method",
    "text": "energy!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.energy!",
    "category": "Method",
    "text": "energy!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N}}",
    "page": "Home",
    "title": "LibXC.energy",
    "category": "Method",
    "text": "energy(func, ρ)\n\n\nComputes the LDA energy density associated with the input density.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.energy",
    "category": "Method",
    "text": "energy(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.energy",
    "category": "Method",
    "text": "energy(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential!-Tuple{LibXC.AbstractLibXCFunctional,DenseArray{Float64,N},DenseArray{Float64,N},DenseArray{Float64,N}}",
    "page": "Home",
    "title": "LibXC.energy_and_potential!",
    "category": "Method",
    "text": "energy_and_potential!(func, ρ, ϵ, ∂ϵ_∂ρ)\n\n\nLDA energy and first derivative\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^2 𝐌 𝐓^-2}, Units:{Eₕ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^5 𝐌 𝐓^-2}, Units:{∂ϵ_∂ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^6 𝐌 𝐓^-2}, Units:{∂ϵ_∂∇ρ}},N}}",
    "page": "Home",
    "title": "LibXC.energy_and_potential!",
    "category": "Method",
    "text": "energy_and_potential!(func, ρ, ∇ρ, ϵ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)\n\n\nGGA energy and potential\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.energy_and_potential!",
    "category": "Method",
    "text": "energy_and_potential!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.energy_and_potential!",
    "category": "Method",
    "text": "energy_and_potential!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.energy_and_potential",
    "category": "Method",
    "text": "energy_and_potential(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.energy_and_potential-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.energy_and_potential",
    "category": "Method",
    "text": "energy_and_potential(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.family-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.family",
    "category": "Method",
    "text": "Whether LDA, GGA, etc \n\n\n\n"
},

{
    "location": "index.html#LibXC.flags-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.flags",
    "category": "Method",
    "text": "Set of flags describing the functional \n\n\n\n"
},

{
    "location": "index.html#LibXC.gga!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Float64,N},DenseArray{Float64,N},DenseArray{Float64,N},Vararg{T<:DenseArray{Float64,N},N}}",
    "page": "Home",
    "title": "LibXC.gga!",
    "category": "Method",
    "text": "gga!(func, ρ, ∇ρ, ϵ, outputs)\n\n\nGGA energy, first, second, and third derivatives, in place. Arrays for the first, second, and third derivatives are optional. They should be given only if available for that particular functional. When requesting higher derivatives, arrays to store the lower derivatives should also be given.\n\n\n\n"
},

{
    "location": "index.html#LibXC.gga-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N}}",
    "page": "Home",
    "title": "LibXC.gga",
    "category": "Method",
    "text": "gga(func, ρ, ∇ρ)\n\n\nComputes the energy and all available derivatives for the given functional\n\n\n\n"
},

{
    "location": "index.html#LibXC.kind-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.kind",
    "category": "Method",
    "text": "Whether Exchange, Correlation, or Exchange-Correlation \n\n\n\n"
},

{
    "location": "index.html#LibXC.lda!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Float64,N},DenseArray{Float64,N},Vararg{T<:DenseArray{Float64,N},N}}",
    "page": "Home",
    "title": "LibXC.lda!",
    "category": "Method",
    "text": "lda!(func, ρ, ϵ, outputs)\n\n\nLDA energy, first, second, and third derivatives. The first, second, and third derivatives are optional. It is an error to request a derivative of the functional that is not implemented in the underlying C library.\n\n\n\n"
},

{
    "location": "index.html#LibXC.lda!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.lda!",
    "category": "Method",
    "text": "lda!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.lda!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.lda!",
    "category": "Method",
    "text": "lda!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.lda-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Float64,N}}",
    "page": "Home",
    "title": "LibXC.lda",
    "category": "Method",
    "text": "lda(func, ρ)\n\n\nComputes the energy and all available derivatives for the given functional\n\n\n\n"
},

{
    "location": "index.html#LibXC.lda-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.lda",
    "category": "Method",
    "text": "lda(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.lda-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.lda",
    "category": "Method",
    "text": "lda(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.libxc_functionals-Tuple{}",
    "page": "Home",
    "title": "LibXC.libxc_functionals",
    "category": "Method",
    "text": "Names of the available libxc functionals \n\n\n\n"
},

{
    "location": "index.html#LibXC.potential!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^5 𝐌 𝐓^-2}, Units:{∂ϵ_∂ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^6 𝐌 𝐓^-2}, Units:{∂ϵ_∂∇ρ}},N}}",
    "page": "Home",
    "title": "LibXC.potential!",
    "category": "Method",
    "text": "potential!(func, ρ, ∇ρ, ∂ϵ_∂ρ, ∂ϵ_∂∇ρ)\n\n\nGGA potential computed in place. The dimensionality of the different arrays are as follows:\n\nGGA unpolarized polarized\nρ any (2, ...)\n∇ρ size(ρ) (3, size(ρ)[2:end]...)\n∂ϵ/∂ρ size(ρ) size(ρ)\n∂ϵ/∂∇ρ size(ρ) (3, size(ρ)[2:end]...)\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^5 𝐌 𝐓^-2}, Units:{∂ϵ_∂ρ}},N}}",
    "page": "Home",
    "title": "LibXC.potential!",
    "category": "Method",
    "text": "potential!(func, ρ, ∂ϵ_∂ρ)\n\n\nComputes the potential in-place for a given LDA functional. For spin-unpolarized functionals, the output array has the dimensions of ρ. For spin-polarized functionals, assuming ndims(ρ) > 1 && size(ρ, 1) == 2, it is size(ρ).\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.potential!",
    "category": "Method",
    "text": "potential!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.potential!",
    "category": "Method",
    "text": "potential!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N}}",
    "page": "Home",
    "title": "LibXC.potential",
    "category": "Method",
    "text": "potential(func, ρ)\n\n\nComputes the LDA potential associated with the input density.\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.potential",
    "category": "Method",
    "text": "potential(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.potential-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.potential",
    "category": "Method",
    "text": "potential(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N}}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative",
    "category": "Method",
    "text": "second_energy_derivative(func, ρ)\n\n\nComputes the LDA second energy derivatives associated with the input density.\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative",
    "category": "Method",
    "text": "second_energy_derivative(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative",
    "category": "Method",
    "text": "second_energy_derivative(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.spin-Tuple{LibXC.CFuncType}",
    "page": "Home",
    "title": "LibXC.spin",
    "category": "Method",
    "text": "Spin polarization \n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N}}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative",
    "category": "Method",
    "text": "third_energy_derivative(func, ρ)\n\n\nComputes the LDA third energy derivatives associated with the input density.\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative-Tuple{Symbol,DenseArray}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative",
    "category": "Method",
    "text": "third_energy_derivative(name, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative",
    "category": "Method",
    "text": "third_energy_derivative(name, spin, ρ)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly specified.\n\n\n\n"
},

{
    "location": "index.html#LibXC.CReal",
    "page": "Home",
    "title": "LibXC.CReal",
    "category": "Constant",
    "text": "Floating points LibXC might now about \n\n\n\n"
},

{
    "location": "index.html#LibXC.FUNCTIONALS",
    "page": "Home",
    "title": "LibXC.FUNCTIONALS",
    "category": "Constant",
    "text": "Functional names and keys \n\n\n\n"
},

{
    "location": "index.html#LibXC.iFUNCTIONALS",
    "page": "Home",
    "title": "LibXC.iFUNCTIONALS",
    "category": "Constant",
    "text": "\" Functional keys and names \n\n\n\n"
},

{
    "location": "index.html#LibXC.AllGGA",
    "page": "Home",
    "title": "LibXC.AllGGA",
    "category": "Type",
    "text": "All outputs from LDA \n\n\n\n"
},

{
    "location": "index.html#LibXC.AllLDA",
    "page": "Home",
    "title": "LibXC.AllLDA",
    "category": "Type",
    "text": "All outputs from LDA \n\n\n\n"
},

{
    "location": "index.html#LibXC.Citation",
    "page": "Home",
    "title": "LibXC.Citation",
    "category": "Type",
    "text": "Holds citation data \n\n\n\n"
},

{
    "location": "index.html#LibXC.GGAEnergyAndPotential",
    "page": "Home",
    "title": "LibXC.GGAEnergyAndPotential",
    "category": "Type",
    "text": "Holds GGA energy and first derivatives \n\n\n\n"
},

{
    "location": "index.html#LibXC.GGAPotential",
    "page": "Home",
    "title": "LibXC.GGAPotential",
    "category": "Type",
    "text": "Energy and potential from LDA \n\n\n\n"
},

{
    "location": "index.html#LibXC.GGASecondDerivative",
    "page": "Home",
    "title": "LibXC.GGASecondDerivative",
    "category": "Type",
    "text": "Second derivative from GGA\n\nInclude the second derivative of the energy with respect to ρ, ∇ρ, and both ρ and ∇ρ.\n\n\n\n"
},

{
    "location": "index.html#LibXC.GGAThirdDerivative",
    "page": "Home",
    "title": "LibXC.GGAThirdDerivative",
    "category": "Type",
    "text": "Third derivative from GGA\n\nInclude the third derivative of the energy with respect to ρ, ∇ρ, and both ρ and ∇ρ.\n\n\n\n"
},

{
    "location": "index.html#LibXC.LDAEnergyAndPotential",
    "page": "Home",
    "title": "LibXC.LDAEnergyAndPotential",
    "category": "Type",
    "text": "Energy and potential from LDA \n\n\n\n"
},

{
    "location": "index.html#LibXC.libkey-Tuple{LibXC.CFuncInfoType}",
    "page": "Home",
    "title": "LibXC.libkey",
    "category": "Method",
    "text": "Integer key of the functional \n\n\n\n"
},

{
    "location": "index.html#LibXC.output_size-Tuple{Bool,Tuple{Vararg{T,N}},Integer}",
    "page": "Home",
    "title": "LibXC.output_size",
    "category": "Method",
    "text": "output_size(polarized, dims, factor)\n\n\nHelps determine size of an output. dims refers to the size of input density array, and factor is the size we expect for the first dimension. It will vary with the functional (LDA, GGA, ...) and the kind of output (energy, potential, ...).\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^8 𝐌 𝐓^-2}, Units:{∂²ϵ_∂ρ²}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^9 𝐌 𝐓^-2}, Units:{∂²ϵ_∂ρ∂∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^10 𝐌 𝐓^-2}, Units:{∂²ϵ_∂∇ρ²}},N}}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative!",
    "category": "Method",
    "text": "second_energy_derivative!(func, ρ, ∇ρ, ∂²ϵ_∂ρ², ∂²ϵ_∂ρ∂∇ρ, ∂²ϵ_∂∇ρ²)\n\n\nSecond derivatives of GGA energy w.r.t. ρ and ∇ρ=|∇ρ|. The dimensionality of the arrays is as follows:\n\nGGA unpolarized polarized\nρ any (2, ...)\n∇ρ size(ρ) (3, size(ρ)[2:end]...)\n∂²ϵ/∂ρ² size(ρ) (3, size(ρ)[2:end]...)\n∂²ϵ/∂ρ∂∇ρ size(ρ) (6, size(ρ)[2:end]...)\n∂²ϵ/∂∇ρ² size(ρ) (6, size(ρ)[2:end]...)\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^8 𝐌 𝐓^-2}, Units:{∂²ϵ_∂ρ²}},N}}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative!",
    "category": "Method",
    "text": "second_energy_derivative!(func, ρ, ∂²ϵ_∂ρ²)\n\n\nComputes the second energy derivative in-place for a given LDA functional. For spin-unpolarized functionals, the output array has the dimensions of ρ. For spin-polarized functionals, assuming ndims(ρ) > 1 && size(ρ, 1) == 2, it is (3, size(ρ)[2:end]...).\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative!",
    "category": "Method",
    "text": "second_energy_derivative!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.second_energy_derivative!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.second_energy_derivative!",
    "category": "Method",
    "text": "second_energy_derivative!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^-4}, Units:{∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^11 𝐌 𝐓^-2}, Units:{∂³ϵ_∂ρ³}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^12 𝐌 𝐓^-2}, Units:{∂³ϵ_∂ρ²∂∇ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^13 𝐌 𝐓^-2}, Units:{∂³ϵ_∂ρ∂∇ρ²}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^14 𝐌 𝐓^-2}, Units:{∂³ϵ_∂∇ρ³}},N}}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative!",
    "category": "Method",
    "text": "third_energy_derivative!(func, ρ, ∇ρ, ∂³ϵ_∂ρ³, ∂³ϵ_∂ρ²∂∇ρ, ∂³ϵ_∂ρ∂∇ρ², ∂³ϵ_∂∇ρ³)\n\n\nThird derivatives of GGA energy w.r.t. ρ and ∇ρ=|∇ρ|. The dimensionality of the arrays is as follows:\n\nGGA unpolarized polarized\nρ any (2, ...)\n∇ρ size(ρ) (3, size(ρ)[2:end]...)\n∂³ϵ/∂ρ³ size(ρ) (4, size(ρ)[2:end]...)\n∂³ϵ/∂ρ²∂∇ρ size(ρ) (9, size(ρ)[2:end]...)\n∂³ϵ/∂ρ∂∇ρ² size(ρ) (10, size(ρ)[2:end]...)\n∂³ϵ/∂∇ρ³ size(ρ) (12, size(ρ)[2:end]...)\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative!-Tuple{LibXC.AbstractLibXCFunctional{Float64},DenseArray{Quantity{Float64, Dimensions:{𝐋^-3}, Units:{ρ}},N},DenseArray{Quantity{Float64, Dimensions:{𝐋^11 𝐌 𝐓^-2}, Units:{∂³ϵ_∂ρ³}},N}}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative!",
    "category": "Method",
    "text": "third_energy_derivative!(func, ρ, ∂³ϵ_∂ρ³)\n\n\nComputes the third energy derivative in-place for a given LDA functional. For spin-unpolarized functionals, the output array has the dimensions of ρ. For spin-polarized functionals, assuming ndims(ρ) > 1 && size(ρ, 1) == 2, it is (4, size(ρ)[2:end]...).\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative!-Tuple{Symbol,DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative!",
    "category": "Method",
    "text": "third_energy_derivative!(name, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is determined from the dimensionality of ρ: ndims(ρ) > 1 && size(ρ, 1) == 2.\n\n\n\n"
},

{
    "location": "index.html#LibXC.third_energy_derivative!-Tuple{Symbol,Union{Bool,LibXC.Constants.SPIN},DenseArray,Vararg{Any,N}}",
    "page": "Home",
    "title": "LibXC.third_energy_derivative!",
    "category": "Method",
    "text": "third_energy_derivative!(name, spin, ρ, args)\n\n\nA simple way to call a functional using it's name. Spin polarization is explicitly requested.\n\n\n\n"
},

{
    "location": "index.html#LibXC.to_markdown-Tuple{Any}",
    "page": "Home",
    "title": "LibXC.to_markdown",
    "category": "Method",
    "text": "Prints functional to markdown, mostly for docs \n\n\n\n"
},

{
    "location": "index.html#LibXC.@check_availability-Tuple{Any,Any}",
    "page": "Home",
    "title": "LibXC.@check_availability",
    "category": "Macro",
    "text": "Adds argument check for energy derivative availability \n\n\n\n"
},

{
    "location": "index.html#LibXC.@check_functional-Tuple{Any,Any}",
    "page": "Home",
    "title": "LibXC.@check_functional",
    "category": "Macro",
    "text": "Adds check for functional type \n\n\n\n"
},

{
    "location": "index.html#LibXC.@check_size-Tuple{Any,Any,Any,Any}",
    "page": "Home",
    "title": "LibXC.@check_size",
    "category": "Macro",
    "text": "Adds argument check for size compatibility \n\n\n\n"
},

{
    "location": "index.html#API-1",
    "page": "Home",
    "title": "API",
    "category": "section",
    "text": "Modules = [LibXC]"
},

{
    "location": "index.html#Available-LDA-functionals-1",
    "page": "Home",
    "title": "Available LDA functionals",
    "category": "section",
    "text": "using Base.Markdown\nusing LibXC\nresult = Base.Markdown.MD()\nfor (name, key) in sort([u for u in LibXC.FUNCTIONALS], by=x -> x[2])\n    if family(name) == LibXC.Constants.lda\n        push!(result.content, LibXC.to_markdown(name))\n    end\nend\nresult"
},

{
    "location": "index.html#Available-GGA-functionals-1",
    "page": "Home",
    "title": "Available GGA functionals",
    "category": "section",
    "text": "using Base.Markdown\nusing LibXC\nresult = Base.Markdown.MD()\nfor (name, key) in sort([u for u in LibXC.FUNCTIONALS], by=x -> x[2])\n	if family(name) == LibXC.Constants.gga\n		push!(result.content, LibXC.to_markdown(name))\n    end\nend\nresult"
},

]}
