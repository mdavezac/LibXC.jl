typealias CReal Union{Cfloat, Cdouble}

immutable CFuncReferenceType
    ref::Cstring
    doi::Cstring
    bibtex::Cstring
end

immutable CFuncInfoType{FLOAT <: CReal}
    """ identifier number """
    number::Cint
    """ XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC """
    kind::Cint

    """ name of the functional, e.g. "PBE" """
    name::Cstring

    """ type of the functional, e.g. XC_FAMILY_GGA """
    family::Cint

    """ index of the references """
    refs::Ptr{Ptr{CFuncReferenceType}}

    flags::Cint

    min_dens::FLOAT
    min_grad::FLOAT
    min_tau::FLOAT
    min_zeta::FLOAT

    """
    ```C
        void (*init)(struct XC(func_type) *p)
    ```
    """
    init::Ptr{Void}
    """
    ```C
        void (*destroy) (struct XC(func_type) *p)
    ```
    """
    destroy::Ptr{Void}
    """
    ```C
        void (*lda) (const struct XC(func_type) *p, int np, const FLOAT *rho, FLOAT *zk,
                     FLOAT *vrho, FLOAT *v2rho2, FLOAT *v3rho3)
    ```
    """
    lda::Ptr{Void}
    """
    ```C
        void (*gga) (const struct XC(func_type) *p, int np, 
                     const FLOAT *rho, const FLOAT *sigma, 
                     FLOAT *zk, FLOAT *vrho, FLOAT *vsigma,
                     FLOAT *v2rho2, FLOAT *v2rhosigma, FLOAT *v2sigma2,
                     FLOAT *v3rho3, FLOAT *v3rho2sigma, FLOAT *v3rhosigma2, FLOAT *v3sigma3)
    ```
    """
    gga::Ptr{Void}
    """
    ```C
        void (*mgga)(const struct XC(func_type) *p, int np, 
                     const FLOAT *rho, const FLOAT *sigma, const FLOAT *lapl_rho, const
                     FLOAT *tau, FLOAT *zk, FLOAT *vrho, FLOAT *vsigma, FLOAT *vlapl_rho,
                     FLOAT *vtau, FLOAT *v2rho2, FLOAT *v2sigma2, FLOAT *v2tau2, FLOAT
                     *v2lapl2, FLOAT *v2rhosigma, FLOAT *v2rhotau, FLOAT *v2rholapl, 
                     FLOAT *v2sigmatau, FLOAT *v2sigmalapl, FLOAT *v2taulapl)
    ```
    """
    mmga::Ptr{Void}
end

immutable CFuncType{FLOAT <: CReal}
    """ all the information concerning this functional """
    info::Ptr{CFuncInfoType{FLOAT}}
    """ XC_UNPOLARIZED or XC_POLARIZED """
    nspin::Cint

    """ number of auxiliary functions """
    n_func_aux::Cint
    """ most GGAs are based on a LDA or other GGAs """
    func_aux::Ptr{Ptr{CFuncType{FLOAT}}}

    """ coefficients for the mixing """
    mix_coef::FLOAT

    """ range-separation parameter for range-separated hybrids """
    cam_omega::FLOAT
    """ fraction of Hartree-Fock exchange for hybrids """
    cam_alpha::FLOAT
    """ fraction of short-range exchange for range-separated hybrids """
    cam_beta::FLOAT

    """ Non-local correlation, b parameter """
    nlc_b::FLOAT
    """ Non-local correlation, C parameter """
    nlc_C::FLOAT

    """ Shortcut in case of several functionals sharing the same interface """
    func::Cint
    """ spin dimensions of the arrays """
    n_rho::Cint
    """ spin dimensions of the arrays """
    n_sigma::Cint
    """ spin dimensions of the arrays """
    n_tau::Cint
    """ spin dimensions of the arrays """
    n_lapl::Cint
    n_zk::Cint

    n_vrho::Cint
    n_vsigma::Cint
    n_vtau::Cint
    n_vlapl::Cint

    n_v2rho2::Cint
    n_v2sigma2::Cint
    n_v2tau2::Cint
    n_v2lapl2::Cint

    n_v2rhosigma::Cint
    n_v2rhotau::Cint
    n_v2rholapl::Cint

    n_v2sigmatau::Cint
    n_v2sigmalapl::Cint
    n_v2lapltau::Cint


    n_v3rho3::Cint
    n_v3rho2sigma::Cint
    n_v3rhosigma2::Cint
    n_v3sigma3::Cint


    """ this allows us to fix parameters in the functional """
    params::Ptr{Void}
end

immutable CFunctionalKey
    """ Name of the functional """
    name::NTuple{256, Cchar}
    """ Key associated with the functional """
    key::Cint
end
