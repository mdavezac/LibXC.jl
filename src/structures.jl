immutable CFuncReferenceType
    ref::Cstring
    doi::Cstring
    bibtex::Cstring
end

immutable CFuncInfoType{FLOAT <: Union{Cfloat, Cdouble}}
    """ identifier number """
    number::CInt
    """ XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC """
    kind::CInt

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

immutable CFuncType
    """ all the information concerning this functional """
    info::Ptr{CFuncInfoType}
                     int nspin;                            /* XC_UNPOLARIZED or XC_POLARIZED  */

                     int n_func_aux;                       /* how many auxiliary functions we need */
                     struct XC(func_type) **func_aux;      /* most GGAs are based on a LDA or other GGAs  */
                     FLOAT *mix_coef;                      /* coefficients for the mixing */

                     FLOAT cam_omega;                      /* range-separation parameter for range-separated hybrids */
                     FLOAT cam_alpha;                      /* fraction of Hartree-Fock exchange for normal or range separated hybrids */
                     FLOAT cam_beta;                       /* fraction of short-range exchange for range-separated hybrids */

                     FLOAT nlc_b;                          /* Non-local correlation, b parameter */
                     FLOAT nlc_C;                          /* Non-local correlation, C parameter */

                     int func;                             /* Shortcut in case of several functionals sharing the same interface */
                     int n_rho, n_sigma, n_tau, n_lapl;    /* spin dimensions of the arrays */
                     int n_zk;

                     int n_vrho, n_vsigma, n_vtau, n_vlapl;

                     int n_v2rho2, n_v2sigma2, n_v2tau2, n_v2lapl2,
                     n_v2rhosigma, n_v2rhotau, n_v2rholapl, 
                     n_v2sigmatau, n_v2sigmalapl, n_v2lapltau;

                     int n_v3rho3, n_v3rho2sigma, n_v3rhosigma2, n_v3sigma3;

                     void *params;                         /* this allows us to fix parameters in the functional */
                    };
