module LibXC
export description, kind, family, flags, citations, spin, energy, energy!
export potential, potential!

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

include("structures.jl")
include("constants.jl")

""" Holds citation data """
immutable Citation
    ref::String
    doi::String
    bibtex::String
end
Base.show(io::IO, x::Citation) = show(io, x.ref)



const LIB_VERSION = let
    version = Cint[0, 0, 0]
    @eval ccall((:xc_version, $libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}),
                Ref($version), Ref($version, 2), Ref($version, 3))
    VersionNumber(version...)
end



""" Functional names and keys """
const FUNCTIONALS = let
    keys = @eval cglobal(
        (:xc_functional_keys, $(LibXC.libxc)), LibXC.CFunctionalKey)
    result = Dict{Symbol, Int32}()

    i, key = 1, unsafe_load(keys)
    while key.key > 0
        n = findfirst(x -> x == UInt32('\0'), key.name) - 1
        push!(result, Symbol(join(Char(key.name[u]) for u in 1:n)) => key.key)
        i += 1
        key = unsafe_load(keys, i)
    end
    result
end

abstract AbstractXCFunctional
abstract AbstractLibXCFunctional{T <: CReal} <: AbstractXCFunctional

type XCFunctional{T <: CReal} <: AbstractLibXCFunctional{T}
    c_ptr::Ptr{CFuncType{T}}
    XCFunctional(ptr::Ptr{CFuncType{T}}) = new(ptr)
end

function XCFunctional(name::Symbol, spin_polarized::Bool=true)
    name ∉ keys(FUNCTIONALS) && error("Functional $name does not exist")
    ptr = ccall((:xc_func_alloc, libxc), Ptr{CFuncType{Cdouble}}, ())
    functional = XCFunctional{Cdouble}(ptr)
    err = ccall((:xc_func_init, libxc), Cint, (Ptr{CFuncType{Cdouble}}, Cint, Cint),
                ptr, FUNCTIONALS[name], spin_polarized ? 2: 1)
    err ≠ 0 && error("Error $err encountered in LibXC")
    finalizer(functional, _delete_libxc_functional)
    functional
end
function XCFunctional(name::Symbol, spin::Constants.SPIN)
    XCFunctional(name, spin == Constants.polarized)
end
XCFunctional(n::Symbol, spin::Integer) = XCFunctional(n, convert(Constants.SPIN, spin))

function _func_info(func::AbstractLibXCFunctional{Cdouble})
    func_type = ccall((:xc_func_get_info, libxc), Ptr{CFuncInfoType{Cdouble}},
                      (Ptr{CFuncType{Cdouble}}, ), func.c_ptr)
    unsafe_load(func_type)
end

function _delete_libxc_functional(func::AbstractLibXCFunctional{Cdouble})
    if func.c_ptr ≠ C_NULL
        ccall((:xc_func_end, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
        ccall((:xc_func_free, libxc), Void, (Ptr{CFuncType},), func.c_ptr)
    end
end

description(info::CFuncInfoType) = unsafe_string(info.name)
""" Integer key of the functional """
libkey(info::CFuncInfoType) = info.number
""" Whether Exchange, Correlation, or Exchange-Correlation """
kind(info::CFuncInfoType) = convert(Constants.KIND, info.kind)
""" Whether LDA, GGA, etc """
family(info::CFuncInfoType) = convert(Constants.FAMILY, info.family)
""" Spin polarization """
spin(func::CFuncType) = convert(Constants.SPIN, func.nspin)

""" Set of flags describing the functional """
function flags(func::CFuncInfoType)
    result = Set{Constants.FLAGS}()
    for flag in Base.instances(Constants.FLAGS)
        (convert(Cint, flag) & func.flags) ≠ 0 && push!(result, flag)
    end
    result
end

""" List of journal references """
function citations(info::CFuncInfoType)
    result = Citation[]
    for ptr in info.refs
        ptr == C_NULL && break
        ref = unsafe_load(ptr)
        push!(result, Citation(unsafe_string(ref.ref), unsafe_string(ref.doi),
                               unsafe_string(ref.bibtex)))
    end
    result
end

description(func::AbstractLibXCFunctional) = description(_func_info(func))
libkey(func::AbstractLibXCFunctional) = libkey(_func_info(func))
kind(func::AbstractLibXCFunctional) = kind(_func_info(func))
family(func::AbstractLibXCFunctional) = family(_func_info(func))
flags(func::AbstractLibXCFunctional) = flags(_func_info(func))
citations(func::AbstractLibXCFunctional) = citations(_func_info(func))
spin(func::AbstractLibXCFunctional) = convert(Constants.SPIN, unsafe_load(func.c_ptr).nspin)

""" Determines the size of the energy from the size of ρ and spin polarization """
function εxc_size(polarized::Bool, dims::NTuple)
    if length(dims) == 0
        throw(ArgumentError("Empty size tuple"))
    elseif !polarized
        dims
    elseif length(dims) == 1 && polarized
        if dims[1] % 2 ≠ 0
            throw(ArgumentError("Odd array size for polarized functional"))
        end
        warn("Spin polarized function, but dimensionality of ρ is 1")
        (dims[1]/2,)
    elseif dims[1] == 2
        dims[2:end]
    else
        throw(ArgumentError("Spin polarization expects size(ρ, 2) == 2"))
    end
end
εxc_size(func::AbstractLibXCFunctional, ρ::DenseArray) = εxc_size(spin(func), size(ρ))
εxc_size(s::Constants.SPIN, dims::NTuple) = εxc_size(s == Constants.polarized, dims)
εxc_size(s::Union{Bool, Constants.SPIN}, ρ::DenseArray) = εxc_size(s, size(ρ))


for (funcname, name, sizer) ∈ [(:xc_lda_exc, :energy, :(εxc_size(func, ρ))),
                               (:xc_lda_vxc, :potential, :(size(ρ)))]
    local name! = Symbol("$(name)!")
    @eval begin
        function $name!(func::AbstractLibXCFunctional{Cdouble}, ρ::DenseArray{Cdouble},
                        $name::DenseArray{Cdouble})
            if family(func) ≠ Constants.lda
                msg = "Incorrect number of arguments: input is not an LDA functional"
                throw(ArgumentError(msg))
            end
            if size($name) ≠ $sizer
                throw(ArgumentError("sizes of ρ and input are incompatible"))
            end

            ccall(($(parse(":$funcname")), libxc), Void,
                  (Ptr{CFuncType}, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
                  func.c_ptr, length(ρ), ρ, $name)
            $name
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
            $name!(func, ρ, similar(ρ, eltype(ρ), $sizer))
        end
    end
end

end # module
