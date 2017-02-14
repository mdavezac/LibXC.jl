module LibXC

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

include("structures.jl")
include("constants.jl")

function _get_version()
    version = Cint[0, 0, 0]
    ccall((:xc_version, libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}), Ref(version),
          Ref(version, 2), Ref(version, 3))
    VersionNumber(version...)
end

const lib_version = let version = Cint[0, 0, 0]
    @eval ccall((:xc_version, $libxc), Void, (Ref{Cint}, Ref{Cint}, Ref{Cint}),
                Ref($version), Ref($version, 2), Ref($version, 3))
    VersionNumber(version...)
end



""" Maximum number of functionals defined in the library """
const NFUNCTIONALS = let current = 1, start = 1, last = 10000
    # Performs bisection to find largest valid N
    # Brute force is a bit slow
    name = @eval ccall(
        (:xc_functional_get_name, $libxc), Cstring, (Cint,), $last)
    while name ≠ Cstring(C_NULL)
        last *= 2
        name = @eval ccall(
            (:xc_functional_get_name, $libxc), Cstring, (Cint,), $last)
    end

    while start + 1 < last
        current = div(start + last, 2)
        name = @eval ccall(
            (:xc_functional_get_name, $libxc), Cstring, (Cint,), $current)
        if name == Cstring(C_NULL)
            last = current
        else
            start = current
        end
    end
    start
end

_xc_key(name::AbstractString) = ccall((:xc_functional_get_number), Cint, (Cstring,), name)
function _xc_name(key::Integer)
    name = ccall((:xc_functional_get_name), Cstring, (Cint,), key)
    unsafe_string(name)
end

end # module
