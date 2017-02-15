module LibXC

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("LiBXC not properly installed. Please run Pkg.build(\"LibXC\")")
end

include("structures.jl")
include("constants.jl")

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

end # module
