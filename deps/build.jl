using BinDeps

@BinDeps.setup

libxc = library_dependency("libxc")
provides(Sources,
         URI("https://github.com/mdavezac/libxc/archive/v3.0.0.tar.gz"),
         libxc, unpacked_dir="libxc-3.0.0")
options = ["--enable-static=no", "--enable-shared=yes", "--disable-fortran"]
autotools = Autotools(libtarget="src/libxc.la", configure_options=options)
provides(BuildProcess, autotools, libxc)
@BinDeps.install Dict(:libxc => :libxc)

if "UnitfulHartree" ∉ keys(Pkg.installed())
    Pkg.clone("https://github.com/mdavezac/UnitfulHartree.jl.git")
end
