# LibXC

[![Build Status](https://travis-ci.org/mdavezac/LibXC.jl.svg?branch=master)](https://travis-ci.org/mdavezac/LibXC.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mdavezac.github.io/LibXC.jl/latest)
[![Coverage Status](https://coveralls.io/repos/github/mdavezac/LibXC.jl/badge.svg?branch=master)](https://coveralls.io/github/mdavezac/LibXC.jl?branch=master)

Wrappers around the library of exchange-correlation functional
[libxc](octopus-code.org/wiki/Libxc)

## Installation

```julia
Pkg.clone("https://github.com/mdavezac/LibXC.jl")
Pkg.build("LibXC")
```

## References

All credit goes to the original C implementation:

* Miguel A. L. Marques, Micael J. T. Oliveira, and Tobias Burnus,
  "Libxc: a library of exchange and correlation functionals for density functional theory"
  [Comput. Phys. Commun. 183, 2272 (2012)](http://dx.doi.org/10.1016/j.cpc.2012.05.007)
  OAI: [arXiv:1203.1739](http://arxiv.org/abs/1203.1739)
