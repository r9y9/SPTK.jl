# SPTK.jl

[![Build Status](https://travis-ci.org/r9y9/SPTK.jl.svg?branch=master)](https://travis-ci.org/r9y9/SPTK.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/a1byruqq7l19puu3/branch/master?svg=true)](https://ci.appveyor.com/project/r9y9/sptk-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/r9y9/SPTK.jl/badge.svg)](https://coveralls.io/r/r9y9/SPTK.jl)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)

SPTK.jl is a light-weight Julia wrapper for the [Speech Signal Processing Toolkit (SPTK)](http://sp-tk.sourceforge.net/). Note that SPTK.jl is based on [the modified version of SPTK](https://github.com/r9y9/SPTK).

## Installation

```julia
Pkg.clone("https://github.com/r9y9/SPTK.jl.git")
Pkg.build("SPTK)
```

## Exported functions

List of exported functions can be found in [src/SPTK.jl](src/SPTK.jl).
