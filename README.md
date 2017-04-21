# SPTK.jl

[![SPTK](http://pkg.julialang.org/badges/SPTK_0.4.svg)](http://pkg.julialang.org/?pkg=SPTK&ver=0.4)
[![SPTK](http://pkg.julialang.org/badges/SPTK_0.5.svg)](http://pkg.julialang.org/?pkg=SPTK&ver=0.5)
[![SPTK](http://pkg.julialang.org/badges/SPTK_0.6.svg)](http://pkg.julialang.org/?pkg=SPTK&ver=0.6)

[![Build Status](https://travis-ci.org/r9y9/SPTK.jl.svg?branch=master)](https://travis-ci.org/r9y9/SPTK.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/a1byruqq7l19puu3/branch/master?svg=true)](https://ci.appveyor.com/project/r9y9/sptk-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/r9y9/SPTK.jl/badge.svg)](https://coveralls.io/r/r9y9/SPTK.jl)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)

SPTK.jl is a Julia wrapper for the [Speech Signal Processing Toolkit (SPTK)](http://sp-tk.sourceforge.net/), which provides a lot of functionalities for speech signal processing such as linear prediction analysis, mel-cepstrum analysis, generalized cepstrum analysis and mel-generalized cepstrum analysis to name a few. See the original project page for more details.


**NOTE**: SPTK.jl is based on a modified version of SPTK ([r9y9/SPTK](https://github.com/r9y9/SPTK)).

## Documentation

A reference manual of the SPTK can be found at http://sp-tk.sourceforge.net/.

## Demonstration notebook

- [Introduction notebook](http://nbviewer.ipython.org/github/r9y9/SPTK.jl/blob/master/examples/Introduction%20to%20SPTK.jl.ipynb): a brief intruduction to SPTK.jl, especially focused on mel-generalized cepstrum analysis

## Getting started

Functions that SPTK.jl provides are basically same as the SPTK, so if you are new to SPTK, please take a look at the original documentation first and then use SPTK.jl for your need. Also notice that there is no function exported, so you should explicitly import functions if you need.

e.g.

```julia
import SPTK: mcep, mgcep, freqt, mc2b
```

## Supported Platforms

- Linux
- Mac OS X
- Windows

## Installation

```julia
Pkg.add("SPTK")
```

## Related packages

There are a few packages on top of SPTK.jl, which provide more Julia-like interface:

- [MelGeneralizedCepstrums.jl](https://github.com/r9y9/MelGeneralizedCepstrums.jl): Mel-Generalized Cepstrum analysis
- [SynthesisFilters.jl](https://github.com/r9y9/SynthesisFilters.jl): Speech waveform synthesis filters
