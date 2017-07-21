__precompile__()

module SPTK

# A thin Julia wrapper for Speech signal Processing Toolkit (SPTK).
#
# The original SPTK: http://sp-tk.sourceforge.net/
# Note that the SPTK.jl is based on the slightly modified SPTK
# https://github.com/r9y9/SPTK
#
# NOTE1: no function is exported; i.e., you should explicitly import functions
# if you need.
#
# NOTE2: Function interfaces might be different betweeen C and Julia,
# as the following convensions.
#
# Conventions for wrapping C interface
#
# 1. Avoid really short names for variables (e.g. a, b, c, aa, bb, dd) .
#
#    Variable names should be informative. If the C functions have such short
#    names, use self-descriptive names instead for Julia interfaces, unless
#    they have clear meanings in their context.
#
# 2. Avoid too many function arguments.
#
#    Less is better. If the C functions have too many function arguments, use
#    keyword arguments with proper default values for optional ones in Julia.
#
# 3. Handle errors in Julia
#
#    Since C functions might `exit` (unfortunately) inside thier functions for
#    unexpected inputs, it should be check if the inputs are supported or not
#    in Julia before `ccall`.

deps = joinpath(dirname(@__FILE__), "..", "deps", "deps.jl")
if isfile(deps)
    include(deps)
else
    error("SPTK not properly installed. Please run Pkg.build(\"SPTK\")")
end

for fname in [
              "common",
              "lib",
              "adaptive",
              "conversions",
              "excite",
              "f0",
              "mfcc",
              "mgcep",
              "synthesis_filters",
              "utils",
              "window",
              "extend"
    ]
    include(string(fname, ".jl"))
end

end # module SPTK
