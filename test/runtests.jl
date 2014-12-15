using SPTK
using Base.Test

for fname in [
              "call",
              "mgcep",
              "mfcc"
    ]
    include(string(fname, ".jl"))
end
