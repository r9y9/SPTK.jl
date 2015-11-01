### F0 estimation ###

const supported_otypes = ["pitch", "f0", "logf0"]

function swipe(x::StridedVector{Cdouble}, fs, hopsize=80;
               min::Float64=60.0,
               max::Float64=240.0,
               threshold::Float64=0.3,
               otype::Union{Int,AbstractString}="pitch")
    otype_int = (-1)::Int
    if isa(otype, AbstractString)
        otype_int = indexin([otype], supported_otypes)[1] - 1
    else
        otype_int = otype
    end

    if otype_int ∉ 0:2
        throw(ArgumentError("""
unsupported otype: $otype, must be (0) pitch, (1) f0, or (2) log(f0)"""))
    end

    if min >= max || max >= fs/2
        throw(ArgumentError("invalid min/max frequency parameters"))
    end

    expectedlen = ceil(Int, length(x) / hopsize)
    f0 = Array{Cdouble}(expectedlen)
    ccall((:swipe, libSPTK), Void,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Cdouble, Cdouble,
           Cdouble, Cint),
          x, f0, length(x), fs, hopsize, min, max, threshold, otype_int)
    f0
end

function rapt(x::StridedVector{Cfloat}, fs, hopsize=80;
              min::Float64=60.0,
              max::Float64=240.0,
              voice_bias::Float64=0.0,
              otype::Union{Int,AbstractString}="pitch")
    otype_int = (-1)::Int
    if isa(otype, AbstractString)
        otype_int = indexin([otype], supported_otypes)[1] - 1
    else
        otype_int = otype
    end

    if otype_int ∉ 0:2
        throw(ArgumentError("""
unsupported otype: $otype, must be (0) pitch, (1) f0, or (2) log(f0)"""))
    end

    if min >= max || max >= fs/2 || min <= fs/10000.0
        throw(ArgumentError("invalid min/max frequency parameters"))
    end

    frame_period = Float64(hopsize / fs)
    frame_period = trunc(0.5 + fs * frame_period) / fs
    if frame_period > 0.1 || frame_period < 1/fs
        throw(ArgumentError("frame period must be between [1/fs, 0.1]"))
    end

    expectedlen = ceil(Int, length(x) / hopsize)
    f0 = Array{Cfloat}(expectedlen)
    ret = ccall((:rapt, libSPTK), Int,
                (Ptr{Cfloat}, Ptr{Cfloat}, Cint, Cdouble, Cint, Cdouble, Cdouble,
                 Cdouble, Cint),
                x, f0, length(x), fs, hopsize, min, max, voice_bias, otype_int)

    if ret == 2
        error("input range too small for analysis by get_f0")
    elseif ret == 3
        error("problem in init_dp_f0()")
    end

    @assert ret == 0

    f0
end
