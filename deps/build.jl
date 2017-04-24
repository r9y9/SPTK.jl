using BinDeps

@BinDeps.setup

sptk = library_dependency("libSPTK", aliases=["libSPTK", "SPTK-3"])

# TODO: needs update for windows
const version = is_windows() ? "3.8.9" : "3.10.2"

github_root = "https://github.com/r9y9/SPTK"
arch = Sys.WORD_SIZE == 64 ? "x86_64" : "i686"
major = version[1]
provides(Binaries,
         URI("$(github_root)/releases/download/v$(version)/sptk-$(major)_mingw$(Sys.WORD_SIZE)_$(arch).zip"),
         sptk, unpacked_dir = "usr/lib", os = :Windows)

provides(Sources,
         URI("$(github_root)/archive/v$(version).tar.gz"),
         sptk,
         unpacked_dir="SPTK-$(version)")

prefix = joinpath(BinDeps.depsdir(sptk), "usr")
srcdir = joinpath(BinDeps.depsdir(sptk), "src", "SPTK-$(version)")

provides(SimpleBuild,
          (@build_steps begin
              GetSources(sptk)
              @build_steps begin
                  ChangeDirectory(srcdir)
                  `./waf configure --prefix=$prefix`
                  `./waf build`
                  `./waf install`
              end
           end), sptk, os = :Unix)

@BinDeps.install Dict(:libSPTK => :libSPTK)
