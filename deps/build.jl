if haskey(ENV, "CODD_BINDIR")
    link = joinpath(ENV["CODD_BINDIR"], "rawcodd.foo")
    target = joinpath(dirname(@__DIR__), "bin", "rawcodd.jl")

    if isfile(link)
        @info "replacing symlink $link -> $target"
        rm(link)
        symlink(target, link)
    elseif !ispath(link)
        @info "creating symlink $link -> $target"
        symlink(target, link)
    else
        @info "$link already exists as non-file"
    end
end
