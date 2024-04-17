#!/bin/bash
#=
export JULIA_PROJECT=$(dirname $(dirname $(readlink -e "${BASH_SOURCE[0]}")))
exec julia --color=yes --startup-file=no --threads=auto \
    -e 'include(popfirst!(ARGS))' "${BASH_SOURCE[0]}" "$@"
=#

using ArgParse, CoherentDedispersion

# Parse command line.
function parse_commandline()
    aps = ArgParseSettings()
    @add_arg_table! aps begin
        "--dm", "-d"
            help = "dispersion measure"
            required = true
            arg_type = Float64
        "--fft", "-f"
            help = "up-channelization FFT length"
            default = 1
            arg_type = Int
        "--int", "-t"
            help = "spectra to integrate"
            default = 4
            arg_type = Int
        "--outdir", "-o"
            help = "output directory"
            arg_type = String
            default = "."
        "RAWFILES"
            help = "GUPPI RAW files to process"
            required = true
            nargs = '+'
            arg_type = String
    end
    return parse_args(aps, as_symbols=true)
end # function parse_commandline

function main()
    args = parse_commandline()

    rawfiles = args[:RAWFILES]
    dm = args[:dm]
    nfpc = args[:fft]
    nint = args[:int]
    outdir = args[:outdir]

    # Ensure all files exist
    if any(!isfile, rawfiles)
        error("some input files seem to be missing (or not files)")
    end

    # Ensure outdir exists
    mkpath(outdir)

    # Sort input file list
    sort!(rawfiles)

    pipeline = create_pipeline(rawfiles, dm; nfpc, nint)
    fbname = run_pipeline(pipeline, rawfiles; outdir)

    @info "saved output to $fbname"
    @info "done"

    return 0
end

main()
