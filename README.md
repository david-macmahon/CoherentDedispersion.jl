# CoherentDedispersion

This package provides coherent dedispersion and full cross-polarization or
Stokes detection of radio telescope voltage data.  The input format is GUPPI
RAW.  The output format is SIGPROC Filterbank.  This package can run with or
without a CUDA-enabled GPU, but use with a CUDA-enabled GPU is recommended for
higher throughput.

# Installation

This package and several of its dependent packages are not yet in the `General`
Julia registry.  They are available in a separate publicly accessible registry.
The easiest way to install this package is to add this registry to your Julia
"depot":

```bash
$ julia -e 'using Pkg; Registry.add(RegistrySpec(url="https://github.com/david-macmahon/MyJuliaRegistry"))'
```

After adding that registry, you can add the `CoherentDedispersion` package
using:

```bash
$ julia -e 'import Pkg; Pkg.add("CoherentDedispersion")'
```

If you would like a symlink to this package's command line utility to be created
in directory `/path/to/bin`, you can export the environment variable
`CODD_BINDIR` before adding the package or by "building" if it has already been
added:

```bash
env CODD_BINDIR="$HOME/bin" julia -e 'import Pkg; Pkg.build("CoherentDedispersion"; verbose=true)'
```

# Command line interface

This package include a bash script, `bin/rawcodd.jl` that provides a convenient
command line interface.  Creating a symlink to the script (see above) in a
directory in your path will allow you to run the script from anywhere without
having to specify the path to it.

Here is the command line help for `rawcodd.jl`:

```bash
$ rawcodd.jl --help
usage: rawcodd.jl -d DM [-f FFT] [-t INT] [-o OUTDIR] [-h] RAWFILES...

positional arguments:
  RAWFILES             GUPPI RAW files to process

optional arguments:
  -d, --dm DM          dispersion measure (type: Float64)
  -f, --fft FFT        up-channelization FFT length (type: Int64,
                       default: 1)
  -t, --int INT        spectra to integrate (type: Int64, default: 4)
  -o, --outdir OUTDIR  output directory (default: ".")
  -h, --help           show this help message and exit
```

Caveats:

* All files given on the command line will be treated as a single sequence of
  contiguous GUPPI RAW files.  No checks are performed to verify that that is
  actually the case.  Be sure to specify GUPPI RAW files accordingly.

* The list of file names will be sorted before being processed.  If you want to
  process a list of files in unsorted order, use the package from within Julia
  (see below).

* No checks are performed on the existence of the output file.  Existing output
  files will be silently overwritten.  If you want to save the output from the
  same data files with different options, be sure to use different output
  directories.

* This script can only process one scan at a time.  Each time it runs, CUDA
  initialization occurs, which imposes several seconds of delay.  If you have
  many similar scans to process with the same dispersion measure (DM), you may
  want to consider using this package from within Julia to be able to amortize
  the CUDA initialization over more scans.

* Output directories will be created as needed.

# Output filename

Currently the output filename is generated from `outdir` (defaults to `"."`) and
the `basename` of the first filename given.  The `.raw` extension of the input
file is removed (if present) and `.rawcodd.0000.fil` is concatenated.  The GUPPI
RAW sequence number (if any) is retained so that individual GUPPI RAW files from
the same scan may be processed individually without overwriting the same output
file.

# Programmatic interface

For more advanced use, one can use the package programmatically from within
Julia.  This allows for more control over the pipeline and the tasks that run
it.

## Reusable pipeline

This package creates data buffers that are sized specifically for a given GUPPI
RAW block size and other parameters.  These buffers and their associated FFT
plans are considered to be a *pipeline*.  A set of contiguous GUPPI RAW data
files (aka a *scan*) is processed by running the pipeline, which creates
asynchronous tasks for that particular set of contiguous GUPPI RAW data files.
It is possible to run a pipeline for multiple sets of contiguous GUPPI RAW data
files, one set after another.  A new set of tasks are started for each set of
contiguous input files.

The sizing of the data buffers depends on multiple factors:

1. The geometry of the GUPPI RAW blocks
2. The requested up-channelization and time integration factors
3. The maximum dispersive delay, which depends on the frequencies and dispersion
   measure being dedispersed

When reusing the pipeline on a new set of files, the new files are first checked
for compatibility with the existing pipeline.  If there is a mismatch an error
message is displayed and no processing will occur for that set of files.

The overall process is:

1. Create the pipeline for the first set of input files (i.e. for a scan)
2. Run the pipeline for each set of input files (i.e. for each scan)

## Creating the pipeline

The first step is to create the pipeline.  This is done using the aptly named
`create_pipeline` function.  The `create_pipeline` function needs information
about the GUPPI RAW data files, the dispersion measure to be dedispersed, the
number of fine channels per coarse channel (if up-channelization is desired),
and the number of time samples to integrate after detection.  Additional
optional keyword arguments can be used to increase buffering within the
pipeline and to explicitly opt out of using CUDA (e.g. for testing).  By
default, CUDA will be used if the local system is determined to be "CUDA
functional".

```julia
create_pipeline(rawinfo, dm; nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())
```

The `rawinfo` parameter can be a `GuppiRaw.Header` object, the name of a GUPPI
RAW file, or a Vector of names of GUPPI RAW files.  The latter two cases the
first GUPPI RAW header of the only/first file will be used.  Generally, a
Vector of GUPPI RAW names will be the most convenient to use.

## (Re-)using the pipeline

The pipeline object returned by `create_pipeline` can be passed to the
`start_pipeline` or `run_pipeline` functions along with a Vector of GUPPI RAW
filenames and an optional output directory (which defaults to the current
directory).  Both `start_pipeline` and `run_pipeline` start the asynchronous
tasks that perform the dedispersion process and create the output file.
Generally, `start_pipeline` is more versatile whereas `run_pipeline` can be
more convenient for single pipeline applications.

Both `start_pipeline` and `run_pipeline` can be called multiple times with the
same pipeline object but different sets of (compatible) input files.  As
described below, `start_pipeline` does not wait for completion so the caller is
responsible for waiting for completion before calling `start_pipeline` again.
`run_pipeline` does wait for completion so the caller may call `run_pipeline`
again as soon as it returns.

By default, `run_pipeline` displays a progress bar as processing progresses, but
`start_pipeline` does not.  This can be controlled explicitly by passing the
keyword argument `progress=true` or `progress=false` to either function.

Both `start_pipeline` and `run_pipeline` support additional keyword arguments
that control the behavior of the detection process:

* `dostokes` - `true`/`false` value indicating whether to compute Stokes
               parameters (`true`) or cross polarization products (`false`,
               default).
* `doconj` - `true`/`false`/`nothing` value indicating whether to negate Stokes
             V/conjugate the cross polarization products (`true`
             negates/conjugates, `false` does not) or to negate/conjugate only
             when `OBSBW` is negative (`nothing`, default).
* `doscale`  - `true`/`false` value indicating whether to scale the outputs to
               match the scaling performed by `rawspec` (another GUPPI RAW to
               Filterbank tool).  `true` (default) scales to match `rawspec`;
               `false` does not.

### Starting the pipeline

`start_pipeline` starts the tasks and returns a NamedTuple of the tasks without
waiting for the tasks to complete (i.e. the tasks will likely still be running
after `start_pipeline` returns).  The tasks in the NamedTuple are in "pipeline
order", so waiting for the last task (i.e. the `outputtask`) will wait for
completion of processing all the input files.  Calling `fetch` on the output
task will wait for completion and return the name of the output Filterbank file.

```julia
start_pipeline(pipeline, rawfiles; outdir=".", progress=false)
```

### Running the pipeline

`run_pipeline` is essentially `start_pipeline` plus a fetch of the last task.
`run_pipeline` returns the name of the output Filterbank file, but only after
processing is complete.

```julia
run_pipeline(pipeline, rawfiles; outdir=".", progress=false)
```

## Putting it all together

Here is a short script that shows how to use `CoherentDedispersion` to
dedisperse a list of GUPPI RAW files obtained from a not shown user-supplied
function) using a dispersion measure of `123.456`, up-channelizing by a factor
of `16`, and integrating `128` time samples (i.e. up-channelized spectra) after
detecting, and outputting a Filterbank file in the current directory.  This is
essentially a simplified version of the `rawcodd.jl` command line interface.

```julia
using CoherentDedispersion

rawnames = your_function_to_get_list_of_raw_files()
dm = 123.456

pipeline = create_pipeline(rawfiles, dm; nfpc=16, nint=64)
fbname = run_pipeline(rawfiles, dm)

@info "saved output to $fbname"
@info "done"
```
