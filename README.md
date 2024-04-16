# CoherentDedispersion

This package provides coherent dedispersion and full Stokes detection of radio
telescope voltage data.  The input format is GUPPI RAW.  The output format is
SIGPROC Filterbank.  This package can run with or without a CUDA-enabled GPU,
but use with a CUDA-enabled GPU is recommended for higher throughput.

# Reusable pipeline

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
2. The requested upchannelization and time integration factors
3. The maximum dispersive delay, which depends on the frequencies and dispersion
   measure being dedispersed

When reusing the pipeline on a new set of files, the new files are first checked
for compatibility with the existing pipeline.  If there is a mismatch an error
message is displayed and no processing will occur for that set of files.

The overall process is:

1. Create the pipeline for the first set of input files (i.e. for a scan)
2. Run the pipeline for each set of input files (i.e. for each scan)

# Creating the pipeline

The first step is to create the pipeline.  This is done using the aptly named
`create_pipeline` function.  The `create_pipeline` function needs information
about the GUPPI RAW data files, the dispersion measure to be dedispersed, the
number of fine channels per coarse channel (if upchannelization is desired), and
the number of time samples to integrate after detection.  Additional optional
keyword arguements can be used to increase buffering within the pipeline and to
explicitly opt out of using CUDA (e.g. for testing).  By default, CUDA will be
used if the local system is determined to be "CUDA functional".

```julia
create_pipeline(rawinfo, dm; nfpc=1, nint=4, N=2, use_cuda=CUDA.functional())
```

The `rawinfo` parameter can be a `GuppiRaw.Header` object, the name of a GUPPI
RAW file, or a Vector of names of GUPPI RAW files.  The the latter two cases the
first GUPPI RAW header of the only/first file will be used.  Generally, a Vector
of GUPPI RAW names will be the most convenient to use.

# (Re-)using the pipeline

The pipeline object returned by `create_pipeline` can be passed to the
`start_pipeline` or `run_pipeline` functions along with a Vector of GUPPI RAW
filenames and an optional output directory (which defaults to the current
directory).  Both `start_pipeline` and `run_pipeline` start the asynchronous
tasks that perform the dedispersion process and create the output file.
Generally, `start_pipeline` is more versatile while `run_pipeline` can be more
convenient for single pipeline applications.

Both `start_pipeline` and `run_pipeline` can be called multiple times with the
same pipeline object but different sets of (compatible) input files.  As
described below, `start_pipeline` does not wait for completion so the caller is
responsible for waiting for completion before calling `start_pipeline` again.
`run_pipeline` does wait for completion so the caller may call `run_pipeline`
again as soon as it returns.

## Starting the pipeline

`start_pipeline` starts the tasks and returns a NamedTuple of the tasks without
waiting for the tasks to complete (i.e. the tasks will likely still be running
after `start_pipeline returns`).  The tasks in the NamedTuple are in "pipeline
order", so waiting for the last task (i.e. the `outputtask`) will wait for
completion of processing all the input files.  Calling `fetch` on the output
task will wait for completion and return the name of the output Filterbank file.

## Running the pipeline

`run_pipeline` is essentially `start_pipeline` plus a fetch of the last task.
`run_pipeline` returns the name of the output Filterbank file, but only after
processing is complete.

# Assembling lists of GUPPI RAW filenames

Creating the list of GUPPI RAW filenames depends on local directory structure
and file naming conventions, but typically GUPPI RAW files for multiple scans
are stored together in the same directory.  The files for a given scan typically
consist of a common *stem* plus a `.####.raw` extension, where the `####` part
is a four digit sequence number.

## Finding all `.0000.raw` files

Given a directory `dir`, the `.0000.raw` file (i.e. the first file of each
scan), can be found by using `filter` and `readdir`:

```julia
zfiles = filter(endswith(".0000.raw"), readdir(dir))
```

## Finding all files for a given `.0000.raw` file:
