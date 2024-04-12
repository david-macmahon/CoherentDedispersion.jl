using CoherentDedispersion, Blio, PoolQueues, BlockArrays, FFTW, LinearAlgebra

#dir = "/datag/collate_mb/PKS_0277_2018-03-21T07:00/blc06"
dir = "/mnt_blpc3/datax/scratch/davidm/rawcodd"
rawfile = joinpath(dir, "guppi_58198_27514_685364_J0835-4510_B3_0001.0000.raw")
rawfiles = [rawfile]
DMVELA = 67.771

# Dispersion measure to dedisperse
dm = DMVELA

# nfpc is Number of Fine channels Per Coarse channel (set to 1 to disable)
# nint is the Number of consecutive fine spectra to INTegrate
#nfpc, nint = 1, 4 # Desired production values
nfpc, nint = 16, 64 # Test values

pqs, tasks = create_pipeline(rawfiles, dm; nfpc, nint)#, use_cuda=true)

fbname = fetch(last(tasks))
println("saved output to $fbname")
println("done")

##

# CPU-Only Pipeline:
#
#     Input Task (pushes overlapping blocks of RAW data)
#     |
#     |
#      > CPU Voltage PQ
#     |
#     V
#     CODD Task (CODD/upchan/detect)
#     |
#     |
#      > CPU Power PQ
#     |
#     V
#     Output Task (writes Filterbank files)

# CPU/GPU Pipeline:
#
#     Input Task (pushes overlapping blocks of RAW data)
#     |
#     |
#      > CPU Voltage PQ
#     |
#     V
#     GPU Copy Task
#                 |
#                 |
#                  > GPU Voltage PQ
#                 |
#                 V
#                 CODD Task (CODD/upchan/detect)
#                 |
#                 |
#                  > GPU Power PQ
#                 |
#                 V
#     GPU Copy Task
#     |
#     |
#      > CPU Power PQ
#     |
#     V
#     Output Task (writes Filterbank files)
