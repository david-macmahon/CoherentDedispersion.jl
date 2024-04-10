using CoherentDedispersion, Blio, PoolQueues, BlockArrays, FFTW, LinearAlgebra

dir = "/datag/collate_mb/PKS_0277_2018-03-21T07:00/blc06"
rawfile = joinpath(dir, "guppi_58198_27514_685364_J0835-4510_B3_0001.0000.raw")
rawfiles = [rawfile]
DMVELA = 67.771

# Dispersion measure to dedisperse
dm = DMVELA

# nfpc is Number of Fine channels Per Coarse channel (set to 1 to disable)
# nint is the Number of consecutive fine spectra to INTegrate
#nfpc, nint = 1, 4 # Desired production values
nfpc, nint = 16, 64 # Test values

# Load data files (reads headers, mmap's data blocks)
# Use explicitly typed Vector for datablocks so we don't have to refine later.
hdrs, blks = GuppiRaw.load(rawfiles; datablocks=Array{Complex{Int8},3}[]);

npol, _, nchan = size(blks[1])
npol == 2 || error("only dual-pol files are supported")
@show length(hdrs) length(blks) size(blks[1])

# Get freq info from header
obsfreq = hdrs[1][:obsfreq]
obsbw = hdrs[1][:obsbw]
# Hack to correct for incorrect sign on obsbw.
if obsbw > 0
    @warn "obsbw > 0, negating to correct"
    obsbw = -obsbw
end

# Concatenate all blocks into a BlockArray "view" along the time (i.e. second)
# dimension
#data = mortar(reshape(convert(Vector{typeof(blks[1])}, blks), (1,:,1)));
data = mortar(reshape(blks, (1,:,1)));

# Compute ntime values
ntpi, ntpo = compute_ntimes(hdrs[1], dm; nfpc, nint)

cvb = CODDVoltageBuffer(ntpi, nfpc, nint, ntpo, nchan)
cpb = CODDPowerBuffer(nfpc, nchan, ntpo)

codd_plan  = plan_fft!(cvb.inputs[1], 1)
upchan_plan = plan_fft!(cvb.upchans[1], 1)

##

# Get parameters for H!
f0j = obsfreq - obsbw/2
dfj = obsbw / (hdrs[1].obsnchan / get(hdrs[1], :nants, 1))
dm = DMVELA

# Copy raw block 3 to cvb
copyraw!(cvb, blks[3])

# Coherently dedisperse
foreach(pol->mul!(pol, codd_plan, pol), cvb.inputs)
H!(cvb, f0j, dfj, dm)
foreach(pol->mul!(pol, inv(codd_plan), pol), cvb.inputs)

# Upchannelize
foreach(pol->mul!(pol, upchan_plan, pol), cvb.upchans)

# Detect
detect!(cpb, cvb);

##

# CPU-Only Pipeline:
#
#     Read Task (reads RAW files)
#     |
#     |
#      > CPU Voltage PQ
#     |
#     V
#     CODD Task (CODD/upchan/detect)
#     |
#     |
#      > GPU Power PQ
#     |
#     V
#     Output Task (writes Filterbank files)

# CPU/GPU Pipeline:
#
#     Read Task (reads RAW files)
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

#=
# Create PoolQueue for input.  Pool holds NTuple{2, Array{ComplexF32,2}}, i.e.
# tuple of 2 matrices (one per polarization).
pqin = CoherentDedispersion.CODDInputArrayPQ(4) do
    ntuple(i->zeros(ComplexF32, ntpb_in, nchan), 2)
end

# Create PoolQueue for output.  Pool holds NTuple{2, Array{ComplexF32,2}}, i.e.
# tuple of 2 matrices (one per polarization).
pqout = CoherentDedispersion.CODDOutputArrayPQ(4) do
    ntuple(i->zeros(ComplexF32, ntpb_out, nchan), 2)
end

# coddbufs are "permanent" (i.e. not PQ-managed) work buffers for coddtask
coddbufs = ntuple(i->zeros(ComplexF32, ntpb_in, nchan), 2)

# outbufs are "permanent" (i.e. not PQ-managed) work buffers for outputtask
outbufs = ntuple(i->zeros(ComplexF32, ntpb_out, nchan), 2)

coddtask = errormonitor(
    Threads.@spawn CoherentDedispersion._coddtask(pqin, pqout; coddbufs)
)

outputtask = errormonitor(
    Threads.@spawn CoherentDedispersion._outputtask(pqout; outbuf)
)

dm = 557

# Produce blocks
for (hdr, blk) in zip(hdrs, blks)
    produce!(pqin) do inputbufs
        # bogus hack to avoid depleting all pqin elements (still needed?)
        sleep(0.1)

        # inputbufs is a tuple of 2 single-pol arrays
        #copyto!.(inputbufs, (@view(blk[i,:,:]) for i in 1:2))
        copyto!.(inputbufs, eachslice(blk, dims=1))
        (; dm, rawname, header=hdr, data=inputbufs)
    end
end

# Signify end of input
produce!(pqin, (;))

#isready(pqout.queue)
#
#nt1 = consume!(pqout)
#recycle!(pqout, nt1.data)
#
#nt2 = consume!(pqout)
#recycle!(pqout, nt2.data)
=#