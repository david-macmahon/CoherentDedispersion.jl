using CoherentDedispersion

#dir = "/datag/collate_mb/PKS_0277_2018-03-21T07:00/blc06" # Gluster
#dir = "/mnt_blpc3/datax/scratch/davidm/rawcodd" # NFS (self mounted)
dir = "/datax/scratch/davidm/rawcodd" # local XFS

rawfile = joinpath(dir, "guppi_58198_27514_685364_J0835-4510_B3_0001.0000.raw")
rawfiles = [rawfile]
DMVELA = 67.771

# nfpc is Number of Fine channels Per Coarse channel (set to 1 to disable)
# nint is the Number of consecutive fine spectra to INTegrate
#nfpc, nint = 1, 4 # Desired production values
nfpc, nint = 16, 64 # Test values

pqs, tasks = create_pipeline(rawfiles, DMVELA; nfpc, nint, outdir=".")

fbname = fetch(last(tasks))
@info "saved output to $fbname"
@info "done"
