
# Continuous Moving Bed Reconstruction 
#
# All data and results are in single-precision (32-bit) float, to
# reduce memory requirements.

using NFFT
using MAT
include("correct.jl")
require("Julia Recon Scripts/coords.jl")
require("Julia Recon Scripts/pnfft.jl")



function coilcompress{T}(x::Array{T,3}; nvirt=6)
  # Geometric coil compression with alignment to middle
  println("compressing data to $nvirt virtual coils")
  nchan, nro, npe = size(x)
  y = Array(T, nvirt, nro, npe)
  n = div(npe,2)  # middle profile, the starting point
  U0 = Array(T, nchan, nvirt)
  Uprev = Array(T, nchan, nvirt)
  for p in n
    t = squeeze(x[:,:,p], 3)
    C = t*t'
    U, s, V = svd(C)
    U0[:] = U[:,1:nvirt]
    y[:,:,p] = U0'*t
  end
  Uprev[:] = U0
  for p in n+1:npe
    t = squeeze(x[:,:,p], 3)
    C = t*t'
    U, s, V = svd(C)
    U = U[:,1:nvirt]
    W, s, V = svd(U'*Uprev)
    U = U*W*V'  # rotate U to align with previous slice
    y[:,:,p] = U'*t
    #Uprev[:] = U  # took this out bc seems unnecessary
  end
  Uprev[:] = U0
  for p in n-1:-1:1
    t = squeeze(x[:,:,p], 3)
    C = t*t'
    U, s, V = svd(C)
    U = U[:,1:nvirt]
    W, s, V = svd(U'*Uprev)
    U = U*W*V' # rotate U to align with previous slice
    y[:,:,p] = U'*t
    #Uprev[:] = U # took this out bc seems unnecessary
  end
  y
end



function test_coilcompress()
  # ADJUSTABLE PARAMETERS
  golden = true # this should correspond to the sequence
  npe_per_z = 128
  sliceΔz = 64
  fcc_phase = 1.4580  # from SIN file
  fcc_scale = 1.1784
  infile = "1305_QBC_GA_12mm_15_3.7.mat"
  #infile = "/Users/dss/Dropbox/1305_QBC_GA_12mm_15_3.7.mat"
  #infile = "/Users/dss/Data/combi/CMT_XTEND_20141216/1219_QBC_GA_12mm_15_3.7_min.mat"
  infile = "/home/dss/Data/combi/CMT_XTEND_20141216/1231_QBC_GA_12mm_15_3.7_min.mat"

  # load some COMBI data
  println("reading data from $infile")
  mat = matread(infile)
  data = mat["data"]
  info = mat["info"]  # in case we want to pull params from INFO
  data = coilcompress(data)
end


function walshcombineworker{T,N}(x::Array{T,N})
  # Performs Walsh adaptive coil combination on arbitrary N-d arrays
  # Assumes:
  #   first dimension is channel dimension
  #   number of dimensions N is at least 2
  @assert N > 1 "input must have two or more dimensions"
  insize = size(x)
  nchan = insize[1]
  nvox = prod(insize[2:end])
  W = Array(T, nchan, nvox)
  tstart = time()
  y = Array(T, nvox)
  for i in 1:nvox
    #if mod(i,div(nvox,20)) == 0
    #  pct = 100i / float(nvox)
    #  telapsed = time() - tstart
    #  trem = 100telapsed/pct - telapsed
    #  @printf "channel combination %.0f%% complete, %.1f s elapsed, %.1f s remaining\n" pct telapsed trem
    #end
    C = x[:,i]*x[:,i]'
    U, s, V = svd(C)
    y[i] = (U[:,1]'*x[:,i])[1]
  end
  return y
end


walshcombine{T}(x::Array{T,2}) = walshcombineworker(x) # single threaded
walshcombine{T}(x::Array{T,3}) = walshcombineworker(x) # single threaded

function walshcombine{T}(x::Array{T,4})  # multi-threaded
  # Wraps Walsh combination code for volume data and manages parallel execution
  # Assumes:
  #   First dimension is the channel dimension.
  #   All workers are available for work.
  nchan, nx, ny, nz = size(x)
  if nchan == 1  # skip fancy combination if only one channel
    return squeeze(x[1,:,:,:], 1)
  end
  y = Array(T, nx, ny, nz)
  nw = nworkers()  # see how many workers available
  if nw > 1
    blas_set_num_threads(1) # don't compete with BLAS threading
    reflist = Any[]
    workerids = workers()
    slices = [[w:nw:nz] for w in 1:nw] # divvy up data by slices
    for w in 1:nw
      println("worker $w gets slices $(slices[w])")
      r = @spawnat workerids[w] walshcombineworker(x[:,:,:,slices[w]])
      push!(reflist, r)   # save the ref for later 'fetch'ing
    end
    y = Array(T, nx, ny, nz)
    for w in 1:nw  # for each worker
      y[:,:,slices[w]] = fetch(reflist[w])  # get result and store
    end
  else  # run single threaded
    y = walshcombine(x)
  end
  y
end

