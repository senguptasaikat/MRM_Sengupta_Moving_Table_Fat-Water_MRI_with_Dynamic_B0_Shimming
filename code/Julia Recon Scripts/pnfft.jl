
# extension of NFFT to parallel workers

using NFFT

function pnfft_combi_worker(myz::Vector{Int}, data::Array{Complex64,3},
    K::Array{Float32,3}, pstart::Vector{Int},
    pend::Vector{Int}, kernwidth=3, oversamp=1.5f0; verbose=false)
  nchan, nro, npe = size(data)
  nz = length(pstart)
  n = int(nro / 2)
  fApprox = Array(Complex64, n, n, nz)
  sliceImgs = Array(Complex64, nchan, n, n)
  for z in myz   # for each "slice", recon 2D image
    if verbose
      println("slice $z/$nz, $nchan channels")
    end

    # initialize NFFT
    Kslice = K[:,:,pstart[z]:pend[z]]
    nfft_plan = NFFTPlan(Kslice[:,:], (n, n), kernwidth, oversamp)

    # SDC
    mynpe = pend[z] - pstart[z] + 1
    weights = radialsdc(nro, mynpe; golden=true, offset=pstart[z])
    #weights = sdc(nfft_plan) # Pipe & Menon SDC, not necessary for golden angle radial

    # apply NFFT plan to each channel
    for c in 1:size(data,1)
      # pull out the relevant profiles
      #fhat = vec(data[c,:,pstart[z]:pend[z]])
      mydata = data[c,:,pstart[z]:pend[z]]
      # radial phase correction
      fhat, slope, intercepts = unshift(mydata; golden=true)
      # density compensate
      fhat = fhat[:] .* weights[:]
      # transform and save the resulting image
      sliceImgs[c,:,:] = nfft_adjoint(nfft_plan, fhat)
    end
    # optimal channel combination
    fApprox[:,:,z] = walshcombine(sliceImgs)
  end
  fApprox
end


function pnfft_combi(data::Array{Complex64,3}, K::Array{Float32,3},
                     pstart::Vector{Int},
                     pend::Vector{Int}, kernwidth=3, oversamp=1.5f0;
                     verbose=false)

  # set up parallel stuff
  nz = length(pstart)
  nchan, nro, npe = size(data)
  n = int(nro / 2)
  nw = nworkers()
  zs = [1:nz]
  zperw = int(nz/nw)
  myz = Any[zs[(w - 1)*zperw + 1 : min(w*zperw, nz)] for w in  1:nw]
  println("processing $nz slices on $nw workers ($zperw slices each)")
  if nw == 1
    fApprox = pnfft_combi_worker(myz[1], data, K, pstart, pend,
                          kernwidth, oversamp; verbose=verbose)
  else
    reflist = Any[]
    wids = workers()
    for w in 1:nw
      println("spawning job on worker $(wids[w])")
      r = @spawnat wids[w] pnfft_combi_worker(myz[w], data, K,
             pstart, pend, kernwidth, oversamp; verbose=verbose)
      push!(reflist, r)
    end
    fApprox = Array(Complex64, n, n, nz)
    for w in 1:nw
      r = fetch(reflist[w])
      println("fetching zs $(myz[w])")
      fApprox[:,:,myz[w]] = r[:,:,myz[w]]
    end
  end
  fApprox
end
