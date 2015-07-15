
# Continuous Moving Bed Reconstruction using Gridding
#
# All data and results are in single-precision (32-bit) float, to
# reduce memory requirements.

using NFFT
using MAT
include("correct.jl")
require("Julia Recon Scripts/coords.jl")
require("Julia Recon Scripts/pnfft.jl")
require("Julia Recon Scripts/coilcomb.jl")


data_path = "../data_input/";
base_name = ARGS[1]
infile = string(data_path,base_name ,".mat")
out_dir = "../data_output/";

println(ARGS)



# ADJUSTABLE PARAMETERS
golden = true # this should correspond to the sequence
npe_per_z = int(ARGS[4])
sliceΔz = int(float(ARGS[3])/(20 * float(ARGS[2])))

nvirt=2  # number of virtual coils to compress to before the recon
fcc_phase = 1.4580  # from SIN file
fcc_scale = 1.1784



# load some COMBI data
println("reading data from $infile")
mat = matread(infile)
data = mat["Data"]
#info = mat["info"]  # in case we want to pull params from INFO
#noise = info["FRC_NOISE_DATA"][1]
#npe = size(data,3)  # reduce data for testing
#data = data[:,:,int(0.4npe):int(0.5npe)]


function recon_combi_ms2d(data, npe_per_z, sliceΔz; kernwidth=3,
                          oversamp=1.5f0, golden=false)

  nchan, nro, npe = size(data)
  #n = int(nro / 2)  # reconstructed image dimension
   n = int(nro)  # reconstructed image dimension
  pstart = [1:sliceΔz:npe-npe_per_z]
  pend = pstart + npe_per_z - 1
  nz = length(pstart)
   
  pcenter = int(pstart+npe_per_z/2-1)
  matwrite(string(out_dir,"pcenter.mat"),["pcenter"=>pcenter])
 

  # combine QBC channels data with Fast Channel Combination
  if nchan == 2
    println("performing FCC with scale=$fcc_scale and phase $fcc_phase")
    data = data[1,:,:] + fcc_scale*exp(-1im*fcc_phase)*data[2,:,:]
  else
    #prewhiten!(data, noise)
    data = coilcompress(data; nvirt=4)
  end

  # other corrections?
  tic()
  data, slope, intercepts = unshift(data; golden=golden)
  toc()

  # Generate coordinates for all profiles up front
  # TODO figure out why I need offset = 1 here
  println("generating radial coordinates")
  K = radialcoords(nro, npe; golden=golden, offset=0)

  # Generate the sample density compensation
  # Should be the same for each slice as long as they have
  # equal numbers of profiles.
  #println("computing SDC")
  #weights = radialsdc(nro, npe_per_z; golden=golden)

  println("reconstructing ... $(size(data))")
  fComb = pnfft_combi(data, K, pstart, pend, kernwidth, oversamp; verbose=true)

  return fComb
end



tic()
fComb = recon_combi_ms2d(data, npe_per_z, sliceΔz; golden=golden)
toc()

#println("saving to $outfile")
#matwrite(outfile, ["data"=>fComb])
matwrite(string(out_dir,base_name,"_JULIA_RECON.mat"),["Image"=>fComb])



#using PyPlot

#nx, ny, nz = size(fComb)
#figure(1)
#clf()
#imshow(abs(squeeze(fComb[:,:,int(0.3nz)],3)), interpolation="nearest", cmap="gray")
#title("axial")
#colorbar()

#figure(2)
#clf()
#aspect_ratio = (0.4/1.5)*nz/nx
#imshow(abs(squeeze(fComb[:,110,:],2)),aspect=aspect_ratio, cmap="gray", interpolation="nearest")
#title("sag or cor")
#colorbar()

#figure(7)
#clf()
#plot(p.windowLUT[1])
#title("kernel")
