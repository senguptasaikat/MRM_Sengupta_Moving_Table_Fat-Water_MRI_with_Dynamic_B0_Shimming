
function randcoords(nx, ny)
  K = zeros(Float32, 2, nx, ny)
  for j in 1:ny, i in 1:nx
    K[1,i,j] = rand() - 0.5
    K[2,i,j] = rand() - 0.5
  end
  K
end


function radialcoords(nr, narms; golden=false, offset=0)
  K = zeros(Float32, 2, nr, narms)
  for j in 1:narms
    θ = (j - 1 + offset) * π / (golden ? φ : narms)
    for i in 1:nr
      r = (i-1) / nr - 0.5
      K[1,i,j] = r*cos(θ)
      K[2,i,j] = r*sin(θ)
    end
  end
  K
end
function radialcoords(nr, pct::Float32=100.0; kwargs...)
  # can call with a radial percentage as well
  narms = iround(nr*pct / 100.0)
  radialcoords(nr, narms; kwargs...)
end


function radialcoords(nr, nprof, nleaves; pct=100.0, golden=false,
                      fid=false, roalt=false)
  nprof = iround(nprof*sqrt(pct/100.0))
  nleaves = iround(nleaves*sqrt(pct/100.0))
  K = zeros(Float32, 3, nr, nprof, nleaves)
  if golden
    α = 0.4656
    β = 0.6823
    ϕ = 0.0
    μ = 0.0
    for k in 1:nleaves, j in 1:nprof
      st = sqrt(1.0 - μ*μ)
      for i in 1:nr
        r = (i-1) / nr - 0.5
        K[1,i,j,k] = r*cos(ϕ)*st
        K[2,i,j,k] = r*sin(ϕ)*st
        K[3,i,j,k] = r*μ
      end
      ϕ = mod(ϕ + 2π*β, 2π)
      μ += 2α
      μ = μ > 1.0 ? μ - 2.0 : μ
    end
  else
    for j in 1:nprof
      if fid
        z = (2j - 1) / float(nprof)  - 1.0 # not sure about this
	z = j % 2 ? z : -z
      else
        z = (j - 0.5) / float(nprof) - 1.0
      end
      for k in 1:nleaves
        ϕ = sqrt(nprof * π / nleaves) * asin(z) + (k-1) * 2π / nleaves
	#from mmiffe_mxg.c phi = sqrt( `MP_radial_nr_angles * PI / `MP_radial_3d_koosh_interleaves ) * asin( z ) + interleaf * 2.0 * PI / `MP_radial_3d_koosh_interleaves;
        st = sqrt(1.0 - z*z)
        for i in 1:nr
          r = (i-1) / nr - 0.5
          K[1,i,j,k] = r*cos(ϕ)*st
          K[2,i,j,k] = r*sin(ϕ)*st
          K[3,i,j,k] = r*z
        end
      end
    end
  end
  reshape(K, 3, nr, nprof*nleaves)
end


function radialsdc(nr, nt; golden::Bool=false, offset::Int=0)
  # First guess radial sample density compensation
  # Uses Ram-Lak.  If golden angle, then the Ram-Lak
  # ramp is adjusted by the average neighbor distance.
  # sort angles in increasing order
  θ = (π / φ)*[offset:nt-1+offset] % (2π)
  Δθavg = 4π / nt
  θhat = sort(θ)
  isort = sortperm(θ)
  # find angular distance to neighbors
  Δθ = zeros(nt)
  Δθ[1] = θhat[2] - θhat[end] + 2π
  for t in 2:nt-1
    Δθ[t] = θhat[t+1] - θhat[t-1]
  end
  Δθ[end] = θhat[1] - θhat[end-1] + 2π
  weights = similar(Δθ)
  weights[isort] = Δθ / Δθavg
  r = linspace(1.0, 1.0/nr, int(nr/2))
  ramp = vcat(r, r[end:-1:1])
  if golden
    sdc = broadcast(*, ramp, weights')
  else
    sdc = repmat(ramp, 1, nt)
  end
  float32(sdc)
end
