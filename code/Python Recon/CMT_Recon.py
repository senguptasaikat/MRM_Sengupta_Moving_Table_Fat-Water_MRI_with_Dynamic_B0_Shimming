#!/usr/bin/env python

from pylab import *
from scipy.io import loadmat, savemat
from scipy.fftpack import fftshift, ifftn, fftn, fft, ifft2
from scipy.linalg import lstsq
from skimage.filter import threshold_otsu
from numpy.linalg import pinv
import numpy as np
import os
import sys
import math
import time
#import pdb    # PYTHON DEBUGGER, uncomment if debugging and add pdb.set_trace() to debug point.

######################################################################################
## CONFIGURATION SECTION
######################################################################################
## Modify these to match your scan parameters
## CHANGE ONLY base_name AND npes_per_slice TO RECONSTRUCT ANY DATA VOLUME TO ANY THICKNESS

# list of reconstructions to perform
recon_list = []

# Subject 1 No Shim
recon_list.append( { 'base_name' : 'Sub1_NS_echo1', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' :  int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_NS_echo2', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_NS_echo3', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_NS_echo4', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )

# Subject 1 Dynamic Shim
recon_list.append( { 'base_name' : 'Sub1_DS_echo1', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' :  int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_DS_echo2', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_DS_echo3', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub1_DS_echo4', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )

# Subject 2 No Shim
recon_list.append( { 'base_name' : 'Sub2_NS_echo1', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' :  int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_NS_echo2', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_NS_echo3', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_NS_echo4', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )

# Subject 2 Dynamic Shim
recon_list.append( { 'base_name' : 'Sub2_DS_echo1', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' :  int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_DS_echo2', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_DS_echo3', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )
recon_list.append( { 'base_name' : 'Sub2_DS_echo4', 'TR': 6.6e-3, 'rad_percentage' : 100.00, 'GA_Flag': True, 'npes_per_slice' : int(128), 'npe_per_180' :  int(256) } )


######################################################################################

######################################################################################
## HELPER FUNCTIONS FOLLOW
######################################################################################

class kernel_kb(object):
    def __init__(self, width=8, oversampling_factor=2, table_length=512):
        self.width = width
        self.radius = width/2
        self.table_length = table_length
        self.alpha = oversampling_factor
        self.beta = pi * sqrt(self.width**2 * (self.alpha - 0.5)**2 / self.alpha**2 - 0.8)  # from Beatty et al.
        self.table = kaiser(2*self.table_length, self.beta).astype(float32)
        self.table = self.table[self.table_length:]
        self.table[-1] = 0.0
        self.u = linspace(0, self.radius, self.table_length)
    def sample(self, K, subscripts):
        """ look up kernel at a given spherical radius """
        ndata, nhood, ndims = subscripts.shape
        rxy = zeros((ndata, nhood), dtype='float32')
        #print 'finding distances'
        for j in range(nhood):
             rxy[:,j] = (K[:,0] - subscripts[:,j,0])**2 +  (K[:,1] - subscripts[:,j,1])**2
        rxy = sqrt(rxy)
        #print 'transforming to data weights'
        rxy[rxy > self.radius] = self.radius
        rxy = (rxy * (self.table_length - 1) / self.radius).astype(int)
        return self.table[rxy]
    def plot(self):
        plot(self.u, self.table, 'k.-')
        title(self)
    def rolloff(self, G):
        ro = G.T * array([1.])   # grid unit datum
        ngrid = int(sqrt(G.shape[1]))
        ro = fftshift(abs(ifft2(reshape(ro, (ngrid, ngrid)))))  # transform to image
        ro = ro / ro.max()  # normalize
        return ro



def sparse_matrix_operator(K, grid_shape, osfactor=2, nhood_shape=None, verbose=False, combi=False, slice_prof_coef=None):
    from scipy.sparse import coo_matrix

    # SANITY CHECKS
    ndata, ndims = K.shape
    assert(type(K) == ndarray)

    if grid_shape is not None:
        n = int(osfactor*grid_shape[0])
    else:
        n = 2*(int(osfactor * abs(K[:,:2]).max()))    # find largest spatial frequency in coordinates and fit grid to that;
                                 # no point in a larger grid, right?
    grid_shape = array([n, n])


    if combi:
        if verbose:
            print 'initializing COMBI gridding kernel'
        kern = kernel_combi(oversampling_factor=osfactor, slice_prof_coef=slice_prof_coef)
    else:
        kern = kernel_kb(oversampling_factor=osfactor)
    radius = 0.5 * kern.width
    if nhood_shape == None:
        nhood_shape = tile(2*radius+1, ndims).astype(int)
    nhood = prod(nhood_shape)

    if verbose:
        print 'reserving memory blocks'
        print 'ndims: %d\nndata: %d\nnhood: %d' % (ndims, ndata, nhood)
        print 'grid shape:', grid_shape
        print 'neighborhood shape:', nhood_shape

    subscripts = zeros((ndata, nhood, ndims), dtype='float32')
    gcols = zeros((ndata, nhood), dtype='int')

    if verbose: print  'shifting frequencies to match Cartesian grid indexing'
    #nhood_shape = nhood_shape.astype(float32)
    kmax = abs(K).max()
    #print 'kmax=',abs(K).max()
    #print 'kmin=',abs(K).min()
    if kmax == 0:
        kmax = 1
    #K = K + (n / 2)
    K = (K / kmax + 1) * (n / 2)
    #print 'KERNEL WIDTH BORDER?'
    #print 'ROLLOFF may be off by 1'
    #print 'kmax=',abs(K).max()
    #print 'kmin=',abs(K).min()

    if verbose: print 'computing neighboring grid subscripts'
    # MAKE A TEMPLATE PATCH that will be added to each center freq (Kctr)
    # duplicate offsets NDIMS times
    offsets = [arange(-nhood_shape[d]/2+1, nhood_shape[d]/2+1) for d in range(ndims)]
    # pass meshgrid split arg version of offsets array using '*'
    patches = array(meshgrid(*offsets, indexing='ij'))
    # stretch out patches into linear arrays for each dimension
    patches = reshape(patches, (ndims, nhood)).T.copy()
    for j in xrange(nhood):
        subscripts[:,j,:] = K + patches[j,:]
    subscripts = floor(subscripts)
    del patches  # recover RAM from patches

    if verbose: print 'computing weights'
    weights = kern.sample(K, subscripts)
    if combi:
        weights[subscripts[:,:,2] < 0] = 0.0
        weights[subscripts[:,:,2] >= grid_shape[2]] = 0.0
    weights = weights.flatten()

    if verbose: print 'defining sparse coordinates'
    subscripts = subscripts % grid_shape    # periodically wrap out-of-bound indices
    # DO WE NEED TO CONJUGATE OOB POINTS?  CORRECT FOR PHASE?
    subscripts = reshape(subscripts, (ndata*nhood, ndims)).astype(int32).T

    # filter out zero weights
    mask = weights > 0.05
    weights = weights[mask]

    gcols = ravel_multi_index(tuple(subscripts), tuple(grid_shape)).flatten()
    gcols = gcols[mask]
    del subscripts
    grows = tile(arange(ndata, dtype='int'), (nhood, 1)).T.flatten()
    grows = grows[mask]

    if verbose:
        print 'filling in matrix ...'
        print 'grows:', grows.shape
        print 'gcols:', gcols.shape
        print 'wgts:', weights.shape

    G = coo_matrix((weights, (grows, gcols)), shape=(ndata, prod(grid_shape)))
    del gcols, grows, weights
    G = G.tocsr()

    if verbose:
        print 'G:' ,G.shape

    return G

def rolloff(grid_shape, nhood_shape=[6,6], ncrop=None, osfactor=2, threshold=0.01, axes=[-1, -2], combi=False, slice_prof_coef=1):
    if combi:
        nz = spatial_dims[-1]
        spatial_dims = grid_shape[:-1]
        ndims = len(spatial_dims)
        ro_data = ones(nz, dtype='complex64')
        ro_K = zeros((nz, ndims), dtype='float32')
        ro_K[:,2] = linspace(-nz/2, nz/2, nz, endpoint=False)
    else:
        spatial_dims = grid_shape
        ndims = len(spatial_dims)
        ro_data = array([1.0], dtype='complex64')
        ro_K = array([[0]*ndims], dtype='float32')
    G = sparse_matrix_operator(ro_K, grid_shape=spatial_dims, nhood_shape=nhood_shape, osfactor=osfactor, combi=combi, slice_prof_coef=slice_prof_coef)
    ro = G.T * ro_data
    n = sqrt(G.shape[1])
    ro = reshape(ro, (n,n))
    #if osfactor > 1 and ncrop == None:
    #    ncrop = int((spatial_dims[0]/osfactor) * (osfactor - 1) / 2)
    ro = fftshift(abs(ifftn(fftshift(ro, axes=axes), axes=axes)), axes=axes)  # transform to image
    if ncrop > 0:
        ro = ro[ncrop:-ncrop, ncrop:-ncrop]
    #print 'rolloff shape:', ro.shape
    ro = ro / ro.max()  # normalize
    ro[ro < threshold] = 1.0
    #ro_max = ro.max()
    #ro[ro < threshold*ro_max] = 1.0
    ro = 1.0 / ro
    #ro = ro**2
    #print 'TOOK OUT SQUARED RO'
    #ro = ro / ro.max()
    return ro

def sample_density(G, niter=1):
    """ iterative sample density compensation """
    weights = ones(G.shape[0], dtype='float32')
    for t in range(niter):
        weights = weights / (G * (G.T * weights))
    return weights

def do_fast_channel_combination(data1, data2, scale=1.0, phase=0.0):
    """ Fast Channel Combination (FCC) for Quadrature Body Coil """
    return data1 + (scale)*exp(-1j * phase) * data2

def do_kspace_shift_correct(kspace, golden_angle=False):
    # calculate sinogram
    sino = fftshift(ifft(fftshift(kspace, axes=1), axis=1), axes=1)
    nKr = kspace.shape[1]
    nAngles = kspace.shape[0]
    kspace_shifted = np.zeros((nAngles,nKr), dtype=np.complex64)
    slope_v = []
    intercept_v = []
    x1 = np.arange(-nKr/2, nKr/2)
    #angles = arange(nAngles)*111.246 % 360.0
    if golden_angle:
       # deltak = 21
        deltak = 34
        print 'angular distance neighbors is %f deg' % (deltak*111.246 % 360.0)
        print 'using deltak of +/- %d' % deltak

    for n in range(nAngles):
        if golden_angle:
            # taking 5th back and forward, should be delta of 16.25 deg
            if n >= deltak:
                a = sino[n-deltak,:]
            else:
                a = None
            b = sino[n,:]
            if n <= nAngles - 1 - deltak:
                c = sino[n+deltak,:]
            else:
                c = None
        else:  # linear
            if n == 0:
                a = sino[nAngles-1,:]
            else:
                a = sino[n-1,:]
            b = sino[n,:]
            if n == nAngles - 1:
                c = sino[0,:]
            else:
                c = sino[n+1,:]

        # average spatial projections a and c
        if a is None:
            a = c
            ac = c
        elif c is None:
            c = a
            ac = a
        else:
            ac = 0.5 * (a + c)
                
        # create a mask for the spatial domain projections
        abc_mag = np.abs(a) + np.abs(b) + np.abs(c)
        abc_mag = abc_mag - abc_mag.min()
        abc_mag = abc_mag / abc_mag.max()
        abc_mag_thresholded = abc_mag > threshold_otsu(abc_mag)

        # calculate masked phase difference
        if golden_angle:
       		phase_difference_ac_b = np.angle(ac[::-1] / b) * abc_mag_thresholded
        else:
            phase_difference_ac_b = np.angle(ac / b) * abc_mag_thresholded

        # perform a weighted least squares fit
        y = phase_difference_ac_b
        W = np.diag(abc_mag_thresholded)
        M = np.vstack([x1, np.ones(nKr)]).T
        # self.log.warn('%s %s %s' % (M.shape, W.shape, y.shape))
        u = pinv(M.T.dot(W.dot(M))).dot(M.T.dot(W.dot(y)))
                
        slope = u[0]
        intercept = u[1]
        slope_v.append(slope)
        intercept_v.append(intercept)

        # apply phase compensation
        phase_ramp_ini = slope * x1.flatten() + intercept
        phase_ramp = phase_ramp_ini / 2.0
        b_phase_compensated = b * np.exp(1j*phase_ramp)
        B_shifted = fftshift(fft(fftshift(b_phase_compensated)))
        kspace_shifted[n,:] = B_shifted

    return kspace_shifted

def do_zero_dc_phase(data):
    npe, nro = data.shape
    for k in range(npe):
        phase = np.angle(data[k,nro/2])
        data[k,:] = data[k,:] * np.exp(-1j*phase)
    return data

def interpolate_radial_data(data, slice_data, slice_pes, npe_per_180):
    
    # Weighted nearest-neighbor interpolation for linear radial COMBI data 
    slice_data_interp = empty_like(slice_data)
    pe_center = len(slice_pes) / 2
    for k, p in enumerate(slice_pes):
    
        if k < pe_center:  # interpolate with future projection
            slice_data_interp[:,k,:] = (pe_center - k)*data[:,p+npe_per_180,::-1] + \
                (k + npe_per_180 - pe_center)*slice_data[:,k,:]
        else:  # interpolate with past projection 
         	slice_data_interp[:,k,:] = (k - pe_center)*data[:,p-npe_per_180,::-1] + \
                abs(pe_center - (k - npe_per_180))*slice_data[:,k,:]
    return slice_data_interp / float(npe_per_180)

def Aop(X, mode, G, W, ro, ncrop, grid_shape, data_shape):
    if mode == 1:   # perform non-uniform sampling  (image/sparse -> data)
        Y = zeros(data_shape, 'complex64')
        X = pad(X, ((ncrop, ncrop), (ncrop, ncrop), (0,0)), mode='constant')
        X = fftshift(fftn(fftshift(X, axes=(0,1)), axes=(0,1)), axes=(0,1))
        npe, nro, nchan = X.shape
        X = reshape(X, (npe*nro, nchan))
        Y = G * X
    elif mode == 2:  # reconstruct images  (data -> image/sparse)
        Y = G.T * (X.T * W).T
        n = sqrt(G.shape[1])
        Y = reshape(Y, (n, n, -1))
        Y = (fftshift(ifftn(fftshift(Y, axes=(0,1)), axes=(0,1)), axes=(0,1)))
        Y = Y[ncrop:-ncrop, ncrop:-ncrop]
        Y = (Y.T * ro.T).T
    return Y

######################################################################################

######################################################################################
## PERFORM RECONSTRUCTIONS
######################################################################################

# DETECT LOCATION OF THIS SCRIPT
pathstr_python_script =  os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))

# Loop over all reconstructions
recon_cnt = 0
for recon in recon_list:

    
    recon_cnt +=1
    print("Performing reconstruction %02d of %02d" % (recon_cnt, len(recon_list) ) )    
    
    
    # DATA_FILE TO BE RECONSTRUCTED (Without .mat extension)
    base_name = recon['base_name']

    # NUMBER OF PROFILES TO USE PER SLICE (We use 128 here)  
    npes_per_slice = recon['npes_per_slice']  
                	 
    # path to input data file              
    data_file = pathstr_python_script + '/' + os.pardir +'/data_input/' + base_name + '.mat'

    # EXAM CARD PARAMETERS
    TR = recon['TR']
    vtable = 20.0   # mm/s
    voxel_size_x_mm = 2.5
    radial_pct = recon['rad_percentage']
    npe_per_180 = recon['npe_per_180']    #  Interpolation profile step, matches npes_per_slice for LA. Is not used for GA	
    nslices = int(720)  # number of slices in z direction
    pe_of_z0 = int(128)     # Definition of z=0 to occur at profile 128
    osfactor = 2   # gridding oversampling factor (1.5-2.0 is best)

    # COIL COMBINATION FILE PARAMETERS
    fcc_phase = 1.4580  
    fcc_scale = 1.1784

    # RECON FLAGS
    plotting = False                 	# do you want a plot of the intermediate results?
    golden_angle = recon['GA_Flag']    # was data golden angle or linear?
    fast_channel_combination = True  	# combine QBC channels before the recon?
    crop_ro_oversampling = False    	# crop the 2X readout oversampling before recon?
    kspace_shift_correct = True      	# correct for k-space shifts?
    zero_dc_phase = False            	# zero the DC phase after k-space shift correction?
    readout_alternation = False     	# are the profiles collected with readout alternation on? 
    apply_phase_ramp = False        	# set to True to apply a 1-pixel phase ramp
                                       # MAT FILE HAS THIS ALREADY CORRECTED, THIS IS NOT NEEDED FOR LINEAR ANGLE!, 
								    
    #
    # DO THE RECON
    #
    
    t = time.time()
    
    # read the data
    mat = loadmat(data_file)
    #try:
    #    mat = loadmat(sys.argv[0])
    #except IOError:
    #    print 'could not open', sys.argv[0]
    data = mat['Data'].astype('complex64').T.copy('C')
    npe, nro, nchan = data.shape
    print ' READING FILE : %s'% (data_file)
    print 'DATA INFO:\nnchan: %d\nnpe: %d\nnro: %d' % (nchan, npe, nro)
    
    
    if fast_channel_combination:
        print 'combining QBC channels with fast channel combination'
        data = do_fast_channel_combination(data[:,:,0], data[:,:,1], scale=fcc_scale, phase=fcc_phase)
        data = reshape(data, (-1, npe, nro))  # add singleton dimension for simplicity of following code
        nchan = 1
    else:
        data = transpose(data, [2,0,1])
    
    if apply_phase_ramp:
        print 'appling one-pixel phase ramp'
        v1 = zeros(nro)
        v2 = zeros(nro)
        v1[0] = 1
        v2[1] = 1
        w = angle(fft(v1) / fft(v2))
        ramp = exp(-1j*w)  # equivalent to 1 pixel spatial shift
        for c in range(nchan):
            data[c,:,:] = data[c,:,:] * ramp
    
    
    if kspace_shift_correct:
        print 'correcting for k-space shifts'
        data[0,:,:] = do_kspace_shift_correct(data[0,:,:], golden_angle=golden_angle)
    
    if zero_dc_phase:
        print 'zeroing DC phase'
        data[0,:,:] = do_zero_dc_phase(data[0,:,:])
    
        
    # CROP RO direction, for easier development
    if crop_ro_oversampling:
        print 'cropping readout oversampling before gridding'
        data = fftshift(ifft(fftshift(data, axes=2), axis=2), axes=2)
        data = data[:,:,nro/4:-nro/4,]
        data = fftshift(fft(fftshift(data, axes=2), axis=2), axes=2).copy('C')
        nro = data.shape[2]
    
    # SET UP K-SPACE COORDINATES
    dz_per_pe = TR * vtable # inter-TR table motion
    zmax = dz_per_pe * npe  # total sc3an z distance in cm, Saikat's 'virtual fov'
    z = linspace(0, zmax, npe, endpoint=False)  # relative table positions at each TR
    
    #kr = linspace(-nro/2, nro/2, nro, endpoint=False)
    kr = linspace(-nro/2, nro/2, nro, endpoint=True)
    if golden_angle:
        theta = arange(npe) * 111.246 * pi / 180.0  # golden angle
    else:
      #  theta = (arange(npe) * pi/ (radial_pct*256 /100.0)) % (2*pi) # linear
         theta = (arange(npe) * pi/ (radial_pct*256 /100.0))  # linear
        
        
    K = zeros((npe, nro, 2))
    kx = outer(cos(theta), kr)
    ky = outer(sin(theta), kr)
    K[:,:,0] = kx
    K[:,:,1] = ky   # don't care about kz for slicewise recon
    
    if readout_alternation:
        print 'flipping every second profile'
        K[1::2,:,:] = -K[1::2,:,:]
    
    
    nx = nro // 4 * 4  # ensure is multiple of 4, for easier cropping
    ny = nx
    n = nx
    ncrop = (osfactor - 1)*nx / 2
    
    # viewing parameters
    grid_shape = array([nx, ny])
    data_shape = array([npe*nro, nchan])
    
    print 'grid shape is', grid_shape
    print 'number of slices:', nslices
    print 'in-plane voxel size:', voxel_size_x_mm
    print 'phase encodes per slice:', npes_per_slice
    print 'padding will be', ncrop
    print 'effective slice thickness:', zmax/nslices
    print 'estimated scan time:', TR*npe

    angle_sweeps = npe* 180.0 / array([256,185,148,92]) 
    print 'angle sweeps:',angle_sweeps
    
    # set up Brian's slice/profile divisions
    slice_idx = 0
    slice_center_z_mm_array = [0]
    slice_center_projection_array = [128]
    while True:
        slice_idx = slice_idx + 1
        next_z_mm = slice_center_z_mm_array[slice_idx-1] + voxel_size_x_mm
        next_projection_idx = int(round( next_z_mm/(vtable * TR) ) + slice_center_projection_array[0])
        if (next_projection_idx+256/2) <= npe:
            slice_center_z_mm_array.append(next_z_mm)
            slice_center_projection_array.append(next_projection_idx)
        else:
            break
    
    nSlices = len(slice_center_projection_array)
    
    im_gridding = zeros((nslices, nx, ny), 'complex64')
    
    t = time.time()
    
    for s in range(nslices): 
            
        pe_center = int(round(voxel_size_x_mm * s / dz_per_pe) + pe_of_z0)
        if npes_per_slice % 2 == 0:
            pe_start = pe_center - npes_per_slice / 2
            pe_end = pe_center + npes_per_slice / 2
        else:
            pe_start = pe_center - (npes_per_slice + 1) / 2
            pe_end = pe_center + (npes_per_slice - 1) / 2
        
        slice_pe_indices = range(pe_start, pe_end)
        
        
        if pe_end > npe-256:
       	 	break
        
        slice_data = data[:, slice_pe_indices, :]    # X = nro x nchan
        if not golden_angle:
            slice_data = interpolate_radial_data(data, slice_data, slice_pe_indices, npe_per_180)
            
        npe_slice = slice_data.shape[1]
        slice_data_shape = (npe_slice*nro, nchan)
        slice_K = K[slice_pe_indices, :, :]
        if s % (nslices/50) == 0:
            print 'processing slice %d/%d: profiles %d -> %d' % (s, nslices, pe_start, pe_end-1)
    
        # compute gridding operator
        slice_K = reshape(slice_K, (npe_slice*nro, 2))
        G = sparse_matrix_operator(slice_K, grid_shape, osfactor=osfactor, verbose=False)
        # compute sample density correction
        W = sample_density(G, niter=5)
    
        # compute rolloff kernel
        spatial_dims = grid_shape
        RO = rolloff(grid_shape, ncrop=ncrop, osfactor=osfactor)
    
        A = lambda X, mode: Aop(X, mode, G, W, RO, ncrop, grid_shape, slice_data_shape)
        slice_data = reshape(slice_data, (-1, nchan))    # [nro*npe_slice, nchan]
        X_gridding = A(slice_data, 2)
        im_gridding[s,:,:] = squeeze(X_gridding)   #sqrt(sum(abs(X_gridding)**2, axis=2))
        
      
    print time.time()-t, 's elapsed'
    
    #mat['im_%d' % npes_per_slice] = im_gridding  # add reconstructed image to existing MAT file and save
    im_file = pathstr_python_script + '/' + os.pardir + '/data_output/' + base_name + '_' + 'PYTHONREC' + '.mat'
    
    # crop axial plane from 512 x 512 to 340 x 340
    # (512 - 340)/2 = 86    
    crop_start = 86
    crop_stop = 86 + 340
    #savemat(im_file, {'img':im_gridding[:,crop_start:crop_stop,crop_start:crop_stop]})
    savemat(im_file, {'img':im_gridding[:,:,:]})
        
    if plotting:
        figure(3)   
        clf()
        x = abs(im_gridding)
        imshow(flipud(x[450,:,:].T), interpolation='nearest', cmap='gray')
        colorbar()
        show()
    
    # end recon_list loop
