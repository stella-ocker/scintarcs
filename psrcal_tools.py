r'''
Useful functions for pulsar scintillation processing.
r'''

from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import interp1d
from astropy.convolution import convolve,Box1DKernel

def mask_region(datarr,freqarr,freq1,freq2,inmask=None): # freqarr should be ordered from low to high freqs
    # masks frequency channels between freq1 and freq2
    if inmask is not None:
        mask = inmask
    else:
        mask = zeros(shape(datarr))
    ind1 = where(freqarr>freq1)[0][0]
    ind2 = where(freqarr<freq2)[0][-1]
    mask[ind1:ind2] = 1
    return ma.masked_array(datarr,mask),mask

def RFImask(data_cut,stdnum):
    # flag RFI
    mask = zeros(shape(data_cut))
    spec_cut = data_cut[:,:]
    stdspec = std(spec_cut)
    meanspec = mean(spec_cut)
    #print(stdspec,meanspec)
    flagrows = []
    for i in range(len(spec_cut[:,0])):
        row = spec_cut[i,:]
        #print(mean(row))
        upthresh = stdnum*stdspec + meanspec
        lowthresh = meanspec - stdnum*stdspec
        if mean(row)>upthresh:
            flagrows.append(i)
        if mean(row)<lowthresh:
            flagrows.append(i)
    flagrows = array(flagrows)
    if len(flagrows)>0:
        mask[flagrows] = 1
        maskspec = ma.masked_array(data_cut,mask=mask)
        return maskspec,mask

    else:
        return data_cut,mask

def RFIzaptime(data_cut,stdnum,pulse_edge):
    # flag RFI
    mask = zeros(shape(data_cut))
    spec_cut = data_cut[:,pulse_edge:] # grab data to right of pulse
    #stdspec = std(spec_cut)
    #meanspec = mean(spec_cut)
    #print(stdspec,meanspec)
    flag_inds = []
    #for i in range(len(spec_cut[:,0])):
    row = mean(spec_cut,axis=0)
    stdspec = std(row)
    meanspec = mean(row)
    #print(mean(row))
    upthresh = stdnum*stdspec + meanspec
    lowthresh = meanspec - stdnum*stdspec
    flags1 = where(row<upthresh)[0] + pulse_edge
    flags2 = where(row>upthresh)[0] + pulse_edge
    flag_inds = array(flag_inds)
    if len(flag_inds)>0:
        mask[:,flag_inds] = 1
        maskspec = ma.masked_array(data_cut,mask=mask)
        return maskspec,mask

    else:
        return data_cut,mask

def readhdu(fhdu,dat_row,pol):
    """
    calculate uncalibrated intensities for one row in the data file (one row = 0.1 s)
    inputs = hdu, data row index, pol (0=aa,1=bb,2=cr,3=ci,-1=total intensity)
    returns data, frequencies, image extent for 2d plotting
    """
    
    data = fhdu[1].data
    freqs = data[dat_row]['DAT_FREQ']
    tsamp = data[dat_row]['TSUBINT']/1024 # 98 microseconds -- sampling time
    tsubint = data[dat_row]['TSUBINT'] # subint time = 0.1 s
    data_arr = squeeze(data[dat_row]['DATA'])
    I2 = swapaxes(data_arr[:,pol,:],1,0)
    if pol==-1: # total intensity, otherwise index is pol
        aabb = data_arr[:,0:2,:]
        I2 = swapaxes(sum(aabb,axis=1),1,0)
    image_extent = [0,tsubint*1000,freqs[0],freqs[-1]] # time in ms, freqs in MHz
    
    return I2,freqs,image_extent

def dm_delay(dm,nu1,nu2):
    """
    calculate time delay for given DM and 2 frequencies
    """
    return 4.15*10**(6) * dm * ((1./nu1)**2. - (1./nu2)**2.) # ms

def shiftit(y, shift):

    """
    for de-dispersion (code originally from Jim)
    shifts array y by amount shift (in sample numbers)
    uses shift theorem and FFT
    shift > 0  ==>  lower sample number (earlier)
    modeled after fortran routine shiftit
    """

    yfft = fft.fft(y)
    constant = (shift*2*pi)/float(size(y))
    theta = constant*arange(size(y))
    c = cos(theta)
    s = sin(theta)
    work = zeros(size(y), dtype='complex')
    work.real = c * yfft.real - s * yfft.imag
    work.imag = c * yfft.imag + s * yfft.real


    # enforce hermiticity:
    nhalf = int(size(y)/2)
    nhalf_arr = arange(nhalf)
    inds1 = size(y)-nhalf_arr-1
    inds2 = nhalf_arr+1
    work.real[inds1] = work.real[inds2]
    work.imag[inds1] = -work.imag[inds2]
    work[nhalf] = 0.+0.j
    workifft = fft.ifft(work)

    return workifft.real

def gaussian(x, mu, sig, A):
    return A*exp(-power(x - mu, 2.) / (2 * power(sig, 2.)))

def calc_secspec(dss,dT,dF,nT,nF):
    ss = abs(fft.fftshift(fft.fft2(dss)))**2.
    find = int(shape(ss)[0]/2)
    sshalf = ss[find:,:]
    ssconjT = fft.fftshift(fft.fftfreq(nT, d=dT))
    ssconjF = fft.fftshift(fft.fftfreq(nF, d=dF))
    ssconjF_half = ssconjF[-find-1:-1] # to get zero point right
    return sshalf,ssconjT,ssconjF_half

def parabola(eta,y):
    return sqrt(y/eta) # operate on indices

def parabfit(secspec,curv_min,curv_max,dt,dlambda,min_fd,max_fd):
    '''
    inputs: secondary spectrum (either pos. or neg. ft half), in log10 units
            min, max curvatures to sum over
            sampling resolution in ft and flambda
            minimum index in flambda to start summing from
            maximum index in flambda to start summing from
    '''
    eta_arr = logspace(log10(curv_min),log10(curv_max),1000)
    eta_arr_datunits = eta_arr*dlambda/((dt*1000)**2)
    y = arange(min_fd,max_fd,1)
    parab_arr = array([parabola(eta_arr[i],y) for i in range(len(eta_arr))],dtype=int)  
    ft_len = len(secspec[0,:])
    sums = []
    for n,pi in enumerate(parab_arr):
        pi = delete(pi,where(pi>=ft_len)) # get rid of delay values larger than data extent
        slic = array([secspec[y[i],pi[i]] for i in range(len(pi))])
        power = ma.mean(slic) # do mean on log10 if ssh1 != caltest
        sums.append(power)
    sums = array(sums)
    return eta_arr_datunits,sums

'''
def dist_eff(x,vp,dist):
    
    #input eta in m^-1 mHz^-2, vp in km/s, dist in kpc
    
    eta2 = x*(1e6) # convert from mHz^-2 to s^2
    vp2 = vp*1000. # convert from km/s to m/s
    dist2 = dist*(3.086e+19) # convert from kpc to m
    fac1 = 2.*eta2*(vp2**2)/dist2
    s = 1./(fac1**(-1) + 1)
    #deff2 = deff/(3.086e+16) # convert from m to pc
    return s

def eta_from_s(x,vp,dist):
    dist2 = dist*(3.086e+19) # m
    vp2 = vp*1000 # m/s
    denom = 2*((1-x)*(vp))**2
    eta = dist*x*(1-x)/denom
    eta2 = eta/1e6 # convert s^2 to mHz^-2
    return eta2
'''


def fill_dynspec(freqs_up,freqs_low,dss_low,obs_time):
	'''
	Fill dynamic spectrum in lower frequency band with zeros to match size of upper band.
	'''
	fsamp = freqs_up[1]-freqs_up[0]
	lowfill = zeros((len(freqs_up)-len(freqs_low),len(obs_time)))
	dsslowfill = concatenate((dss_low,lowfill))
	freqs = arange(min(freqs_low),min(freqs_low) + len(freqs_up)*fsamp,fsamp)
	return dsslowfill,freqs


def secspec_lambdagrid(dss,freqs_up,freqs_low,obs_time,lowband=False,lowband_fill=False,cutnoise=False):

	c = 2.998e+8 # m/s

	if lowband:
		freqs = copy(freqs_low)
	else:
		freqs = copy(freqs_up)

	# applying hanning window before resampling works better
	window2d = sqrt(outer(hanning(len(freqs)),hanning(len(obs_time))))
	dss = dss*window2d

	if lowband_fill:
		dss,freqs = fill_dynspec(freqs_up,freqs_low,dss,obs_time)

	lambdas = c/(freqs*1e6) # freqs_up originally MHz
	lambdas_grid = linspace(max(lambdas),min(lambdas),len(freqs)) # equi-spaced wavelength grid
	ds_f = interp1d(lambdas,dss,kind='zero',axis=0)

	ds_regrid = ds_f(lambdas_grid)

	dt = abs(obs_time[1]-obs_time[0])
	dl = abs(lambdas_grid[1]-lambdas_grid[0])
	nt = len(obs_time)
	nl = len(lambdas_grid)

	ss,ft,fd = calc_secspec(ds_regrid,dt,dl,nt,nl)
	ss = ss[1:,:]

	ss_offarc1 = log10(mean(ss[:,:50],axis=1))
	ss_offarc2 = log10(mean(ss[:,-50:],axis=1))
	ss_offarc = mean((ss_offarc1,ss_offarc2),axis=0)

	ss_offarc_smooth = convolve(ss_offarc,Box1DKernel(6),boundary='extend')
	ss_offarc_smooth[0] = ss_offarc[0]

	mean_offarc = mean(ss_offarc_smooth[500:])
	print(mean_offarc)

	if cutnoise:
		# remove noise
		caltest = zeros(shape(ss))
		for i in range(len(ss[0,:])):
			col = log10(ss[:,i])
			col -= ss_offarc_smooth
			caltest[:,i] = col + mean_offarc # to compare directly with pre-subtraction
		return caltest,ft,fd

	else:
		return ss,ft,fd




