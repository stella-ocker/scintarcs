r'''
Useful functions for pulsar scintillation processing.
r'''

from numpy import *

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

def readhdu(data,dat_row,pol): # data = hdu[1].data
    """
    calculate uncalibrated intensities for one row in the data file (one row = 0.1 s)
    inputs = hdu data, data row index, pol (0=aa,1=bb,2=cr,3=ci,-1=total intensity)
    returns data, frequencies, image extent for 2d plotting
    """
    
    #data = fhdu[1].data
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

