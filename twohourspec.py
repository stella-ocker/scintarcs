from numpy import *
from matplotlib.pyplot import *
from astropy.io import fits
from astropy.convolution import convolve,Box1DKernel
from glob import glob
from scipy.signal import decimate


def readhdu(fhdu,dat_row):
    """
    calculate total intensity (uncalibrated) for one row in the data file (one row = 0.1 s)
    inputs = hdu, data row index, DM
    returns data, frequencies, image extent for 2d plotting
    """

    data = fhdu[1].data
    freqs = data[dat_row]['DAT_FREQ']
    tsamp = data[dat_row]['TSUBINT']/1024 # 98 microseconds -- sampling time
    tsubint = data[dat_row]['TSUBINT'] # subint time = 0.1 s
    data_arr = squeeze(data[dat_row]['DATA'])
    aabb = data_arr[:,0:2,:]
    totalI2 = swapaxes(sum(aabb,axis=1),1,0)

    image_extent = [0,tsubint*1000,freqs[0],freqs[-1]] # time in ms, freqs in MHz

    return totalI2,freqs,image_extent

path = 'B1929+10/20211028/'
f = sort(glob(path+'*.fits'))

# first and last 10 files (with noise injection)
f1 = concatenate((f[:10],f[-10:])).flatten()

# files without noise injection
f2 = f[10:-10]

# 1098 files w/o noise injection -- divisible by factors of 3

# run through data w/o noise injection, calculate dynamic spectrum
# (64 0.1 s subintegrations per file)
alldat = []
for n,fi in enumerate(f2):
    print(fi)
    hdu = fits.open(fi)
    hdr = hdu[1].header
    nrows = hdr['NAXIS2']
    print('nrows:',nrows,'dim:',hdr['TDIM17'])
    full_data = zeros((8192,1024*nrows))
    for j in range(nrows):
        i2_arr,freqs,extent = readhdu(hdu,j)
        full_data[:,j*1024:(j+1)*1024] = i2_arr
    hdu.close()
    freq1d = mean(full_data,axis=1)
    #freq1d_down = decimate(freq1d,2) # downsample to 4k frequency channels
    alldat.append(freq1d)
alldat = array(alldat)
alldat = swapaxes(alldat,0,1) # make frequency first axis

print(shape(alldat))
savez('twohourspec',dat=alldat,freqs=freqs,time=arange(0,len(f2)*6.4,len(f2)))

r'''
# average every three spectra
alldat_mean = []
for i in range(len(alldat_mean)/3):
    alldat_mean.append(mean(alldat[i:i+3],axis=0)
#savez('twohourspec',dat=alldat_mean,freqs=freqs,time=arange(0,1098*6.4,3*6.4))

print(shape(alldat_mean))
r'''

imshow(alldat,origin='lower',aspect='auto',extent=[0,len(f2)*6.4,freqs[0],freqs[-1]])
xlabel('Time (s)')
ylabel('Frequency (MHz)')
title('Total Intensity (38 s Resolution)')
savefig('twohourspec.pdf')
    



