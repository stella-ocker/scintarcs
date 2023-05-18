r''' 
Program to make 1D time series at 98 us (for period determination).
Command-line usage: time_series.py [file start number (in raw file units)] [file end number (in raw file units)]
r'''

from numpy import *
from matplotlib.pyplot import *
from psrcal_tools import *
from astropy.io import fits
import glob,os
import sys

path = '/home/fast02/B1957+20/20211019/'
fstart = int(sys.argv[1])
fend = int(sys.argv[2])
fnums = arange(fstart,fend,1)
data_files = [path+'B1957+20_tracking-M01_'+str(fnums[i]).zfill(4)+'.fits' for i in range(len(fnums))]

dm = 29.1066
for i in range(len(data_files)):
    print('loading '+data_files[i])
    hdu = fits.open(data_files[i])
    data = hdu[1].data
    nrows = len(data)
    tsamp = data[0]['TSUBINT']/1024
    full_data = zeros((8192,1024*nrows))
    print('     stitching subintegrations')
    for j in range(nrows):
        #print('row ',j)
        i2_arr,freqs,extent = readhdu(hdu,j,-1)
        full_data[:,j*1024:(j+1)*1024] = i2_arr # not dedispersed, no RFI mask
    hdu.close()

    # cut upper and lower ends of bandpass
    badfreqs1 = where(freqs<1050.)[0][-1]
    badfreqs2 = where(freqs>1450.)[0][0]
    fulldat_cut = full_data[badfreqs1:badfreqs2,:]
    freqcut = freqs[badfreqs1:badfreqs2] 

    fulldat_cut,fulldatmask = mask_region(fulldat_cut,freqcut,1138,1300)

    # de-disperse
    print('     de-dispersing')
    ref_freq = min(freqs) + ((max(freqs) - min(freqs))/2.) # central frequency to de-disperse to
    dmdelays = dm_delay(dm,freqs,ref_freq) # milliseconds
    delay_nums = (dmdelays*10**(-3.))/tsamp
    dedm_arr = array([shiftit(fulldat_cut[j,:],delay_nums[j]) for j in range(len(fulldat_cut[:,0]))])

    tdat = mean(fulldat_cut[freqcut>1250.,:],axis=0)
    savez('/home/fast04/processed_data/B1957+20_out/single_pulses/fft_out/time_series'+str(fnums[i]).zfill(4),tdat=tdat,tsamp=tsamp)

    #plot(tdat)
    #show()
    #break
