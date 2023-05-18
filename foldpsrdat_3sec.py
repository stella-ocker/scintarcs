r'''
Program to fold filterbank data at pulse period. Loads data, de-disperses, folds, downsamples. Save output for bandpass calibration.
Run for AA and BB separately. Loads one file at a time and splits into two subintegrations (3.2 s = integration time, about 15 pulses).
Command-line usage: python foldpsrdat.py [start file (raw file number)] [end file number (raw file number)] [polnum]
r'''

from numpy import *
from matplotlib.pyplot import *
from astropy.io import fits
import glob,os
import sys
from psrcal_tools import *
from astropy.convolution import Box1DKernel
import astropy.convolution as ap
from scipy.signal import decimate

fstart = int(sys.argv[1])
fend = int(sys.argv[2])
polnum = int(sys.argv[3])

path = 'B1929+10/20211028/B1929+10_tracking-M01_'
flist = arange(fstart,fend+1,1)
flist = [str(flist[i]).zfill(4) for i in range(len(flist))]
data_files = [path+flist[i]+'.fits' for i in range(len(flist))]
print('begin with '+data_files[0]+' end with '+data_files[-1])

finds = arange(0,len(flist),1) # process one file at a time

def form_hdu(hdu,polnum,nrow1,nrow2):
    print('     loading '+hdu)
    hdu1 = fits.open(hdu)
    data = hdu1[1].data
    tsamp = data[0]['TSUBINT']/1024
    nrows = nrow2-nrow1
    full_data = zeros((8192,1024*nrows))
    print('     stitching subintegrations')
    for j in range(nrow1,nrow2,1):
        i2_arr,freqs,extent = readhdu(data,j,polnum)
        full_data[:,(j-nrow1)*1024:(j+1-nrow1)*1024] = i2_arr
    hdu1.close()
    return full_data,freqs,tsamp

def dedisperse(dat,freq,tsamp,dm):
    print('     de-dispersing')
    ref_freq = min(freq) + ((max(freq) - min(freq))/2.) # central frequency to de-disperse to
    dmdelays = dm_delay(dm,freq,ref_freq) # milliseconds
    delay_nums = (dmdelays*10**(-3.))/tsamp
    dedm_arr = array([shiftit(dat[j,:],delay_nums[j]) for j in range(len(dat[:,0]))])
    return dedm_arr

def main(i,nrow1,nrow2):
    
    fulldat,freqs,tsamp = form_hdu(data_files[i],polnum,nrow1,nrow2)
    dedm_arr = dedisperse(fulldat,freqs,tsamp,dm)
    del fulldat

    # cut upper and lower ends of bandpass
    badfreqs1 = where(freqs<1050.)[0][-1]
    badfreqs2 = where(freqs>1450.)[0][0]
    fulldat_cut = dedm_arr[badfreqs1:badfreqs2,:]
    freqcut = freqs[badfreqs1:badfreqs2]

    print('     masking RFI')
    fulldatm,fulldatmask = mask_region(fulldat_cut,freqcut,1145,1210) # change back to dedm_arr if switching order
    fulldatm,fulldatmask = mask_region(fulldatm,freqcut,1138,1142.5,fulldatmask)
    fulldatm,fulldatmask = mask_region(fulldatm,freqcut,1222,1235,fulldatmask)

    #imshow(fulldatm,aspect='auto',origin='lower')
    #title('pre-folding')
    #show()

    #savez('folding_testdat',dat=fulldatm,freqs=freqcut,tsamp=tsamp)
    #break

    # remove unnecessary data from memory
    del dedm_arr,fulldat_cut

    print('     folding data')
    bins = 2304
    tdim = len(fulldatm[0,:]) # number of time samples
    fdim = len(fulldatm[:,0]) # number of frequency samples
    m = int(tdim // bins)
    nsamp_eff = m*bins
    phase_inds = arange(0,nsamp_eff,bins)
    tfold = array([fulldatm[:,k:k+2304] for k in phase_inds])

    del fulldatm

    #tfold = reshape(fulldatm[:,:nsamp_eff],(m,fdim,bins))
    tfold_mask = zeros(shape(tfold))
    minds = where(fulldatmask[:,0]==1)[0]
    tfold_mask[:,minds,:] = 1
    tfold = ma.masked_array(tfold,tfold_mask)

    print('     new data shape',shape(tfold))

    #imshow(tfold[10,:,:],aspect='auto',origin='lower')
    #title('single pulse')
    #show()

    # shift main pulse peak to fixed phase (500 samples)
    tmean = ma.mean(tfold,axis=0) # average over pulses
    pulse_prof = ma.mean(tmean,axis=0) # calculate mean temporal profile to determine phase shift
    pulse_prof2 = ap.convolve(pulse_prof,Box1DKernel(6)) # smooth temporal profile to avoid noise spikes
    mp_center = 500
    mp_window = int(166/2) # half window for main pulse
    shift = int(bins-argmax(pulse_prof2)+mp_center)
    tmean = roll(tmean,shift,axis=1)
    pulse_prof = roll(pulse_prof,shift)
    pulse_prof2 = roll(pulse_prof2,shift)
    ip_center = argmax(pulse_prof2[1000:])+1000-5 # correction for asymmetry 
    ip_w1 = int(55/2) # half window for interpulse
    ip_w2 = int(55/2)
    ip_window = ip_w2+ip_w1

    mppinds = arange(mp_center-mp_window,mp_center+mp_window,1)
    ippinds = arange(ip_center-ip_w1,ip_center+ip_w2,1)
    pinds = append(mppinds,ippinds)

    print('     saving on/off-pulse I(f)')
    mainpulse = ma.mean(tmean[:,mppinds],axis=1)
    interpulse = ma.mean(tmean[:,ippinds],axis=1)
    off_mp = ma.mean(tmean[:,50:50+2*mp_window],axis=1)
    off_ip = ma.mean(tmean[:,50:50+ip_window],axis=1)
    offpulse = delete(tmean,pinds,axis=1)
    bandpass = ma.mean(offpulse,axis=1)
    savez('/home/fast04/processed_data/B1929+10_out/folded/B1929+10_1D_ON_OFF_'+flist[i]+'_'+str(nrow1)+'_POL'+str(polnum),mainpulse=mainpulse,interpulse=interpulse,off_mp=off_mp,
            off_ip = off_ip,bandpass=bandpass,freqs=freqcut)


    #print(int(shape(tmean)[1]/4)-40,shape(datcal_resamp))
    print('     downsampling folded 2D spectrum')
    datcal_resamp = zeros((fdim,int(bins/4)-40))
    datcal_mask = zeros(shape(datcal_resamp))
    for l in range(len(datcal_resamp[:,0])):
        trow = tmean[l,:]
        if mean(fulldatmask[l])==1: # preserve masked frequency channels
            tdec = zeros(len(datcal_resamp[l,:]))
            datcal_mask[l,:] = 1
        else:
            tsmooth = ap.convolve(trow,Box1DKernel(4)) # smoothing before decimating
            tdec = decimate(tsmooth,4)[20:-20] # cuts off windowing
        datcal_resamp[l,:] = tdec
    times = linspace(0,shape(datcal_resamp)[1]*tsamp,int(shape(datcal_resamp)[1]))

    datcal_resamp = ma.masked_array(datcal_resamp,datcal_mask)

    print('     saving folded 2D spectrum')
    datcal_reshape = reshape(swapaxes(datcal_resamp,0,1),(1,shape(datcal_resamp)[0],shape(datcal_resamp)[1]))
    hdr = fits.Header()
    primary_hdu = fits.PrimaryHDU(header=hdr)
    c1 = fits.Column(name='FREQS',format=str(len(freqcut))+'E',unit='MHz',array=reshape(freqcut,(1,len(freqcut))))
    c2 = fits.Column(name='TIMES',format=str(len(times))+'E',unit='sec',array=reshape(times,(1,len(times))))
    c3 = fits.Column(name='DATA',format=str(len(datcal_resamp.flatten()))+'B',dim=str(shape(datcal_reshape)),unit='Counts',array=datcal_reshape)
    hdu2 = fits.BinTableHDU.from_columns([c1,c2,c3])
    hdu_new = fits.HDUList([primary_hdu,hdu2])
    hdu_new.writeto('/home/fast04/processed_data/B1929+10_out/folded/B1929+10_'+flist[i]+'_'+str(nrow1)+'_POL'+str(polnum)+'.fits',overwrite=True)

    del datcal_reshape,datcal_resamp,tfold,tmean

    return 

dm = 3.18321
for i in finds:
    
    # splitting each file into two subintegrations
    main(i,0,32)
    main(i,32,64)



