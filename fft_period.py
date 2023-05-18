r'''
Program to calculate period from FFT of time series.
r'''

from numpy import *
from matplotlib.pyplot import *
import glob,os
from scipy.signal import find_peaks

fnums = arange(11,411,1)
data_files = ['/home/griz03/stella/B1929+10_out/single_pulses/fft_out/time_series'+str(fnums[i]).zfill(4)+'.npz' for i in range(len(fnums))]

# load and string together time series
print('loading data files')
all_dat = []
for i in range(len(data_files)):
    #print('     loop',i)
    df = load(data_files[i])
    tdat = df['tdat']
    tsamp = df['tsamp']
    all_dat.append(tdat)

all_dat = array(all_dat).flatten()

#plot(all_dat)
#show()

print('calculating FFT')
tspec_fft = abs(fft.fft(all_dat))**2.
psfreqs = fft.fftfreq(all_dat.size,tsamp)
psdcut = tspec_fft[:int(len(psfreqs)/2)]
psfcut = psfreqs[:int(len(psfreqs)/2)]

peaks,pprops = find_peaks(psdcut,distance=10000)
f0 = psfcut[peaks[50]]/50.
p0 = 1./f0
#print(shape(peaks))
print('pulse frequency (Hz)',psfcut[peaks[50]]/50.)
print('FFT frequency resolution (Hz)',psfcut[1]-psfcut[0])
print('pulse period (s)',p0)
print('number of samples per pulse phase',p0/tsamp)

plot(psfcut,psdcut)
plot(psfcut[peaks],psdcut[peaks],'o')
#yscale('log')
xlim(0,150)
ylim(0,1e16)
xlabel('Frequency (Hz)')
ylabel('PSD')
show()

#savez('B1929_periodogram',psd=tspec_fft,psfreqs=psfreqs)
