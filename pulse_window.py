r'''
Program to calculate optimal pulse windows for main pulse and interpulse (assuming mean profile for 57 pulses, 12.8 s of integration time).
r'''

from numpy import *
from matplotlib.pyplot import *
import glob,os
from scipy.signal import find_peaks

fnums = arange(11,211,1)
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

peaks,pprops = find_peaks(psdcut,distance=5000)
f0 = psfcut[peaks[50]]/50.
p0 = 1./f0
#print(shape(peaks))
print('pulse frequency (Hz)',psfcut[peaks[50]]/50.)
print('FFT frequency resolution (Hz)',psfcut[1]-psfcut[0])
print('pulse period (s)',p0)
print('number of samples per pulse phase',p0/tsamp)

print('folding time series')
bins = int(p0/tsamp)
m = int(len(all_dat) // bins)
nsamp_eff = m*bins
tfold = reshape(all_dat[:nsamp_eff],(m,bins))
print(shape(tfold))

#plot(tfold[0])
#plot(tfold[25])
#plot(tfold[57])
#show()

# joy division plot for 57 phase bins ~ equivalent to 12.8 s (2 files)
#amp = 0
#for i in range(57):
#    plot(tfold[i]+amp,color='k',lw=1)
#    amp+=5
#yticks([])
#show()

tmean = mean(tfold[:57,:],axis=0)
tmean = tmean - min(tmean)
thresh = 0.05*max(tmean) # experimenting with 90% intensity window
w1 = where(tmean>thresh)[0][0]
w2 = where(tmean>thresh)[0][-1]
pwidth = w2-w1
print('main pulse width in samples',pwidth)

ip_cut = tmean[:1000]
ip_thresh = 0.25*max(ip_cut)
ip_w1 = where(ip_cut>ip_thresh)[0][0]
ip_w2 = where(ip_cut>ip_thresh)[0][-1]
ip_width = ip_w2 - ip_w1
print('interpulse width in samples',ip_width)

ip_mp_dist = argmax(tmean) - argmax(tmean[:1000])
print('distance between interpulse and main pulse in samples',ip_mp_dist)

plot(tmean)
plot(w1,tmean[w1],'o')
plot(w2,tmean[w2],'o')
plot(ip_w1,ip_cut[ip_w1],'o')
plot(ip_w2,ip_cut[ip_w2],'o')
xlabel('Samples')
ylabel('Mean pulse profile (50 pulses)')
show()



