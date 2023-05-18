from numpy import *
from matplotlib.pyplot import *

fnums = arange(11,15,1)
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

bins = 2304
tdim = len(all_dat) # number of time samples
#fdim = len(fulldatm[:,0]) # number of frequency samples
m = int(tdim // bins)
print(m,'pulses')
nsamp_eff = m*bins
phase_inds = arange(0,nsamp_eff,bins)
tfold = array([all_dat[k:k+2304] for k in phase_inds])

print(shape(tfold))

amp = 0
z = 0
phase_arr = linspace(0,1,bins)
fig,ax = subplots(figsize=(5,20))
for i in range(m):
    x = phase_arr
    y = tfold[i,:]-amp
    ax.plot(x,y,color='tab:blue',lw=1,zorder=z)
    #ax.fill_between(x,full(len(y),y.min()),y,where=(y>min(y)),color='white',zorder=z)#,edgecolors=('white','tab:blue'),linewidth=1)
    amp+=5
    z+=1
ax.set_yticks([])
ax.set_xlabel('Phase')
savefig('joy_division_B1929_test.pdf')
#show()

