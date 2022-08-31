#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 2 16:53:00 2022

@author: mozef

This code does:
- Make a plot of correlogram gather
- Make FK-transform over several cross-correlograms/traces
- Pick highest amplitude manually if needed
- If highest amplitude is picked, will plot phase velocity

"""

#export HDF5_USE_FILE_LOCKING=FALSE
import numpy as np
import matplotlib.pyplot as plt
from xdas.io.febus import read as read_das
from icoords import InterpolatedDataArray
from xdas.io.febus import correct_gps_time
from scipy import signal


#========================================================== Functions =========================================================

# from Alister
def fk_analysis(x, dt, ds, dbscale=True, onesided=True):                        
    """Compute the FK analysis of x with expected dims (time, offset)"""        
    fk = np.fft.fft2(x)                                                         
    f = np.fft.fftshift(np.fft.fftfreq(fk.shape[0], dt))                        
    k = np.fft.fftshift(np.fft.fftfreq(fk.shape[1], ds))                        
    fk = np.fft.fftshift(fk)                                                    
    if onesided:                                                                
        mask = (f >= 0)                                                         
        fk = fk[mask]                                                           
        f = f[mask]                                                             
    if dbscale==True:                                                                 
        fk = 20 * np.log(np.abs(fk)) 
    else:
        fk = abs(fk)                                     
    return f, k, fk  

def filter_bandpass(trace, lowcut, highcut, fs, order=2):
	nyq = 0.5 * fs
	low = lowcut/nyq
	high = highcut/nyq
	order = 2
	b,a = signal.butter(order, [low, high], 'bandpass', analog=False)
	trace_filtered = signal.filtfilt(b, a, trace, axis=0)
	return trace_filtered

def normalize(trace):	# normalize trace
		return trace / max(trace)

def stack(folders):
	stacked = np.zeros(len(read(folders[0])[0]))
	for file in folders:
		a = read(file)
		stacked = stacked + a[0].data
	return stacked

# gather the tf-pws
def collect_tfpws(folder,x1,x2):
	files = glob.glob(f'{folder}/Correlogram_{str(x1)}_{str(x2)}_phase/ts_pws_tspws_pcc.sac') # only consider the phase
	#files = glob.glob(f'{folder}/Correlogram_{str(x1)}_{str(x2)}_phase/ts_pws_tspws_pcc_half_zeromid.sac') # only consider the phase
	print(files)
	try:
		stacked = read(files[0])[0].data
	except:
		files = glob.glob(f'{folder}/Correlogram_{str(x1)}_{str(x2-1)}_phase/ts_pws_tspws_pcc_half_zeromid.sac') # try except works when using half stack
		stacked = read(files[0])[0].data
	return stacked


#=================================================== Correlogram gather plot ====================================================
from obspy import read
import glob

#folder = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/01-10Hz_Window300_mid100_gl10'
folder = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/Auto-correlation/test/'
# to plot half stack, change: collect_tfpws, file, and plt.plot, and cc_fk

# parameters
start = 100
channel = 200  # define the number of channel needed to be gathered, depend on start_position. not used when start_position = 'mid_all'
	       # because automatically use all channels
window = 300
start_position = 'mid_all'
dt = 0.02
fs = 1 / dt

x1, x2 = start, 100
file = f'{folder}/Correlogram_{str(x1)}_{str(x2)}_phase/ts_pws_tspws_pcc.sac'
#file = f'{folder}/Correlogram_{str(x1)}_{str(x2)}_phase/ts_pws_tspws_pcc_half_zeromid.sac' 
a = read(file)[0].data

cut1 = 100   # start and end of the channel to be plotted
cut2 = 102
if start_position == 'mid':
	cc = np.zeros((a.shape[0],2*channel+1))
elif start_position == 'mid_all':
	cc = np.zeros((a.shape[0],cut2-cut1+1))	
else:
	cc = np.zeros((a.shape[0],channel+1))

t = np.arange(0,a.shape[0]*dt,dt)

if start_position == 'mid':
	ranges = np.arange(start-channel,start+channel+1,1)
elif start_position == 'mid_all':
	ranges = np.arange(cut1,cut2+1,1)
elif start_position == 'early':
	ranges = np.arange(start,start+channel+1,1)
elif start_position == 'end':
	ranges = np.arange(start,start-channel-1,-1)

lowcut = 0.1  # filter if needed
highcut = 10
j = 0

# collect all tf-pws of each channel pair and store it in an array
for i in ranges:
	x1,x2 = start, i
	#folders = glob.glob(f'{folder}/Correlogram_{str(x1)}_{str(x2)}_{mode}/*.sac')
	#stacked = stack(folders)
	stacked = collect_tfpws(folder,x1,x2)
	#stacked_filtered = filter_bandpass(stacked, lowcut, highcut, fs, order=2)
	cc[:,j] = normalize(stacked)
	j = j + 1

# set time boundary
time_lim1 = -10
time_lim2 = 10
boundary1 = int(((time_lim1 + window) / (window * 2)) * cc.shape[0])
boundary2 = int(((time_lim2 + window) / (window * 2)) * cc.shape[0])
#boundary1 = int(((time_lim1) / (window)) * cc.shape[0]) # use this for the half stack
#boundary2 = int(((time_lim2) / (window)) * cc.shape[0])
t_cut = np.arange(time_lim1,time_lim2,dt)

lag = len(cc[:,0])/2 - np.where(cc[:,0]==max(cc[:,0]))[0][0]
t_cut_zeromid = np.arange(time_lim1+lag*dt,time_lim2+lag*dt,dt)
cc_cut = cc[boundary1:boundary2,:]

# functions for plotting secondary axis (offset) in plot
def get_channel(x):
	return x

def get_offset(x):
	return 2.4 * (x)

scaler = 1
remove= []

plt.figure(figsize=(7,12))
j = 0
for i in ranges:
	if i not in remove:
		#plt.plot(t, normalize(cc[:,j]) + i*1, color='black',alpha=0.6)
		plt.plot(t_cut,cc_cut[:,j] + (i)*scaler, color='black',alpha=0.6) # use t_cut_zeromid if needed
		'''try:
			plt.plot(t_cut,normalize(cc_cut[:,j]) + (i)*scaler, color='black',alpha=0.6) 
		except:
			plt.plot(t_cut,normalize(np.append(cc_cut[:,j],0)) + (i)*scaler, color='black',alpha=0.6) 
		'''
	j = j + 1

ax = plt.gca()
y2 = ax.secondary_yaxis('right', functions=(get_offset, get_channel))
y2.set_ylabel('Offset (m)')
plt.title('Filter 0.1-10 Hz')
plt.xlabel('Lag Time (s)')
plt.ylabel('Channel')
#plt.xlim(0,5)
plt.show()


# np.savetxt("west_256_430.txt", cc_fk.T,fmt="%s")
# content = np.genfromtxt('west_256_430.txt', delimiter=',', dtype=None,encoding = 'UTF-8')

'''
# for plotting cross-correlogram gather without certain channels
#remove = [50,51,52,53,54,55,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250]
remove= []
mask = np.isin(ranges,remove)
mask = [not mask[i] for i in range(len(mask))]
cc_cut = cc[boundary1:boundary2,ranges[mask]-50]
'''

# ========================================================== F-K computation =========================================================
dx = 2.4
dt= 0.02
print('dx=',dx,' dt=',dt,'  data= ',str(2020))

pick = False    # change to True if pick is needed
cc_fk = cc_cut  # use 'cc' if half stack is used
#cc_fk = cc
f, k, fk = fk_analysis(cc_fk, dt, dx, dbscale=True, onesided=True)    

fig = plt.figure(figsize=(8,8))                                           
ax = plt.gca()                                                    
fkmax = np.nanmax(fk.flatten())                                             
fkmin = np.nanmin(fk.flatten())                                       
im = ax.pcolor(k,f,fk,vmin = -100,vmax=fkmax,cmap=plt.cm.viridis)
ax.set_xlabel('Wavenumber (1/m)')                                                                 
ax.set_ylabel('Frequency (Hz)')                                             
ax.set_title(f'DAS Stromboli 2020 \nF-K Transform of Correlogram between Channel {str(ranges[0])} and {str(ranges[-1])}' ) 
colorbar = plt.colorbar(im, aspect=50)
colorbar.set_label(label='Amplitude (dB)', rotation=90)
plt.tight_layout()
plt.ylim(0.1,10)
plt.xlim(-0.05,0.05)           

# pick value
if pick:
	print("Pick the strong amplitude. Left button to pick, right button to cancel last pick, middle button to stop (max=100 picks)")
	x = plt.ginput(100)
	print(x)
	k_pick = np.array([i[0] for i in x])
	f_pick = np.array([i[1] for i in x])
	plt.plot(k_pick,f_pick,'--r',lw=1)
	plt.plot(-k_pick,f_pick,'--b',lw=1)

plt.show()


# ==================================================== phase velocity plotting ======================================================
v = np.array(f_pick) / np.array(k_pick)
period = 1 / np.array(f_pick)

# plotting in frequency
plt.figure()
plt.plot(f_pick,abs(v))
plt.title(f'DAS Stromboli 2020 \nPhase velocity between Channel {str(ranges[0])} and {str(ranges[-1])}' ) 
#plt.xlabel('Frequency (Hz)')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Velocity (m/s)')
plt.show()

# plotting in period
plt.figure()
plt.plot(period,abs(v))
plt.title(f'DAS Stromboli 2020 \nPhase velocity between Channel {str(ranges[0])} and {str(ranges[-1])}' ) 
#plt.xlabel('Frequency (Hz)')
plt.xlabel('Period (s)')
plt.ylabel('Velocity (m/s)')
plt.show()


'''
# Auto-pick fk plot
# still hard to use because of the maximum amplitude of FK plot is highly scattered (maybe smoothing is needed?)

from scipy.ndimage import gaussian_filter
from scipy.interpolate import make_interp_spline, BSpline
#fk_int = gaussian_filter(fk, sigma=2)

cut = 0
f_peak = [] # getting the index of v based on the highest stransform value 
f_low = []
f_high = []

for idx_k in range(fk.T.shape[0]):
	#if k[idx_k] < 0:
	print('yes')
	a = np.where(fk[cut:,idx_k] == np.max(fk[cut:,idx_k]))
	bound = np.where((fk[cut:,idx_k] > 0.3 * np.max(fk[cut:,idx_k])) & (fk[cut:,idx_k] < np.max(fk[cut:,idx_k])))
	b = f[a[0][0]+cut]
	bound_high = f[bound[0][-1]+cut]
	bound_low = f[bound[0][0]+cut]
	f_peak.append(b)
	f_low.append(bound_low)
	f_high.append(bound_high)

def interpolate(f,v):
	f_int = np.linspace(f[0], f[-1], 200)
	interpolate = make_interp_spline(sorted(f), v, k=3)
	v_int = interpolate(f_int)
	return f_int, v_int

k_int, f_peak = interpolate(k,np.array(f_peak))
k_int, f_low = interpolate(k,np.array(f_low))
k_int, f_high = interpolate(k,np.array(f_high))
'''



#=================================================== FK plot using traces data =========================================================

iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/SR_Stromboli_2020-09-22_19-34-16_UTC_decimated50Hz.nc") # many small events and one big event
Das_nc = iDas_nc.load_icoords()

dt = float(Das_nc.time[1]-Das_nc.time[0])/1.e9
time = np.arange(0,Das_nc.shape[0]*dt,dt)
starttime = Das_nc.time.data[0].astype('datetime64[ms]')
t1 = 2280
t2 = 2340
time_cut = np.arange(t1,t2,dt)

lowcut = 5
highcut = 10
fs = 1/dt
apply_filter = True

if apply_filter:
	for i in range(Das_nc.shape[1]):
		Das_nc[:,i] = filter_bandpass(Das_nc[:,i], lowcut, highcut, fs, order=2)


start = 0
channel = 463 - start
das = np.flip(Das_nc[int(t1/dt):int(t2/dt),start:start+channel].T,0) # the purpose of this weird line is to rearrange the 2D array for imshow so that the starting point of offset and time is in the bottom left
starttime = Das_nc.time.data[0].astype('datetime64[ms]')
extent = [time[int(t1/dt)],time[int(t2/dt)],start,start+channel]

plt.figure(figsize=(10,8))
a = plt.imshow(das, aspect='auto', vmin=np.min(das.data)/10, vmax=np.max(das.data)/10, extent=extent, cmap='seismic')
colorbar = plt.colorbar(a,orientation="vertical", pad=0.05, location='right', aspect=100)
colorbar.set_label(label='Amplitude ($10^{-8}$ strain $s^{-1}$)', rotation=90)
plt.xlabel('Time (s)')
plt.ylabel('Channel')
#plt.title(f'Start time: {starttime} \nTime Window: {str(t1)} to {str(t2)} s' ) 
plt.title('5-10 Hz')
plt.clim(-1600,2200)
plt.show()

# F-K computation    
dx = float(Das_nc.offset[1] - Das_nc.offset[0])
dt= float(Das_nc.time[1]-Das_nc.time[0])/1.e9

print('dx=',dx,' dt=',dt,'  start= ',str(time_cut[0]))

f, k, fk = fk_analysis(das.T, dt, dx, dbscale=True, onesided=True)    

plt.figure(figsize=(8,8))                                                   
ax = plt.gca()                                                                          
fkmax = np.nanmax(fk.flatten())                                             
fkmin = np.nanmin(fk.flatten())                                       
im = ax.pcolor(k,f,fk,vmin = fkmin,vmax=fkmax,cmap=plt.cm.viridis)                                        
#ax.set_xlim(-0.5,0.5)     
#ax.set_ylim(0.05,0.5)      
ax.set_xlabel('Wavenumber (1/m)')                                              
ax.set_ylabel('Frequency (Hz)')                                             
ax.set_title(f'Start time: {starttime} \nTime Window: {str(t1)} to {str(t2)} s' ) 
colorbar = plt.colorbar(im, aspect=50)
colorbar.set_label(label='Amplitude (dB)', rotation=90)
plt.tight_layout()
plt.show()



