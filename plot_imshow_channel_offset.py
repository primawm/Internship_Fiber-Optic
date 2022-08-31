#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:09:00 2022

@author: mozef
"""
#export HDF5_USE_FILE_LOCKING=FALSE
from xdas.io.febus import read as read_das
from icoords import InterpolatedDataArray
import matplotlib.pylab as plt
import numpy as np
from xdas.io.febus import correct_gps_time
from scipy import signal
import pandas as pd
from stack_functions import *

iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/DAS_2020_200Hz_nc/SR_Stromboli_2020-09-22_19-34-16_UTC_decimated200Hz.nc") # many small events and one big event
#iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/SR_Stromboli_2020-09-22_19-34-16_UTC_decimated50Hz.nc") # many small events and one big event
Das_nc = iDas_nc.load_icoords()

plt.figure()
plt.plot(Das_nc[:,100])
plt.show()


#a = pd.read_csv('/home/mozef/Documents/station/Geolocation_DAS/2020/DAS_channels_2020_UTM.csv')  # offset information
#offset = a['offset']
#channel = a['channel']

def get_channel(x):
	return x

def get_offset(x):
	return 2.4 * x

def get_fft(data, fs):
	time = np.arange(0, data.shape[0]/fs, 1/fs)
	fourier_transform = np.fft.rfft(data)
	abs_fourier_transform = np.abs(fourier_transform)
	power_spectrum = np.square(abs_fourier_transform)
	power_spectrum = 10 * np.log10(power_spectrum)
	frequency = np.linspace(0,fs/2, len(power_spectrum))
	return power_spectrum, frequency

dt = 0.005
t1 = 2330 #2330 # 2290
t2 = 2370
#time = np.arange(t1,t2,dt)
time = np.arange(0,Das_nc.shape[0]*dt,dt)

start = 0
channel = 463 - start
das = np.flip(Das_nc[int(t1/dt):int(t2/dt),start:start+channel].T,0) # the purpose of this weird line is to simplify the imshow so that the starting point of offset and time is in the bottom left (trust me)
starttime = Das_nc.time.data[0].astype('datetime64[ms]')

# frequency
index = 200
fs = 1/dt
power_spectrum, frequency = get_fft(normalize(das[index,:]), fs)

#trace_filtered = filter_bandpass(das[200,:], 1, 5, 1/dt, order=2)
fig, ax = plt.subplots(3, figsize=(7,10), sharex=False, gridspec_kw={'height_ratios': [0.3, 0.3, 2]})
fig.suptitle('Start time: ' + str(starttime + np.timedelta64(t1,'s')))

ax[0].plot(frequency,power_spectrum,color='black',linewidth=0.5)
#ax[0].set_title('Spectrum')
ax[0].set_xlabel('Frequency (Hz)')
ax[0].set_ylabel('Power (dB)')
ax[0].set_xlim(0,25)
ax[0].grid()

ax[1].plot(time[int(t1/dt):int(t2/dt)], das[index,:], c='k', linewidth=0.5) #trace_filtered)
ax[1].set_ylabel('Amplitude\n[$10^{-9}$ strain $s^{-1}$]')
ax[1].text(t1 + 1.5, 300,f'Offset = {index*2.4} m', bbox={'facecolor':'white','pad':5})
ax[1].set_xlim(t1,t2)
ax[1].grid()
ax[1].set_xlabel('Time (s)')


extent = [time[int(t1/dt)],time[int(t2/dt)],start,start+channel]
a = ax[2].imshow(das, aspect='auto', vmin=np.min(das.data)/100, vmax=np.max(das.data)/100, extent=extent, cmap='seismic')
y2 = ax[2].secondary_yaxis('right', functions=(get_offset, get_channel))
y2.set_ylabel('Offset (m)')

#ax[1].plot([time[int(t1/dt)],time[int(t2/dt)]], [200,200], 'k-', lw=1)
ax[2].axhline(y=index, c='k')
ax[2].set_ylabel('Channel')
ax[2].set_xlabel('Time (s)')
colorbar = plt.colorbar(a,orientation="horizontal", pad=0.09, location='bottom')
colorbar.set_label(label='Amplitude [$10^{-8}$ strain $s^{-1}$]', rotation=0)
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #pad=0.4,w_pad=0.5,h_pad=1.0)
#fig.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

