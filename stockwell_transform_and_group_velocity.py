#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#export HDF5_USE_FILE_LOCKING=FALSE

"""
Created on Wed Feb 30 16:08:00 2022

@author: mozef

This code does:
- Compute stockwell transform of a signal
- Auto-pick the group velocity

"""

import numpy as np
from scipy.signal import chirp
import matplotlib.pyplot as plt
from stockwell import st
from scipy import signal
from scipy.interpolate import make_interp_spline, BSpline

# ============================================================ Functions ========================================================

def S_transform(trace,fmin,fmax,df):
	#dt = Das_nc[1,200].time - Das_nc[0,200].time
	#dt = float(dt) * 1e-9
	#df = 1/dt
	fmin_samples = int(fmin/df)
	fmax_samples = int(fmax/df)
	Stransform = st.st(trace, fmin_samples, fmax_samples)
	return Stransform

def inverse_S_transform(trace,fmin,fmax,df):
	#dt = Das_nc[1,200].time - Das_nc[0,200].time
	#dt = float(dt) * 1e-9
	#df = 1/dt
	fmin_samples = int(fmin/df)
	fmax_samples = int(fmax/df)
	iStransform = st.ist(trace, fmin_samples, fmax_samples)
	return iStransform

# ========================================= Test stockwell transform using synthetic signal =====================================

# S-Transform of chirp signal
import numpy as np
from scipy.signal import chirp
import matplotlib.pyplot as plt
from stockwell import st

t = np.linspace(0, 10, 5001)
w = chirp(t, f0=12.5, f1=2.5, t1=10, method='linear')

fmin = 0  # Hz
fmax = 25  # Hz
df = 1./(t[-1]-t[0])  # sampling step in frequency domain (Hz)
fmin_samples = int(fmin/df)
fmax_samples = int(fmax/df)
stock = st.st(w, fmin_samples, fmax_samples)
extent = (t[0], t[-1], fmin, fmax)

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(t, w)
ax[0].set(ylabel='amplitude')
ax[1].imshow(np.abs(stock), origin='lower', extent=extent)
ax[1].axis('tight')
ax[1].set(xlabel='time (s)', ylabel='frequency (Hz)')
plt.show()

# ====================================================== Preparing tf-PWS results ================================================

# Time-frequency phase weighted stack
window = 400	# no need to change
dt = 0.02	# 2020 data
t = np.arange(-window,window,dt)
time_lim1 = -10   # edit here for the time boundary
time_lim2 = 10	# edit here for the time boundary
start = 100
index = 150	# pair channel index
distance = (index - start) * 2.4 # in meters

from obspy import read
file = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-10Hz_Window400_mid100_gl10/Correlogram_{start}_{index}_phase/ts_pws_tspws_pcc.sac'
#data = read(file)[0].data
data = np.flip(read(file)[0].data)
plt.figure()
plt.plot(t,data)
plt.xlim(-10,10)
plt.show()


# shifting zero lagtime to a certain mid point
zero_mid = 100 # the id of channel that will be zero point
distance = abs((index - zero_mid)) * 2.4 # in meters
file_lag = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-10Hz_Window400_mid100_gl10/Correlogram_{start}_{zero_mid}_phase/ts_pws_tspws_pcc.sac'
data_lag = read(file_lag)[0].data
lag = len(data_lag)/2 - np.where(data_lag==max(data_lag))[0][0] # basically midpoint subtracted by lcoation of the peak in sampling point domain
#lag = 150

boundary1 = int(((time_lim1 + window) / (window * 2)) * data.shape[0])
boundary2 = int(((time_lim2 + window) / (window * 2)) * data.shape[0])

# applying zero lagtime shift to the data
trace = data[int(boundary1-lag):int(boundary2-lag)]
t_cut = t[boundary1:boundary2]

time_lim1_transform = 0.25
time_lim2_transform = 8.3
boundary1 = int(((time_lim1_transform + 10) / (10 * 2)) * trace.shape[0])
boundary2 = int(((time_lim2_transform + 10) / (10 * 2)) * trace.shape[0])
trace_transform = trace[boundary1:boundary2]
t_cut_transform = t_cut[boundary1:boundary2]

v_group = distance / t_cut

plt.figure()
plt.plot(t_cut_transform,trace_transform)
#plt.xlim(-10,10)
plt.show()


v_group_plot = v_group[520:920]
v = [round(i,2) for i in v_group_plot]
t = [distance/i for i in v]

plt.figure()
plt.plot(t,v)
plt.title(f'Velocity over time for distance = {distance} m')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.show()


# ============================================= Stockwell transform using tf-PWS results ========================================

fmin = 1
fmax = 10
df = 1. / (t_cut_transform[-1] - t_cut_transform[1])

Stransform = S_transform(trace_transform,fmin,fmax,df)
iStransform = inverse_S_transform(Stransform,fmin,fmax,df)

f = np.linspace(fmin,fmax,Stransform.shape[0])
t_plot = np.linspace(t_cut_transform[-1],t_cut_transform[0],Stransform.shape[1])
v = distance / t_plot


# =================================================== Auto-picking group velocity ===============================================

cut = 100 # picking is only conducted above this value in y-axis. the number depends on the length of y-axis, i.e. velocity in s-transform matrix
Stransform_plot = np.flip(np.abs(Stransform).T,axis=0)

v_peak = [] # preparing an array to save the index at stransform matrix that has the highest value
v_low = []  # saving the index of lower boundary velocity
v_high = [] # saving the index of higher boundary velocity

for idx_f in range(Stransform_plot.shape[1]):
	a = np.where(Stransform_plot[cut:,idx_f] == np.max(Stransform_plot[cut:,idx_f]))
	bound = np.where((Stransform_plot[cut:,idx_f] > 0.95 * np.max(Stransform_plot[cut:,idx_f])) & (Stransform_plot[cut:,idx_f] < np.max(Stransform_plot[cut:,idx_f])))
	if idx_f == 30:
		print(bound)
		print(a)
	b = v[a[0][0]+cut]
	bound_high = v[bound[0][-1]+cut]
	bound_low = v[bound[0][0]+cut]
	v_peak.append(b)
	v_low.append(bound_low)
	v_high.append(bound_high)

def interpolate(f,v):
	f_int = np.linspace(f[0], f[-1], 200)
	interpolate = make_interp_spline(sorted(f), v, k=3)
	v_int = interpolate(f_int)
	return f_int, v_int

f_int, v_peak = interpolate(f,np.array(v_peak))
f_int, v_low = interpolate(f,np.array(v_low))
f_int, v_high = interpolate(f,np.array(v_high))

import pandas as pd
left = 40      # add left and right boundary for the picking. this value depends on the length of the velocity array (v_peak, v_low, v_high)
right = 170
d = {'v_peak':v_peak[left:right], 'v_low':v_low[left:right], 'v_high':v_high[left:right],'f':f_int[left:right]}
df = pd.DataFrame(data=d)
#df.to_csv(f'/home/mozef/Documents/scripts/group_vel/{zero_mid}_{index}_start{start}_{fmin}-10Hz.txt', index=None, sep=' ', mode='a')

# =============================================================== Plotting =================================================================

def normalize(trace):	# normalize trace
		return trace / max(trace)

extent = (fmin, fmax, t_cut_transform[-1], t_cut_transform[0])

from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

a = Stransform_plot.max()
Stransform_plot_norm = np.divide(Stransform_plot, a)

fig, ax = plt.subplots(2, 1, sharex=False, gridspec_kw={'height_ratios': [1.5, 2]},figsize=(10,10))
plt.suptitle(f'Stockwell Transform of TF-PWS Correlogram\n Channel {zero_mid} and {index}\nFilter 1-10 Hz $-$ Distance = {round(distance,2)} m')
ax[0].plot(t_cut, normalize(trace), color='black')
ax[0].set(ylabel='normalized amplitude')
ax[0].set(xlabel='time (s)')
ax[0].set_xlim(0,time_lim2)
a = ax[1].pcolor(f,v,Stransform_plot_norm, norm=colors.LogNorm(vmin=Stransform_plot_norm.min()*10000, vmax=Stransform_plot_norm.max()))
#a = ax[1].pcolor(f,v,Stransform_plot_norm)
ax[1].axis('tight')
ax[1].set(xlabel='frequency (Hz)', ylabel='Group Velocity (m/s)')
colorbar = plt.colorbar(a,orientation="horizontal", pad=0.2, location='bottom', aspect=100)
colorbar.set_label('normalized amplitude', rotation=0)
ax[1].plot(f_int[left:right],v_peak[left:right], color='red')
ax[1].plot(f_int[left:right],v_low[left:right],'--', color='orange')

ax[1].plot(f_int[left:right],v_high[left:right],'--', color='orange')

a = time_lim1_transform
b = time_lim2_transform
ax[0].axvline(a, color='red')
ax[0].axvline(b, color='red')
#ax[0].text(a, 1, a, size=8)
#ax[0].text(b, 1, b, size=8)

# pick value (in case we want to pick it manually instead) 
#print("Pick the strong amplitude:\nLeft button to pick\nRight button to cancel last pick\nMiddle button to stop (max=100 picks)")
#x = plt.ginput(100)
#print(x)
#v_pick = [i[0] for i in x]
#f_pick = [i[1] for i in x]
#plt.plot(v_pick,f_pick,'--r')
plt.show()



