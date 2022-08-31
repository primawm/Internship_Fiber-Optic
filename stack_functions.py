#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#export HDF5_USE_FILE_LOCKING=FALSE

"""
Created on Tue Mar 1 11:21:00 2022

@author: mozef

"""

from xdas.io.febus import read as read_das
from icoords import InterpolatedDataArray
import matplotlib.pylab as plt
import numpy as np
from xdas.io.febus import correct_gps_time
from scipy import signal
from tools import *
import os
from obspy import Stream,Trace,UTCDateTime
from obspy.core import Stats
from obspy import read


# import data
#iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/DAS_2020_200Hz_nc/SR_Stromboli_2020-09-20_13-41-24_UTC_decimated200Hz.nc") # many small events and one big event
#iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/DAS_2020_200Hz_nc/SR_Stromboli_2020-09-23_10-35-49_UTC_decimated200Hz.nc") # one event?

# decimated 50Hz
iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/SR_2020-09-17_13-47-47_UTC_decimated50Hz.nc")

Das_nc = iDas_nc.load_icoords()
dt = Das_nc[1,200].time - Das_nc[0,200].time
dt = float(dt) * 1e-9

def normalize(trace):	# normalize trace
		return trace / max(trace)

def cut_traces(trace, start, duration): # netcdf trace afrom tools import *s input
	dt = trace[1,0].time - trace[0,0].time
	dt = float(dt) * 1e-9
	traces_cut = trace[int(start/dt):int(start/dt)+int(duration/dt),:]
	return traces_cut

def filter_bandpass(trace, lowcut, highcut, fs, order=2):
	nyq = 0.5 * fs
	low = lowcut/nyq
	high = highcut/nyq
	#order = 2
	b,a = signal.butter(order, [low, high], 'bandpass', analog=False)
	print(b,a)
	trace_filtered = signal.filtfilt(b, a, trace, axis=0)
	return trace_filtered

def one_bit(trace):
	trace_onebit = np.where(trace < 0, -1, trace)
	trace_onebit = np.where(trace_onebit > 0, 1, trace_onebit)
	return trace_onebit

def make_date(traces):
	date = UTCDateTime(str(traces[0].time.data))
	year = str(date.year)
	month = str(date.month).zfill(2)
	day = str(date.day).zfill(2)
	hour = str(date.hour).zfill(2)
	minute = str(date.minute).zfill(2)
	fdate = f'{year}{month}{day}_{hour}{minute}'
	return fdate

def make_dir(traces, pcc, parent_dir, folder, stats, start,i, time, window):
	try:
		os.mkdir(parent_dir + folder)
	except:
		None
		# print(f'Folder {parent_dir + folder} already exists. Continue writing ...')

	stats.npts = len(pcc)
	correlogram = Trace(data=pcc, header=stats)
	fdate = make_date(traces)
	fname = f'correlogram_{start}_{i}_{int(time*dt)}_{window}_{fdate}'
	correlogram.write(parent_dir + folder + fname + '.sac', format='sac')		

# define function for cross correlation: one-bit cc, normal cc, and phase cc
def cross_correlate(cc,traces,parent_dir,window,lag,length,time,channel,start,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	dt = Das_nc[1,200].time - Das_nc[0,200].time
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# looping cross correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([channel,cc.shape[0],cc.shape[1]])
	stats = Stats()
	stats.delta = dt
	stats.sampling_rate = 1/stats.delta
	stats.network = 'Correlogram'
	stats.station = 'STRMBL'
	stats.location = 'Italy'

	while time < length:
		if iteration%decimate == 0:
			j = 0
			print(f'Start of the {iteration}th time window. Starttime = {time*dt} s') 

			# cross correlate all selected channel in time window
			if start_position == 'mid':
				ranges = np.arange(start-channel,start+channel+1,1)
			elif start_position == 'early':
				ranges = np.arange(start, start+channel+1, 1)
			elif start_position == 'end':
				ranges = np.arange(start, start-channel-1, -1)
			elif start_position == 'mid_all':
				ranges = np.arange(57, 430+1, 1)

			#for i in ranges:
			for i in range(57,60): # use this when only interested in some channels
				if onebit == True:
					x1 = one_bit(traces[time:time+int(window/dt),start])
					x2 = one_bit(traces[time:time+int(window/dt),i])
				else:
					x1 = traces[time:time+int(window/dt),start]
					x2 = traces[time:time+int(window/dt),i]
				
				if phase_cc == True:
					try:
						_t , pcc = pcc2(x1, x2, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc
						
						folder = f'Correlogram_{start}_{i}_phase/'
						make_dir(traces, pcc, parent_dir, folder, stats, start,i, time, window)
							
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break
				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x2)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start}_{i}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start,i, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt)
		iteration = iteration + 1
	return cc, t

# define function for auto correlation: one-bit cc, normal cc, and phase cc
def auto_correlate(cc,traces,parent_dir,window,lag,length,time,channel,start,start_position,
			onebit,phase_cc,decimate=1,resample=1):
	dt = Das_nc[1,200].time - Das_nc[0,200].time
	dt = float(dt) * 1e-9
	#dt = dt*resample

	t = np.arange(-window, window, dt)

	# looping auto correlation over varying time window
	iteration = 1
	#pcc_all = np.zeros([channel,cc.shape[0],cc.shape[1]])
	stats = Stats()
	stats.delta = dt
	stats.sampling_rate = 1/stats.delta
	stats.network = 'Correlogram'
	stats.station = 'STRMBL'
	stats.location = 'Italy'

	while time < length:
		if iteration%decimate == 0:
			j = 0
			print(f'Start of the {iteration}th time window. Starttime = {time*dt} s') 

			# cross correlate all selected channel in time window
			if start_position == 'mid':
				ranges = np.arange(start-channel,start+channel+1,1)
			elif start_position == 'early':
				ranges = np.arange(start, start+channel+1, 1)
			elif start_position == 'end':
				ranges = np.arange(start, start-channel-1, -1)
			elif start_position == 'mid_all':
				ranges = np.arange(57, 430+1, 1)

			#for i in ranges:
			for i in range(57,60):
				x1 = traces[time:time+int(window/dt),i]
				
				if phase_cc == True:
					try:
						_t , pcc = apcc2(x1, dt, -window, window)
						cc[:,j] = cc[:,j] + pcc # linear stack
						
						folder = f'Correlogram_{start}_{i}_phase/'
						make_dir(traces, pcc, parent_dir, folder, stats, start,i, time, window)
							
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break
				else:
					try: # to account for cutted trace in last time window
						corr = signal.correlate(x1,x2)
						corr = np.append(corr,0) # pad with zero to similarize the array length
						cc[:,j] = cc[:,j] + corr  
						folder = f'Correlogram_{start}_{i}_normal/'
						make_dir(traces, corr, parent_dir, folder, stats, start,i, time, window)
					except:
						# print(len(corr),len(cc[:,j])) # show length
						print(f'Time window starting at {time*dt} s is shorter. Finishing the program ...')
						break

				j = j + 1
		time = int(time + lag/dt)
		iteration = iteration + 1
	return cc, t
