#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#export HDF5_USE_FILE_LOCKING=FALSE

"""
Created on Wed Mar 2 16:53:00 2022

@author: mozef

This code does:
- Cross-correlate DAS strain data using PCC or aPCC
- Stack cross-correlogram using tf-PWS

"""

#export HDF5_USE_FILE_LOCKING=FALSE
from xdas.io.febus import read as read_das
from icoords import InterpolatedDataArray
from xdas.io.febus import correct_gps_time
import matplotlib.pylab as plt
import numpy as np
from scipy import signal
from tools import *
from stack_functions import *
#from stockwell import st
import glob



# =================================================== Cross-correlation ========================================================
# define variable for correlation
window = 300 			# duration of time window in second
lag = window/2 			# shift of starting time on the next iteration
time = 0 			# time initialization

# filter boundary (Hz)
lowcut = 1
highcut = 10
fs = 50 # Hz

start = 100		# starting/reference point. currently used: 100, 200, 320, 400
channel = 150		# number of channel used			 
start_position = 'mid_all' 	# 'mid', 'mid all', 'early', or 'end'. If mid, channel is defining number of channel after the start channel and 					   before the start channel, thus the number of channel will be 2*n+1. If mid_all, all channel is 					   cross-correlated to the channel indicated by start_position (the number of channel is not used because it's 					   already considering all channel)

files = glob.glob("/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/*.nc")
folder_destination = '/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/Auto-correlation/test/' # please make the folder before running

change_gauge_length = False
gauge_length = 50

# open one file to get sampling time information
iDas_nc = InterpolatedDataArray.from_netcdf("/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/SR_Stromboli_2020-09-21_12-02-26_UTC_decimated50Hz.nc") 
Das_nc = iDas_nc.load_icoords()
dt = Das_nc[1,200].time - Das_nc[0,200].time
dt = float(dt) * 1e-9

# make output array
if start_position == 'mid':
	cc_normal = np.zeros((2*int(window/dt),2*channel+1))
	cc_onebit = np.zeros((2*int(window/dt),2*channel+1))
	cc_phase = np.zeros((2*int(window/dt),2*channel+1))
elif start_position == 'mid_all':
	cc_normal = np.zeros((2*int(window/dt),430+1))
	cc_onebit = np.zeros((2*int(window/dt),430+1))
	cc_phase = np.zeros((2*int(window/dt),430+1)) # change to 100-57+1 if only want from 57 to 100
else:
	cc_normal = np.zeros((2*int(window/dt),channel+1))
	cc_onebit = np.zeros((2*int(window/dt),channel+1))
	cc_phase = np.zeros((2*int(window/dt),channel+1))


f = open('error.txt','w') # make file containing files with few channels

#for i in range(len(files)):
for i in range(20):
	print(f'\nFile {i+1} / {len(files)}') 
	print(f'\nOpening the file {files[i]} ...')
	iDas_nc = InterpolatedDataArray.from_netcdf(files[i]) 
	Das_nc = iDas_nc.load_icoords()
	if Das_nc.shape[1] < 463: 
		print(f'This file only has {Das_nc.shape[1]} channels. Skipping this file')
		f.write(f'{files[i]}	{Das_nc.shape[1]} channels\n')
		continue	# skipping file if the channel is not complete

	
	
	print(f'Filtering using bandpass filter between {lowcut} and {highcut} Hz ...')

	if change_gauge_length == False:
		traces_filt = Das_nc.copy()
	else:
		
		# to try GL = 10
		traces_filt = Das_nc.copy()
		traces_filt[:,0] = Das_nc[:,0] + Das_nc[:,1] # the first channel with gl = 10 is defined as channel 0 + channel 1
		traces_filt[:,Das_nc.shape[1]-1] = Das_nc[:,Das_nc.shape[1]-2] + Das_nc[:,Das_nc.shape[1]-1] # the last channel with gl = 10 is defined as channel N-1 + channel N
		print(f'Converting data to using gauge length = 10')
		#for i in range(1,traces_filt.shape[1]-2): # each end is skipped because it has been calculated in the previous line
		for i in range(57,60):
			traces_filt[:,i] = Das_nc[:,i-1] + Das_nc[:,i] + Das_nc[:,i+1]
		#traces_filt_gl10 = traces_filt[:,::2]
		

		"""
		# to try GL = 20
		traces_filt = Das_nc.copy()
		print(f'Converting data to using gauge length = 20')
		for i in range(0,traces_filt.shape[1]-1): # each end is skipped because it has been calculated in the previous line
			try:
				traces_filt[:,i] = Das_nc[:,i-3] + Das_nc[:,i-2] + Das_nc[:,i-1] + Das_nc[:,i] + Das_nc[:,i+1] + Das_nc[:,i+2] + Das_nc[:,i+3]
			except:
				None
		#traces_filt_gl20 = traces_filt[:,::3]
		"""

	for j in [57,430]:
		try:
			#traces_filt[:,j] = filter_bandpass(Das_nc[:,j], lowcut, highcut, fs,order=4) # scipy.signal filter
			traces_filt[:,j] = bandpass(Das_nc[:,j],lowcut,highcut,fs, zerophase=True) # obspy filter (faster)
		except:
			None
		
	length = traces_filt.shape[0]	# total length of data

	# cross correlate

	print('\nStarting the phase auto-correlation  ...')
	cc_phase_auto, t = auto_correlate(cc_phase,traces_filt,folder_destination,window,lag,length,time,channel,
			start,start_position,onebit=False,phase_cc=True,decimate=1,resample=1)

	#print('\nStarting the phase cross-correlation  ...')
	#cc_phase, t = cross_correlate(cc_phase,traces_filt,folder_destination,window,lag,length,time,channel,
	#		start,start_position,onebit=False,phase_cc=True,decimate=1,resample=1)

	#print('\nStarting the normal cross-correlation  ...')
	#cc_normal, t = cross_correlate(cc_normal,traces_filt,window,lag,length,time,channel,
	#		start,start_position,onebit=False,phase_cc=False,decimate=decimate,resample=resample)

	#cc_onebit, t = cross_correlate(cc_onebit,traces_filt,window,lag,length,time,channel,
	#		start,start_position,onebit=True,phase_cc=False,decimate=decimate,resample=resample)



f.close()


# =================================================== tf-PWS ========================================================

# tf-pws using eleonore's code
from obspy import read
from obspy import UTCDateTime

start = 100
channel = 300
start_position = 'mid_all'

folder = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/Auto-correlation/test/'
#folder = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/01-10Hz_Window300_mid100_gl10'

for i in range(57,430):
#for i in range(100,103):
	os.system(f'cp /home/mozef/Documents/scripts/ts_pws {folder}/Correlogram_{start}_{i}_phase/')
	os.chdir(f'{folder}/Correlogram_{start}_{i}_phase/')
	os.system('find . -type f -empty -print -delete') # delete empty file
	print(f'tf-pws for correlogram {start} and {i}')
	os.system('ls *.sac > tf-pws.in')
	os.system('./ts_pws tf-pws.in osac="tspws_pcc"')
	# there will be a warning that says the header is corrupted, just ignore it since we just need the cross-correlogram time array, the header information is not important



