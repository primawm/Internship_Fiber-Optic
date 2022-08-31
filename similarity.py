#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 10:09:00 2022

@author: mozef
"""

#export HDF5_USE_FILE_LOCKING=FALSE
import matplotlib.pylab as plt
import numpy as np
from scipy import signal
import glob
from obspy import read
from scipy.stats import pearsonr
import random

# reading 1 file
fname1 = '/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-5Hz_Window400_start100/Correlogram_100_100_normal/correlogram_100_100_0_400_20200919_1657.sac'
fname2 = '/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-5Hz_Window400_start100/Correlogram_100_100_phase/correlogram_100_100_0_400_20200916_1256.sac'
file_normalcc = read(fname1)
file_phasecc = read(fname2)

file_normalcc.plot()
file_phasecc.plot()

# stack the content of folder
def stack_write(folder,x1,x2,mode='normal'):
	# to account for the station pair input
	files = glob.glob(f'{folder}/Correlogram_{str(x1)}_{str(x2)}_{mode}/*.sac')
	print(files)
	stacked = np.zeros(len(read(files[0])[0].data))
	for file in files:
		a = read(file)
		stacked = stacked + a[0].data
	a = read(files[0])
	a[0].data = stacked
	#a.write(f'{folder}/Correlogram_{str(x1)}_{str(x2)}_{mode}/stack_{mode}.sac', format='SAC')
	return stacked

def stack(folders):
	stacked = np.zeros(len(read(folders[0])[0]))
	for file in folders:
		try:
			a = read(file)
			stacked = stacked + a[0].data
		except:
			None
	return stacked

def normalize(trace):	# normalize trace
		return trace / max(trace)

x1 = 100
x2 = 100
folder = f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-5Hz_Window400_start100/'
stack_all_normal = stack_write(folder,x1,100,mode='normal')
stack_all_phase = stack_write(folder,100,100,mode='phase')
stack_tfpws_all = read(f'/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-5Hz_Window400_start100/Correlogram_100_100_phase/ts_pws_tspws_pcc.sac')

window = 400
dt = 0.02
t = np.arange(-window,window,dt)

plt.figure()
plt.plot(t, normalize(stack_all_normal), label='normal')
plt.plot(t, normalize(stack_all_phase), label='phase')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Normalized Amplitude')
plt.xlabel('Time lag (s)')
plt.ylabel('Normalized Amplitude')
plt.show()

# cut trace
time_lim1 = -20
time_lim2 =20
boundary1 = int(((time_lim1 + window) / (window * 2)) * len(stack_all_normal))
boundary2 = int(((time_lim2 + window) / (window * 2)) * len(stack_all_normal))

stack_cut_normal_all = stack_all_normal[boundary1:boundary2]
stack_cut_phase_all = stack_all_phase[boundary1:boundary2]
stack_cut_tfpws_all = stack_tfpws_all[0][boundary1:boundary2]

t_cut = t[boundary1:boundary2]

files_phase = glob.glob('/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-5Hz_Window400_start100/Correlogram_100_100_phase/*')
a = read(files_phase[0]) # check one trace
b = a[0].data[boundary1:boundary2]

plt.figure(figsize=(10,8))
plt.plot(t_cut, normalize(stack_cut_normal_all), label='Normal CC - Linear stack')
plt.plot(t_cut, normalize(stack_cut_phase_all), label='Phase CC - Linear stack')
plt.plot(t_cut, normalize(stack_cut_tfpws_all), label='Phase CC - TF-PWS')
plt.plot(t_cut, normalize(b), label='Phase CC - 1 Random trace')
plt.legend()
plt.title(f'Stack of 2020 DAS Stromboli data\nFilter 1-5 Hz\nChannel {x1} and {x2}')
plt.xlabel('Lag time (s)')
plt.ylabel('Normalized Amplitude')
plt.show()

# asses SNR
signal = stack_all_phase
def signalnoiseratio(signal):
	window = 400
	# signal boundary
	time_lim1 = -5
	time_lim2 =5
	boundary1 = int(((time_lim1 + window) / (window * 2)) * len(stack_all_normal))
	boundary2 = int(((time_lim2 + window) / (window * 2)) * len(stack_all_normal))

	# noise boundary
	time_noise1 = -300
	time_noise2 = -200 
	boundary_noise1 = int(((time_noise1 + window) / (window * 2)) * len(stack_all_phase))
	boundary_noise2 = int(((time_noise2 + window) / (window * 2)) * len(stack_all_phase))

	# from sabra2005
	#sd = np.std(np.gradient(stack_all_phase[boundary_noise1:boundary_noise2]))
	#max_amp = np.max(stack_all_phase)

	rms_signal = np.sqrt(np.mean(signal[boundary1:boundary2]**2))
	rms_noise = np.sqrt(np.mean(signal[boundary_noise1:boundary_noise2]**2))
	snr = 20 * np.log(rms_signal/rms_noise)
	return snr

snr = signalnoiseratio(signal)


# assess similarity of stack given several number of signal
import os
import shutil

# signal boundary
window = 400

time_lim1 = -5 # in second. this is the time boundary to cut the correlogram
time_lim2 = 5
boundary1 = int(((time_lim1 + window) / (window * 2)) * len(stack_all_normal)) # get the position of time boundary in sampling domain
boundary2 = int(((time_lim2 + window) / (window * 2)) * len(stack_all_normal))

# correlogram folder
folder_all = '/net/runic/moby/data/projets/monidas/biagioli/Prima/Correlogram/1-10Hz_Window400_mid100'

num_signal1 = 1
num_signal2 = 901
step = 100
x = np.arange(num_signal1,num_signal2+step,step)
x = [1, 21, 41, 61, 81, 101, 201, 301, 401, 501, 601, 701, 801, 901]
iteration = 20

# initialization of plotting variable
coeffs_normal = np.zeros((iteration,len(x)))
coeffs_phase = np.zeros((iteration,len(x)))
coeffs_tfpws = np.zeros((iteration,len(x)))

snr_normal = np.zeros((iteration,len(x)))
snr_phase = np.zeros((iteration,len(x)))
snr_tfpws = np.zeros((iteration,len(x)))

start = 100
pair = 200

for i in range(iteration):
	for k in range(len(x)):
		folder_phase = f'{folder_all}/Correlogram_{start}_{pair}_phase'
		files_phase = glob.glob(f'{folder_phase}/*.sac')

		#folder_normal = f'{folder_all}/Correlogram_{start}_{pair}_normal'
		#files_normal = glob.glob(f'{folder_phase}/*.sac')
		#stack_all_normal = stack(files_normal)
		#stack_cut_normal_all = stack_all_normal[boundary1:boundary2]
		#files_normal_sample = random.sample(files_normal, k=x[k])
		#stack_normal = stack(files_normal_sample)
		#stack_cut_normal = stack_normal[boundary1:boundary2]
		#snr_phase[i,k] = signalnoiseratio(stack_phase)
		#snr_normal[i,k] = signalnoiseratio(stack_normal)
		#coeff_normal = np.correlate(normalize(stack_cut_normal), normalize(stack_normal))
		#coeff_normal = pearsonr(stack_cut_normal_all, stack_cut_normal)
		#coeffs_normal[i,k] = coeff_normal[0]
		
		print(f'Phase linear stack iteration {i}\n{folder_all}')
		stack_all_phase = stack(files_phase)
		stack_tfpws_all = read(f'{folder_phase}/ts_pws_tspws_pcc.sac')[0] # this 0 index is for accessing the trace from stream
		stack_cut_phase_all = stack_all_phase[boundary1:boundary2]
		stack_cut_tfpws_all = stack_tfpws_all[boundary1:boundary2]

		files_phase_sample = random.sample(files_phase, k=x[k])
		
		stack_phase = stack(files_phase_sample)
		stack_cut_phase = stack_phase[boundary1:boundary2]

		#coeff_phase = np.correlate(normalize(stack_cut_phase), normalize(stack_phase))
		coeff_phase = pearsonr(stack_cut_phase_all, stack_cut_phase)

		#print(coeff_phase[0], snr_phase)

		coeffs_phase[i,k] = coeff_phase[0]
		print(f'Phase TF-PWS stack iteration {i}')

		# tf-pws
		try:
			os.mkdir(f'{folder_phase}/tmp')
		except:
			None

		os.system(f'cp /home/mozef/Documents/scripts/ts_pws {folder_phase}/tmp')
		
		for file in files_phase_sample:
				shutil.copy(file, f'{folder_phase}/tmp')
			
	
		os.chdir(f'{folder_phase}/tmp')

		print(f'tf-pws of {x[k]} correlogram files')
		os.system('ls *.sac > tf-pws.in')
		os.system('./ts_pws tf-pws.in osac="tspws_pcc"')

		stack_tfpws = read(f'{folder_phase}/tmp/ts_pws_tspws_pcc.sac')[0] # this 0 index is for accessing the trace from stream
		stack_cut_tfpws = stack_tfpws[boundary1:boundary2]
		coeff_tfpws = pearsonr(stack_cut_tfpws_all, stack_cut_tfpws)
		coeffs_tfpws[i,k] = coeff_tfpws[0]

		#snr_tfpws[i,k] = signalnoiseratio(stack_tfpws)

		os.system(f'rm {folder_phase}/tmp/*')


#coeffs_normal_plot = abs(coeffs_normal)
coeffs_phase_plot = abs(coeffs_phase)
coeffs_tfpws_plot = abs(coeffs_tfpws)

#avr_normal = [np.average(coeffs_normal_plot[:,i]) for i in range(coeffs_normal_plot.shape[1])] 
#avr_phase = [np.average(coeffs_phase_plot[:,i]) for i in range(coeffs_phase_plot.shape[1])] 
#avr_tfpws = [np.average(coeffs_tfpws_plot[:,i]) for i in range(coeffs_tfpws_plot.shape[1])] 
avr_phase = [np.median(coeffs_phase_plot[:,i]) for i in range(coeffs_phase_plot.shape[1])] 
avr_tfpws = [np.median(coeffs_tfpws_plot[:,i]) for i in range(coeffs_tfpws_plot.shape[1])] 

plt.figure(figsize=(8,6))
for i in range(iteration):
	#plt.scatter(x,coeffs_normal_plot[i,:], s=1, color='cyan')
	plt.scatter(x,coeffs_phase_plot[i,:], s=6, color='red')	
	plt.scatter(x,coeffs_tfpws_plot[i,:],s=6,color='green')		

#plt.plot(x,np.round(avr_normal,2), color='cyan', label ='Normal CC')
plt.plot(x,avr_phase, color='red', label ='Phase CC_Linear stack')
plt.plot(x,avr_tfpws, color='green', label ='Phase CC_TF-PWS')
plt.xlabel('Num. of signals')
plt.ylabel('Similarity')
plt.legend(loc='lower right')
plt.ylim(-0.05,1.05)
plt.xlim(-30,950)
plt.xticks(np.arange(0,950,100))
plt.title(f'DAS Stromboli 2020\nChannel {start} and {pair} $-$ Filter 1-10 Hz')
plt.show()

'''
# assess the snr over number of data
plt.figure(figsize=(8,6))
for i in range(iteration):
	#plt.scatter(x,snr_normal[i,:], s=1, color='cyan')
	plt.scatter(x,snr_phase[i,:], s=1, color='red')	
	plt.scatter(x,snr_tfpws[i,:], s=1, color='green')
	
avr_normal = [np.average(snr_normal[:,i]) for i in range(snr_normal.shape[1])] 
avr_phase = [np.average(snr_phase[:,i]) for i in range(snr_phase.shape[1])] 
avr_tfpws = [np.average(snr_tfpws[:,i]) for i in range(snr_tfpws.shape[1])] 

#plt.plot(x,np.round(avr_normal,2), color='cyan', label ='Normal CC')
plt.plot(x,np.round(avr_phase,2), color='red', label ='Phase CC_Linear stack')
plt.plot(x,np.round(avr_tfpws,2), color='green', label ='Phase CC_TF-PWS')
plt.xlabel('Num. of signals')
plt.ylabel('SNR (dB)')
plt.legend(loc='lower right')
#plt.ylim(0.6,1.05)
plt.title(f'DAS Stromboli 2020 - SNR over number of channel\n3-8 Hz\nStart channel = {start}, End channel = {pair}')
plt.show()
'''

