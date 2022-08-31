#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 15:09:00 2022

@author: mozef
"""
#export HDF5_USE_FILE_LOCKING=FALSE
from xdas.io.febus import read as read_das
from icoords import InterpolatedDataArray
import matplotlib.pylab as plt
import numpy as np
from scipy import signal
import glob
from obspy import read

files = glob.glob('/net/runic/moby/data/projets/monidas/biagioli/Prima/DAS_2020_50Hz_nc/*2020*')
 

# import data
iDas_nc = InterpolatedDataArray.from_netcdf(files[0]) 
Das_nc = iDas_nc.load_icoords()

time0 = Das_nc.time.data[0].astype('datetime64[ms]')
time1 = Das_nc.time.data[-1].astype('datetime64[ms]')

duration = time1 - time0
duration = duration.astype('timedelta64[s]')

# gather all starttime information
duration = []
starttime = []
for file in files:
	iDas_nc = InterpolatedDataArray.from_netcdf(file) 
	Das_nc = iDas_nc.load_icoords()

	time0 = Das_nc.time.data[0].astype('datetime64[ms]')
	time1 = Das_nc.time.data[-1].astype('datetime64[ms]')

	d = time1 - time0
	d = d.astype('timedelta64[s]')

	duration.append(d)
	starttime.append(str(time0))

starttime, duration = zip(*sorted(zip(starttime, duration)))

y = np.array(duration)
x = np.arange(0,len(starttime))
x_ticks_labels = starttime

fig, ax = plt.subplots(1,1, figsize=(8,6)) 
ax.bar(x,y,color='black')
ax.set_xticks(x)
ax.set_xticklabels(x_ticks_labels, rotation='vertical', fontsize=10)
#plt.xlabel('Start time of File')
plt.ylabel('Duration (s)')
plt.title('DAS Stromboli 2020\nData Availability',size=12)
plt.tight_layout()
plt.locator_params(axis='x', nbins=20)
plt.show()



# brick view
import datetime

# to check if the time is between two other
def time_in_range(start, end, current):
    """Returns whether current is in the range [start, end]"""
    return start <= current <= end

dt = np.timedelta64(60,'s')
current = np.datetime64('2020-09-16T00:00:00.00')

file_idx = 0

day_num = 24-16
num = int(np.timedelta64(24,'h')/dt) # 24 hour divided by 600 s
matrix = np.zeros((num+1,day_num))


for j in range(matrix.shape[1]):
	start = current
	minutes = 0

	while current - start < np.timedelta64(24,'h'):
	#while current <= np.datetime64('2020-09-24T24:00:00.00'):
		print(minutes)
		try:
			file_start = np.datetime64(starttime[file_idx])
			file_end = file_start + duration[file_idx]
			file_start_next = np.datetime64(starttime[file_idx+1])
		except:	
			matrix[minutes,j] = False
			

		if current < np.datetime64(starttime[0]):
			print('not yet')
			matrix[minutes,j] = False
		else:
			if time_in_range(file_start,file_end,current) == True:
				print(file_end,file_start_next,current)
				print('yes')
				matrix[minutes,j] = True
			else:
				if time_in_range(file_end,file_start_next,current) == True:
					#print(file_end,file_start_next,current)
					print('no data')
					matrix[minutes,j] = False
				else:
					matrix[minutes,j] = False
					print('changing file')
					file_idx = file_idx + 1
	
		
		minutes = minutes+1
		current = current + dt


plt.figure(figsize=(8,4))
plt.pcolor(matrix[:,2:].T, cmap='Greys')
#plt.colorbar()
day_axis = np.arange(0,7,1)
plt.gca().set_yticks(day_axis[:-1]+0.5)
plt.gca().set_yticklabels(['18 Sep','19 Sep','20 Sep','21 Sep','22 Sep','23 Sep'])
time_axis = np.arange(0,24,1)
plt.gca().set_xticks(np.linspace(0,60*24,24))
plt.gca().set_xticklabels(time_axis)
plt.xlabel('Hour')
plt.title('DAS Stromboli 2020\nData Availability')
#plt.gca().set_xlim(0,24)
plt.grid(color='black')
plt.show()


