#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 11:16:00 2022

@author: mozef
"""
#export HDF5_USE_FILE_LOCKING=FALSE

import matplotlib.pylab as plt
import numpy as np
import pandas as pd


a = pd.read_csv('/home/mozef/Documents/station/GPS_UTM.csv')  # offset information
offset = a['channel']
offset_new = ['d' + str(int(offset[i][1:4]) + 10) for i in range(len(offset))]
#offset_label = offset_new[start+57:start+channel+57+1]

lon = a['Lon']
lat = a['Lat']

plt.figure()
plt.plot(lon,lat)
plt.xlabel('UTM X (m)')
plt.ylabel('UTM Y (m)')
plt.show()

plt.figure()
plt.scatter(lon,lat, s=0.5)
plt.xlabel('UTM X (m)')
plt.ylabel('UTM Y (m)')
plt.show()

n = offset_new
fig, ax = plt.subplots(figsize=(15,15))
ax.scatter(lon, lat)
for i, txt in enumerate(n):
    ax.annotate(txt, (lon[i], lat[i]), size=8)

#ax.set_aspect('equal')
#plt.xlim(519450,519650)
#plt.ylim(200+4.294e6,400+4.294e6)
plt.xlabel('UTM X (m)')
plt.ylabel('UTM Y (m)')
plt.show()
