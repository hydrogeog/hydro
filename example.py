#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd, numpy as np, pyproj, matplotlib.pyplot as plt
import hydro 

# calculate rating curve for give flows
stage = np.array([0.96, 1.02, 1.218, 1.313, 1.186, 1.215,
                  2.121, 1.638, 2.163, 1.902, 2.008])
discharge = np.array([0.43325, 0.5345, 0.4995, 1.26475, 0.5385, 0.9645,
                  25.391, 8.6795, 30.348, 17.604, 22.2395])
r_curve = hydro.RC(stage, discharge)
# rating curve plot and coefficients
r_curve.plot(log=False)
coef = r_curve.popt

# read flow data csv
flow = pd.read_csv("flow.csv")
flow.TimeStamp = pd.to_datetime(flow.TimeStamp)

# convert stage to discharge using rating curve equation from above
flow['discharge_cfs'] = r_curve.Q(flow.Level_ft)
# create Discharge object
flowdata = hydro.Discharge(flow.TimeStamp, flow.discharge_cfs, flow.Rainfall_in)
# find daily mean discharge
dailyMean = flowdata.dailyQ()

# plot daily mean data with 15-minute data
# this is only intended as a quick plot, not save-quality
flowdata.plot([dailyMean.day, dailyMean.meanQ])

# Richards-Baker Flashiness Index
RB_index = flowdata.RB_Flashiness()
print(f'Richards-Baker Flashiness: {RB_index:.4f}')

# Flow duration curve
# if plot=True, you need to pip install probscale
flowdur = flowdata.flow_duration(plot=True)

# Baseflow separation using Lyne_Hollick method
# 3 forward passes with the digital filter
flow['LH_bflow'] = flowdata.Lyne_Hollick(direction='fff')
flowdata.plot([flow.TimeStamp, flow.LH_bflow], log=False)

# Baseflow separation using Eckhardt method
# re-initialize discharge object to erase baseflow from above
flowdata = hydro.Discharge(flow.TimeStamp, flow.discharge_cfs, flow.Rainfall_in)
flow['E_bflow'] = flowdata.Eckhardt(re=3)
flowdata.plot([flow.TimeStamp, flow.E_bflow], log=False)

# flow DataFrame now has discharge and base flow calculated with two different 
# methods. Plot that DataFrame if you want to compare the methods.

################################################################################

## Geographic functions
data = pd.read_csv("stream.csv")

# Longitudintal profile creation and plot
adjusted_elevation = hydro.Profile_smoothing(data.ELEVATION, data.DISTANCE_FROM_MOUTH,
                                             plot=True)

# Reproject coordinates
NAD83_GAfeet = pyproj.Proj("+init=EPSG:2239")
def trans(lon, lat):
    lon2, lat2 = NAD83_GAfeet(lon, lat, preserve_units=True)
    return (lon2/0.3048006096012192, lat2/0.3048006096012192)
data['Easting'] = np.zeros(len(data))
data['Northing'] = np.zeros(len(data))
for i in range(len(data)):
    data.Easting.values[i], data.Northing.values[i] = trans(data.LONGITUDE.iloc[i], data.LATITUDE.iloc[i])

# Calculate sinuosity 1000 ft on either side of points that are 4 ft apart
sin = hydro.sinuosity(data.Easting, data.Northing, 500, 4)

# plot the figure
from matplotlib.ticker import NullFormatter
fig = plt.figure(figsize=(22,5))
ax1 = plt.subplot2grid((1,6), (0,0), colspan=5) # gridspec to create frequency histogram
ax1.plot(data['cng_(Feet)'], sin, lw=3, alpha=.6, label='Sinuosity')

ax1.set_ylim(1, 1.6)
plt.grid(which='both')
plt.xlabel('Distance from Downstream')
plt.ylabel('Sinuosity')
plt.title('Creek Sinuosity')
plt.legend()

ax2 = plt.subplot2grid((1,6), (0,5))
ax2.hist(sin, bins=60, orientation='horizontal', alpha=.6)
ax2.xaxis.set_major_formatter(NullFormatter())
ax2.yaxis.set_major_formatter(NullFormatter())
ax2.set_ylim(1,1.6)

plt.tight_layout()
plt.show()

# Sinuosity of the entire stream
sin_total = hydro.sinuosity([data.Easting.iloc[0], data.Easting.iloc[-1]],
                            [data.Northing.iloc[0], data.Northing.iloc[-1]],
                            data['cng_(Feet)'].max()*2, data['cng_(Feet)'].max())
print(f'Sinuosity of entire stream: {sin_total:.4f}')
