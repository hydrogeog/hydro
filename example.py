#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 14:13:47 2017

@author: carson
"""
import pandas as pd, numpy as np
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
flow['LH_bflow'] = flowdata.Lyne_Hollick()
flow['LH_bflow'] = flowdata.Lyne_Hollick()
flow['LH_bflow'] = flowdata.Lyne_Hollick()
flowdata.plot([flow.TimeStamp, flow.LH_bflow], log=False)

# Baseflow separation using Eckhardt method
# re-initialize discharge object to erase baseflow from above
flowdata = hydro.Discharge(flow.TimeStamp, flow.discharge_cfs, flow.Rainfall_in)
flow['E_bflow'] = flowdata.Eckhardt()
flow['E_bflow'] = flowdata.Eckhardt()
flow['E_bflow'] = flowdata.Eckhardt()
flowdata.plot([flow.TimeStamp, flow.E_bflow], log=False)

# flow DataFrame now has discharge and base flow calculated with two different 
# methods. Plot that DataFrame if you want to compare the methods.

