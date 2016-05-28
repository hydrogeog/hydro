import matplotlib.pyplot as plt, numpy as np, pandas as pd
from matplotlib.ticker import FuncFormatter  # used in formatting log scales
import mpl_toolkits.basemap.pyproj as pyproj
import hydro
#%matplotlib inline

data = pd.read_csv("stream.csv")

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
print("Sinuosity of entire stream: ", sin_total)

# Load flow data
flow = pd.read_csv("flow.csv")
flow.TimeStamp = pd.to_datetime(flow.TimeStamp)

# Compute rating curve from given flows
stage = np.array([0.96, 1.02, 1.218, 1.313, 1.186, 1.215,
                  2.121, 1.638, 2.163, 1.902, 2.008])
discharge = np.array([0.43325, 0.5345, 0.4995, 1.26475, 0.5385, 0.9645,
                  25.391, 8.6795, 30.348, 17.604, 22.2395])
coeff, r_2 = hydro.ratingCurve(discharge, stage)

# Plot rating curve
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1, axisbg='white')
plt.grid(True, which='both', color='k', ls='-', alpha=0.2)
ax1.scatter(stage, discharge, marker='^', color='k')
ax1.set_ylabel('Discharge, cfs')
ax1.set_xlabel('Stage, ft')
ax1.set_yscale('log'); ax1.set_xscale('log')      # log scale x and y
ax1.yaxis.set_major_formatter(FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
ax1.xaxis.set_major_formatter(FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
plt.title('Rating Curve')
ax1.set_axisbelow(True)         # puts grid below plot
ax1.text(0.05, 0.9, 'y = ' + str(coeff[0])[0:5] + 'x^' + str(coeff[1])[0:5],
             transform=ax1.transAxes)
ax1.text(0.05, 0.85, 'r^2 = ' + str(r_2)[0:5],
             transform=ax1.transAxes)
line = np.linspace(min(stage), max(stage), 100)
ax1.plot(line, hydro.exp_curve(line, coeff[0], coeff[1]), color='#4a2628')
plt.show()

# Calculate discharge from rating curve
flow['Q'] = hydro.exp_curve(flow.Level_ft, coeff[0], coeff[1])

# Daily mean discharge
dailyQ, day = hydro.dailyMean(flow.Q, flow.TimeStamp, 15)

# Richards-Baker Flashiness
print("R-B Flashiness: ", hydro.RB_Flashiness(dailyQ))

# Baseflow separation
# 3 forward passes with the digital filter
flow['bflow'] = hydro.Lyne_Hollick(flow.Q, .925, 'f')
flow.bflow = hydro.Lyne_Hollick(flow.bflow, .925, 'f')
flow.bflow = hydro.Lyne_Hollick(flow.bflow, .925, 'f')

# plot baseflow
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1, axisbg='white')
plt.title("Lyne and Hollick Filter")
plt.ylabel("Discharge (cfs)")
plt.fill_between(flow.TimeStamp.values, flow.bflow, flow.Q, color='lightblue', label="Discharge")
plt.fill_between(flow.TimeStamp.values, 0, flow.bflow, color='steelblue', label="Baseflow")
plt.legend()
# rainfall
ax2 = ax1.twinx()
ax2.step(flow.Timestamp, flow.Rainfall_in, color='b', label='Rainfall', where='post')
ax2.set_ylim(0,15)
ax2.invert_yaxis()
plt.show()

# Flow duration curve
import probscale # flow duration curves use a probability scale for the x axis
duration = hydro.flow_duration(flow.Q)
durationb = hydro.flow_duration(flow.bflow)
durationc = hydro.flow_duration(flow.Q - flow.bflow)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1, axisbg='white')
ax1.plot(duration, duration.index, 'x', ls='', color='k', label='Total Flow')
ax1.plot(durationb, durationb.index, '.', markeredgecolor='none', ls='', color='steelblue', label='Baseflow')
ax1.plot(durationc, durationc.index, '.', markeredgecolor='none', ls='', color='lightblue', label='Runoff')

# set y axis to log scale and x axis to probability scale
ax1.set_yscale('log')
ax1.set_xscale('prob')
plt.xticks([.01,.1,.5,1,2,5,10,20,30,40,50,60,70,80,90,95,98,99,99.5,99.9,99.99],
           rotation='vertical')
plt.legend()
plt.title('Flow Duration Curve')
plt.ylabel('Flow (cfs)')
plt.xlabel('Percentage of time flow was equaled or exceeded')
plt.show()
