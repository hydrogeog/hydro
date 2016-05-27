import matplotlib.pyplot as plt, numpy as np, pandas as pd
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

# Calculate sinuosity
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

sin_total = hydro.sinuosity([data.Easting.iloc[0], data.Easting.iloc[-1]],
                            [data.Northing.iloc[0], data.Northing.iloc[-1]],
                            data['cng_(Feet)'].max()*2, data['cng_(Feet)'].max())
print("Sinuosity of entire stream: ", sin_total)
