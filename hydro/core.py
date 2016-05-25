import numpy as np, pandas as pd

def flow_duration(series):
    """ Creates the flow duration curve for a flow dataset """
    fd = pd.Series(series).value_counts()               # frequency of unique values
    fd.sort_index(inplace=True)                         # sort in order of increasing discharges
    fd = fd.cumsum()                                    # cumulative sum of frequencies
    fd = fd.apply(lambda x: 100 - x/fd.max() * 100)     # normalize
    return fd

def Lyne_Hollick(series, alpha=.925, direction='f'):
    """Recursive digital filter for baseflow separation. Based on Lyne and Hollick, 1979.
    series = array of discharge measurements
    alpha = filter parameter
    direction = (f)orward or (r)everse calculation
    """
    series = np.array(series)
    f = np.zeros(len(series))
    if direction == 'f':
        for t in np.arange(1,len(series)):
            f[t] = alpha * f[t-1] + (1 + alpha)/2 * (series[t] - series[t-1])
            if series[t] - f[t] > series[t]:
                f[t] = 0
    elif direction == 'r':
        for t in np.arange(len(series)-2, 1, -1):
            f[t] = alpha * f[t+1] + (1 + alpha)/2 * (series[t] - series[t+1])
            if series[t] - f[t] > series[t]:
                f[t] = 0
    return np.array(series - f)
 
def Eckhardt(series, alpha=.98, BFI=.80):
    """Recursive digital filter for baseflow separation. Based on Eckhardt, 2004.
    series = array of discharge measurements
    alpha = filter parameter
    BFI = BFI_max (maximum baseflow index)
    """
    series = np.array(series)
    f = np.zeros(len(series))
    f[0] = series[0]
    for t in np.arange(1,len(series)):
        f[t] = ((1 - BFI) * alpha * f[t-1] + (1 - alpha) * BFI * series[t]) / (1 - alpha * BFI)
        if f[t] > series[t]:
            f[t] = series[t]
    return f

def sinuosity(Easting, Northing, length, distance):
    """Calculates sinuosity at each data point. Easting and Northing are lat/longs
    projected into measureable units. Length is distance to calculate the
    sinuosity on either side of points. Distance is the stream distance between data points.
    """
    pnts = int(length/distance) # number of points for each reach 
    East = np.array(Easting)
    North = np.array(Northing)
    
    # Calculate sinuosity: stream length / straight line distance
    sin = np.zeros(len(East))
    l = len(East)
    b = pnts * 2 * distance           # calculates stream distance for pnts in middle
    for i in range(int(l)):
        if i < pnts: # first few points
            a = (i + pnts) * distance # calculates stream distance
            sin[i] = a / np.sqrt(np.abs(East[i+pnts] - East[0])**2 
                     + np.abs(North[i+pnts] - North[0])**2)
        elif len(sin)-i < pnts +1: # last few points
            a = (len(sin)-i + pnts) * distance
            sin[i] = a /np.sqrt(np.abs(East[len(sin)-1] - East[i-pnts])**2
                     + np.abs(North[len(sin)-1] - North[i-pnts])**2)
        else: # most points are evaluated here
            sin[i] = b / np.sqrt(np.abs(East[i+pnts] - East[i-pnts])**2
                     + np.abs(North[i+pnts] - North[i-pnts])**2)                           
    return sin
