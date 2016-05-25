import numpy as np

def sinuosity(Easting, Northing, length, distance):
    """ Calculates sinuosity at each data point. Easting and Northing are lat/longs
    projected into measureable units. Length is distance to calculate the
    sinuosity on either side of points. Distance is the distance between data points.
    """
    pnts = int(length/distance) # number of points for each reach 
    East = np.array(Easting)
    North = np.array(Northing)
    
    # Calculate sinuosity
    sin = np.zeros(len(East))
    l = len(East)
    b = pnts * 2 * distance           # calculates straight line distance for pnts in middle
    for i in range(int(l)):
        if i < pnts: # first few points
            a = (i + pnts) * distance # calculates straight line distance
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
