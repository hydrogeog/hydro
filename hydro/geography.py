#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np, matplotlib.pyplot as plt

def sinuosity(Easting, Northing, length, distance):
    """Calculates sinuosity at each data point. Easting and Northing are lat/longs
    projected into measureable units. Length is distance to calculate the
    sinuosity on either side of points. Distance is the stream distance between data points.

    To calculate sinuosity for an entire stream, use the start and end points in
    'Easting' and 'Northing', set length to 2 * stream length, and distance to
    stream length.
    """
    pnts = int(length/distance) # number of points for each reach
    East = np.array(Easting)
    North = np.array(Northing)

    if pnts / 2 == 1:  # if the calculation is only for two points
        return distance / np.sqrt(np.abs(East[0] - East[1])**2
                 + np.abs(North[0] - North[1])**2)
    else:

        # Calculate sinuosity: stream length / straight line distance
        sin = np.zeros(len(East))
        l = len(East)
        b = pnts * 2 * distance  # calculates stream distance for pnts in middle of dataset
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

def Profile_smoothing(elevation, distance=None, plot=False):
    """Removes the 'bumps' present in an elevation profile caused by roads &
     imperfections in DEMs. Data must be arranged from highest elevation to lowest.
     
     Set plot=True and define distance to graph the resulting profile
    """
    elevation = np.array(elevation)
    output_elevation = np.zeros(len(elevation)); output_elevation[0]=elevation[0]
    elevation[-1] = min(elevation)
    i=1
    while i < len(elevation):                    # loops through elevation dataset
        if elevation[i] > elevation[i-1]:             # if the current value is greater than the previous,
            j=i
            while j < len(elevation)-1:          # loop through the remaining data...
                j = j+1
                if elevation[j] <= elevation[i-1]:    # to find the next point that is lower or equal to the last datum.
                    for indx,k in enumerate(range(i, j)):       # changes the incorrect data...
                        adj = np.linspace(elevation[i-1], elevation[j], num=len(elevation[i:j]), endpoint=False)   # to a linear interp.
                        output_elevation[k] = adj[indx]
                    break       # exits if statement
            i=j-1               # starting point for the next iteration
        else:
            output_elevation[i] = elevation[i]   # if the elevation is less than or equal to the previous, that point kept
        i = i+1                     # keeps the while loop interating to next datum
        
    if plot:
            fig = plt.figure(figsize=(20,3))
            ax1 = fig.add_subplot(111)
            ax1.plot(distance, output_elevation, label='Elevation Profile')
            ax1.invert_xaxis()
            ax1.set_ylabel('Elevation')                    # y label
            ax1.set_xlabel('Distance from mouth')          # x label
            plt.title("Longitudinal Profile")       # title
            plt.show()

    return output_elevation

def IDW(x, y, z, xi, yi, power=1):
    '''
    Inverse distance weighted interpolation
    
    x, y, and z are the data arrays\n
    xi and yi are the x and y grid over which to calculate IDW\n
    power is an int which is multiplied by the weight given to each grid point.
    A low power leads to a greater weight towards a grid point value of rainfall 
    from remote rain gauges. 
    '''
    xi, yi = xi.flatten(), yi.flatten()
    dist = distance_matrix(x, y, xi, yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / (dist)**power

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point by all observed Z-values
    zi = np.dot(weights.T, z)
    s = int(np.sqrt(len(zi)))
    zi = zi.reshape((s, s))
    return zi

def distance_matrix(x0, y0, x1, y1):
    '''Distance matrix for IDW calculation'''
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])
    return np.hypot(d0, d1)
