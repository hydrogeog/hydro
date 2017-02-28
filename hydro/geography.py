#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

def IDW(x, y, z, xi, yi):
    '''
    Inverse distance weighted interpolation
    
    x, y, and z are the data arrays\n
    xi and yi are the x and y grid over which to calculate IDW
    '''
    xi, yi = xi.flatten(), yi.flatten()
    dist = distance_matrix(x, y, xi, yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / dist

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
