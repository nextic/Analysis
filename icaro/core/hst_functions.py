#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 20:21:59 2017

Decorate histogram function in IC hist_functions module

add statitics into the plots
define canvas to subdivide figure in subplots

@author: hernando
"""

import numpy             as np
import scipy.stats       as stats
import matplotlib.pyplot as plt
import invisible_cities.icaro.hst_functions as hst
from invisible_cities.icaro.hst_functions import shift_to_bin_centers
from invisible_cities.icaro.hst_functions import resolution, gausstext
from invisible_cities.icaro.hst_functions import labels
from invisible_cities.icaro.hst_functions import save_to_folder

class Canvas:
    """ class object to store subplot in figure
    """
    def __init__(self, nx, ny, nxsize=5.4, nysize=6.2):
        """ divide the figure in (nx, ny) subplots
        inputs:
            nx, ny : int, int
                     number of subplots in (x, y)
            nxsize, nysize: float, flot
                    x, y size of the subplots
        """
        self.nx = nx
        self.ny = ny
        self.n  = 1
        plt.figure(figsize=(nysize*ny, nxsize*nx))
        plt.subplot(nx, ny, self.n)

    def __call__(self, i):
        """ locate in subplot i in figure
        """
        plt.subplot(self.nx, self.ny, i)
        return True

    def __add__(self, i):
        """ locate to next i plot in figure
        """
        self.n += i
        plt.subplot(self.nx, self.ny, n)
        return True


def _xypos(xs, ys, xf=0.1, yf=0.7):
    x0, dx = min(xs), max(xs)-min(xs)
    y0, dy = min(ys), max(ys)-min(ys)
    xp = x0 + xf*dx
    yp = y0 + yf*dy
    return (xp, yp)

def _hist(xs, bins=100, range=None, stats=('entries', 'mean', 'rms'),
         xylabels = (), stats_xypos=(0.1, 0.7),
         *args, **kargs):
    """
    hist function extended to plot stattistics
    inputs:
        stats: (string,) or None. Default: ('entries', 'mean', 'rms')
            name of the statistics parameters to be into the histogram
            if none use None
            possibilities: 'total entries'
        stats_xypos: (float, float) (default (0.1, 0.7))
            (x, y) relative positions in the plot,
            x, y values in rage (0., 1.)
    outputs:
        none
    """
    if (range==None):
        range = (np.min(xs), np.max(xs))
    cc = hst.hist(xs, bins=bins, range=range, *args, **kargs);
    if (not stats):
        return cc
    ys, xedges = np.histogram(xs, bins, range=range)
    ns = len(xs)
    sel = np.logical_and(xs >= range[0], xs <= range[1])
    nos, mean, rms = len(xs[sel]), np.mean(xs[sel]), np.std(xs[sel])
    epsilon = (1.*nos)/(1.*ns)
    ss = ''
    if ('total entries') in stats:
        ss += 'total entries  {0:d} \n'.format(ns)
    if ('entries') in stats:
        ss += 'entries {0:d} \n'.format(nos)
    if ('mean') in stats:
        ss += 'mean {0:.3f} \n'.format(mean)
    if ('rms') in stats:
        ss += 'rms  {0:.3f} \n'.format(rms)
    xp, yp = _xypos(xedges, ys, xf=stats_xypos[0], yf=stats_xypos[1])
    plt.text(xp, yp, ss)
    return cc


def _decorate(fun):
    def _decorator(*args, **kargs):
        if ('canvas' in kargs):
            kargs['new_figure'] = False
            del kargs['canvas']
        has_xylabels, xylabels = 'xylabels' in kargs, ()
        if (has_xylabels):
            xylabels = kargs['xylabels']
            del kargs['xylabels']
        response = fun(*args, **kargs)
        if (len(xylabels)>0):
            hst.labels(*xylabels)
        return response
    return _decorator

hist                = _decorate(_hist)
plot                = _decorate(hst.plot)
hist2d              = _decorate(hst.hist2d)
errorbar            = _decorate(hst.errorbar)
pdf                 = _decorate(hst.pdf)
scatter             = _decorate(hst.scatter)
profile_and_scatter = _decorate(hst.profile_and_scatter)
hist2d_profile      = _decorate(hst.hist2d_profile)
display_matrix      = _decorate(hst.display_matrix)
