"""
Define plotting functions that return a Figure (rather than using pyplot)
Used by GUI applications
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from invisible_cities.core.mpl_functions import set_plot_labels
from invisible_cities.core.core_functions import define_window
from invisible_cities.core.core_functions import loc_elem_1d
from invisible_cities.core.system_of_units_c import units


def _hist_outline(dataIn, *args, **kwargs):
    (histIn, binsIn) = np.histogram(dataIn, bins='auto', *args, **kwargs)

    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

def hist_1d(data, xlo, xhi):
    """Returns a Figure corresponding to 1d histogram"""
    (bins, n) = _hist_outline(data)
    ylo = 0
    yhi = max(n) * 1.1

    fig = Figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(bins, n, 'k-')
    ax1.axis([xlo, xhi, ylo, yhi])

    return fig


def fplot_vector(v, figsize=(10,10)):
    """Plot vector v and return figure """
    fig = Figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(v)
    return fig


def plot_vector(v, figsize=(10,10)):
    """Plot vector v and return figure """
    plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(v)


def fplot_xy(x, y, figsize=(10,10)):
    """Plot y vs x and return figure"""
    fig = Figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(x,y)
    return fig


def plot_xy(x, y, figsize=(10,10)):
    """Plot y vs x and return figure"""
    plt.figure(figsize=figsize)
    ax1 = plt.subplot(1, 1, 1)
    plt.plot(x,y)


def fplot_pmt_waveforms(pmtwfdf, zoom=False, window_size=800, figsize=(10,10)):
    """plot PMT wf and return figure"""
    fig = Figure(figsize=figsize)
    for i in range(len(pmtwfdf)):
        first, last = 0, len(pmtwfdf[i])
        if zoom:
            first, last = define_window(pmtwfdf[i], window_size)

        ax = fig.add_subplot(3, 4, i+1)
        set_plot_labels(xlabel="samples", ylabel="adc")
        ax.plot(pmtwfdf[i][first:last])
    return fig


def plot_pmt_waveforms(pmtwfdf, zoom=False, window_size=800, figsize=(10,10)):
    """plot PMT wf and return figure"""
    plt.figure(figsize=figsize)
    for i in range(len(pmtwfdf)):
        first, last = 0, len(pmtwfdf[i])
        if zoom:
            first, last = define_window(pmtwfdf[i], window_size)

        ax = plt.subplot(3, 4, i+1)
        set_plot_labels(xlabel="samples", ylabel="adc")
        plt.plot(pmtwfdf[i][first:last])


def fplot_signal_vs_time_mus(signal,
                            t_min      =    0,
                            t_max      = 1200,
                            signal_min =    0,
                            signal_max =  200,
                            figsize=(10,10)):
    """Plot signal versus time in mus (tmin, tmax in mus). """
    fig = Figure(figsize=figsize)
    tstep = 25 # in ns
    PMTWL = signal.shape[0]
    signal_t = np.arange(0., PMTWL * tstep, tstep)/units.mus
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_xlim([t_min, t_max])
    ax1.set_ylim([signal_min, signal_max])
    set_plot_labels(xlabel = "t (mus)",
                    ylabel = "signal (pes/adc)")
    ax1.plot(signal_t, signal)
    return fig


def plot_signal_vs_time_mus(signal,
                            t_min      =    0,
                            t_max      = 1200,
                            signal_min =    0,
                            signal_max =  200,
                            figsize=(6,6)):
    """Plot signal versus time in mus (tmin, tmax in mus). """
    plt.figure(figsize=figsize)
    tstep = 25 # in ns
    PMTWL = signal.shape[0]
    signal_t = np.arange(0., PMTWL * tstep, tstep)/units.mus
    ax1 = plt.subplot(1, 1, 1)
    ax1.set_xlim([t_min, t_max])
    ax1.set_ylim([signal_min, signal_max])
    set_plot_labels(xlabel = "t (mus)",
                    ylabel = "signal (pes/adc)")
    plt.plot(signal_t, signal)


def fplot_pmt_signals_vs_time_mus(pmt_signals,
                                 pmt_active,
                                 t_min      =    0,
                                 t_max      = 1200,
                                 signal_min =    0,
                                 signal_max =  200,
                                 figsize=(10,10)):
    """Plot PMT signals versus time in mus  and return figure."""

    tstep = 25
    PMTWL = pmt_signals[0].shape[0]
    signal_t = np.arange(0., PMTWL * tstep, tstep)/units.mus
    fig = Figure(figsize=figsize)

    for j, i in enumerate(pmt_active):
        ax1 = fig.add_subplot(3, 4, j+1)
        ax1.set_xlim([t_min, t_max])
        ax1.set_ylim([signal_min, signal_max])
        set_plot_labels(xlabel = "t (mus)",
                        ylabel = "signal (pes/adc)")

        ax1.plot(signal_t, pmt_signals[i])

    return fig


def plot_pmt_signals_vs_time_mus(pmt_signals,
                                 pmt_active,
                                 t_min      =    0,
                                 t_max      = 1200,
                                 signal_min =    0,
                                 signal_max =  200,
                                 figsize=(10,10)):
    """Plot PMT signals versus time in mus  and return figure."""

    tstep = 25
    PMTWL = pmt_signals[0].shape[0]
    signal_t = np.arange(0., PMTWL * tstep, tstep)/units.mus
    plt.figure(figsize=figsize)

    for j, i in enumerate(pmt_active):
        ax1 = plt.subplot(3, 4, j+1)
        ax1.set_xlim([t_min, t_max])
        ax1.set_ylim([signal_min, signal_max])
        set_plot_labels(xlabel = "t (mus)",
                        ylabel = "signal (pes/adc)")

        plt.plot(signal_t, pmt_signals[i])


def fplot_s12(S12, figsize=(10,10)):
    """Plot the peaks of input S12 and return figure..

    S12 is a dictionary
    S12[i] for i in keys() are the S12 peaks
    """
    fig = Figure(figsize=figsize)
    set_plot_labels(xlabel = "t (mus)",
                    ylabel = "S12 (pes)")
    xy = len(S12)
    if xy == 1:
        t = S12[0][0]
        E = S12[0][1]
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.plot(t, E)
    else:
        x = 3
        y = xy/x
        if y % xy != 0:
            y = int(xy/x) + 1
        for i in S12.keys():
            ax1 = fig.add_subplot(x, y, i+1)
            t = S12[i][0]
            E = S12[i][1]
            ax1.plot(t, E)
    return fig


def plot_s12(S12, figsize=(6,6)):
    """Plot the peaks of input S12.

    S12 is a dictionary
    S12[i] for i in keys() are the S12 peaks
    """
    plt.figure(figsize=figsize)

    set_plot_labels(xlabel = "t (mus)",
                    ylabel = "S12 (pes)")
    xy = len(S12)
    if xy == 1:
        t = S12[0][0]
        E = S12[0][1]
        ax1 = plt.subplot(1, 1, 1)
        ax1.plot(t, E)
    else:
        x = 3
        y = xy/x
        if y % xy != 0:
            y = int(xy/x) + 1
        for i in S12.keys():
            ax1 = plt.subplot(x, y, i+1)
            t = S12[i][0]
            E = S12[i][1]
            plt.plot(t, E)
