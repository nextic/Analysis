import numpy as np
from   invisible_cities.core.core_functions import loc_elem_1d
from   invisible_cities.core.system_of_units_c import units
from   collections import namedtuple

KrCoordinates       = namedtuple('KrCoordinates',
                      'Q x y z r phi')

S12Features         = namedtuple('S12Features',
                       'tmin tmax tpeak etot epeak width')

def print_s12(S12):
    """Print peaks of input S12.

    S12 is a dictionary
    S12[i] for i in keys() are the S12 peaks
    """
    print('number of peaks = {}'.format(len(S12)))
    for i in S12:
        print('S12 number = {}, samples = {} sum in pes ={}'
              .format(i, len(S12[i][0]), np.sum(S12[i][1])))
        print('time vector (mus) = {}'.format(S12[i][0]/units.mus))
        print('energy vector (pes) = {}'.format(S12[i][1]/units.pes))


def print_s2si(S2Si):
    """Scan the S2Si objects."""
    for peak, sipm_set in S2Si.items():
        print('S2Si for peak number = {}'.format(peak))
        for sipm, e_array in sipm_set.items():
            print('sipm number = {}, energy = {}'.format(sipm,
                                                         np.sum(e_array)))


def compare_S1(S1, PMT_S1, peak=0, tol=0.25*units.mus):
    """Compare sum S1 with S1 in individual PMT

    input:
    S1 computed with the sum
    PMT_S1 computed with individual PMTs.
    tol is the matching tolerance.

    Return number of matches

    """
    n_match_s1 = 0
    t = S1[peak][0]
    E = S1[peak][1]
    for pmt in PMT_S1:
        if len (PMT_S1[pmt]) > 0:
            for peak, (t2,E2) in PMT_S1[pmt].items():
                diff = abs(t2[0] - t[0])
                if diff < tol:
                    n_match_s1 +=1
                    break  # if one peak is matched look no further
    return n_match_s1


def s12_features(S12, peak_number=0):
    """Return a named tuple of S12 features."""

    t = S12[peak_number][0]
    E = S12[peak_number][1]

    epeak = np.max(E)
    etot  = np.sum(E)
    i_t = loc_elem_1d(E, epeak)
    tpeak = t[i_t]
    width = t[-1] - t[0]

    return S12Features(tmin  = t[0],
                       tmax  = t[-1],
                       tpeak = tpeak,
                       etot  = etot,
                       epeak = epeak,
                       width = width)


def charge_and_position_from_S2Si(S2Si, xs, ys, peak_number=0):
    """
    Charge and position from S2Si
    """
    s2si = S2Si[peak_number]
    xsipm = []
    ysipm = []
    Q = []
    for key, value in s2si.items():
        xsipm.append(xs[key])
        ysipm.append(ys[key])
        Q.append(np.sum(value))
    return np.array(xsipm), np.array(ysipm), np.array(Q)


def drift_time(s1f_peak, s2f, peak_number):
    """
    drift time from EventPmaps
    """

    ts1 = s1f_peak.tpeak
    ts2 = s2f[peak_number].tpeak

    return ts2 -ts1

def kr_coordinates(xsipm, ysipm, Q, drift_time, drift_velocity):
    """
    find kr coordinates
    """
    x    = np.average(xsipm, weights=Q) if np.any(Q) else -9999
    y    = np.average(ysipm, weights=Q) if np.any(Q) else -9999
    r    = (x**2 + y**2)**0.5 if np.any(Q) else -9999
    phi   = np.arctan2(y, x) if np.any(Q) else -9999
    z    = drift_time * drift_velocity
    #z    = drift_time
    Q_tp = np.sum(Q)
    return KrCoordinates(Q=Q_tp, x=x, y=y, z=z, r=r, phi=phi)


def integrate_charge(d):
    """
    Integrate charge from a SiPM dictionary.
    """
    newd = dict((key, np.sum(value)) for key, value in d.items())
    return map(np.array, list(zip(*newd.items())))
