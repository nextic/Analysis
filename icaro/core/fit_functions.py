import numpy as np

import invisible_cities.core .fit_functions  as     fitf
from   invisible_cities.core .core_functions import in_range
from   invisible_cities.evm  .ic_containers  import Measurement
from   invisible_cities.icaro.hst_functions  import shift_to_bin_centers
from   invisible_cities.icaro.hst_functions  import labels


def conditional_labels(with_title=True):
    """
    Wrapper around the labels function that allows removing
    the title by changing one parameter.
    """
    with_    = lambda *args: labels(*args)
    without_ = lambda *args: labels(*args[:2], "")
    return with_ if with_title else without_


def gauss_seed(x, y, sigma_rel=0.05):
    """
    Estimate the seed for a gaussian fit to the input data.
    """
    y_max  = np.argmax(y) # highest bin
    x_max  = x[y_max]
    sigma  = sigma_rel * x_max
    amp    = y_max * (2 * np.pi)**0.5 * sigma * np.diff(x)[0]
    seed   = amp, x_max, sigma
    return seed


def expo_seed(x, y, eps=1e-12):
    """
    Estimate the seed for a exponential fit to the input data.
    """
    x, y  = zip(*sorted(zip(x, y)))
    const = y[0]
    slope = (x[-1] - x[0]) / np.log(y[-1] / (y[0] + eps))
    seed  = const, slope
    return seed


def relative_errors(values, errors, default=0, percentual=False):
    """
    Compute relative errors from input values with safety checks.
    If the relative error cannot be computed, a default value is
    used. The errors can be in percent if the `percentual`
    argument is True.
    """
    ok         = values != 0
    scale      = 100 if percentual else 1
    rel_e      = np.empty_like(values)
    rel_e[ ok] = errors[ok] / np.abs(values[ok]) * scale
    rel_e[~ok] = default
    return rel_e


def to_relative(data, *args, **kwargs):
    """
    Produce another Measurement instance with relative instead of
    the absolute ones.
    """
    return Measurement(data.value, relative_errors(*data, *args, **kwargs))


def quick_gauss_fit(data, bins):
    """
    Histogram input data and fit it to a gaussian with the parameters
    automatically estimated.
    """
    y, x  = np.histogram(data, bins)
    x     = shift_to_bin_centers(x)
    seed  = gauss_seed(x, y)
    f     = fitf.fit(fitf.gauss, x, y, seed)
    assert np.all(f.values != seed)
    return f


def fit_profile_1d_expo(xdata, ydata, nbins, *args, **kwargs):
    """
    Make a profile of the input data and fit it to an exponential
    function with the parameters automatically estimated.
    """
    x, y, yu     = fitf.profileX(xdata, ydata, nbins, *args, **kwargs)
    valid_points = yu > 0

    x    = x [valid_points]
    y    = y [valid_points]
    yu   = yu[valid_points]
    seed = expo_seed(x, y)
    f    = fitf.fit(fitf.expo, x, y, seed, sigma=yu)
    assert np.all(f.values != seed)
    return f


def fit_slices_1d_gauss(xdata, ydata, xbins, ybins, min_entries=1e2):
    """
    Slice the data in x, histogram each slice, fit it to a gaussian
    and return the relevant values.

    Parameters
    ----------
    xdata, ydata: array_likes
        Values of each coordinate.
    xbins: array_like
        The bins in the x coordinate.
    ybins: array_like
        The bins in the y coordinate for histograming the data.
    min_entries: int (optional)
        Minimum amount of entries to perform the fit.

    Returns
    -------
    mean: Measurement(np.ndarray, np.ndarray)
        Values of mean with errors.
    sigma: Measurement(np.ndarray, np.ndarray)
        Values of sigma with errors.
    chi2: np.ndarray
        Chi2 from each fit.
    valid: boolean np.ndarray
        Where the fit has been succesfull.
    """
    nbins  = np.size (xbins) - 1
    mean   = np.zeros(nbins)
    sigma  = np.zeros(nbins)
    meanu  = np.zeros(nbins)
    sigmau = np.zeros(nbins)
    chi2   = np.zeros(nbins)
    valid  = np.zeros(nbins, dtype=bool)

    for i in range(nbins):
        sel = in_range(xdata, *xbins[i:i + 2])
        if np.count_nonzero(sel) < min_entries: continue

        try:
            f = quick_gauss_fit(ydata[sel], ybins)
            mean  [i] = f.values[1]
            meanu [i] = f.errors[1]
            sigma [i] = f.values[2]
            sigmau[i] = f.errors[2]
            chi2  [i] = f.chi2
            valid [i] = True
        except:
            pass
    return Measurement(mean, meanu), Measurement(sigma, sigmau), chi2, valid


def fit_slices_2d_gauss(xdata, ydata, zdata, xbins, ybins, zbins, min_entries=1e2):
    """
    Slice the data in x and y, histogram each slice, fit it to a gaussian
    and return the relevant values.

    Parameters
    ----------
    xdata, ydata, zdata: array_likes
        Values of each coordinate.
    xbins, ybins: array_likes
        The bins in the x and y coordinates.
    zbins: array_like
        The bins in the z coordinate for histograming the data.
    min_entries: int (optional)
        Minimum amount of entries to perform the fit.

    Returns
    -------
    mean: Measurement(np.ndarray, np.ndarray)
        Values of mean with errors.
    sigma: Measurement(np.ndarray, np.ndarray)
        Values of sigma with errors.
    chi2: np.ndarray
        Chi2 from each fit.
    valid: boolean np.ndarray
        Where the fit has been succesfull.
    """
    nbins_x = np.size (xbins) - 1
    nbins_y = np.size (ybins) - 1
    nbins   = nbins_x, nbins_y
    mean    = np.zeros(nbins)
    sigma   = np.zeros(nbins)
    meanu   = np.zeros(nbins)
    sigmau  = np.zeros(nbins)
    chi2    = np.zeros(nbins)
    valid   = np.zeros(nbins, dtype=bool)

    for i in range(nbins_x):
        sel_x = in_range(xdata, *xbins[i:i + 2])
        for j in range(nbins_y):
            sel_y = in_range(ydata, *ybins[i:i + 2])
            sel   = sel_x & sel_y
            if np.count_nonzero(sel) < min_entries: continue

            try:
                f = quick_gauss_fit(zdata[sel], zbins)
                mean  [i, j] = f.values[1]
                meanu [i, j] = f.errors[1]
                sigma [i, j] = f.values[2]
                sigmau[i, j] = f.errors[2]
                chi2  [i, j] = f.chi2
                valid [i, j] = True
            except:
                pass
    return Measurement(mean, meanu), Measurement(sigma, sigmau), chi2, valid


def fit_slices_2d_expo(xdata, ydata, zdata, tdata,
                       xbins, ybins, nbins_z, zrange=None,
                       min_entries = 1e2):
    """
    Slice the data in x and y, make the profile in z of t,
    fit it to a exponential and return the relevant values.

    Parameters
    ----------
    xdata, ydata, zdata, tdata: array_likes
        Values of each coordinate.
    xbins, ybins: array_like
        The bins in the x coordinate.
    nbins_z: int
        The number of bins in the z coordinate for the profile.
    zrange: length-2 tuple (optional)
        Fix the range in z. Default is computed from min and max
        of the input data.
    min_entries: int (optional)
        Minimum amount of entries to perform the fit.

    Returns
    -------
    const: Measurement(np.ndarray, np.ndarray)
        Values of const with errors.
    slope: Measurement(np.ndarray, np.ndarray)
        Values of slope with errors.
    chi2: np.ndarray
        Chi2 from each fit.
    valid: boolean np.ndarray
        Where the fit has been succesfull.
    """
    nbins_x = np.size (xbins) - 1
    nbins_y = np.size (ybins) - 1
    nbins   = nbins_x, nbins_y
    const   = np.zeros(nbins)
    slope   = np.zeros(nbins)
    constu  = np.zeros(nbins)
    slopeu  = np.zeros(nbins)
    chi2    = np.zeros(nbins)
    valid   = np.zeros(nbins, dtype=bool)

    if zrange is None:
        zrange = np.min(zdata), np.max(zdata)
    for i in range(nbins_x):
        sel_x = in_range(xdata, *xbins[i:i + 2])
        for j in range(nbins_y):
            sel_y = in_range(ydata, *ybins[j:j + 2])
            sel   = sel_x & sel_y
            if np.count_nonzero(sel) < min_entries: continue

            try:
                f = fit_profile_1d_expo(zdata[sel], tdata[sel], nbins_z, xrange=zrange)
                const [i, j] = f.values[0]
                constu[i, j] = f.errors[0]
                slope [i, j] = f.values[1]
                slopeu[i, j] = f.errors[1]
                chi2  [i, j] = f.chi2
                valid [i, j] = True
            except:
                pass
    return Measurement(const, constu), Measurement(slope, slopeu), chi2, valid
