"""Early-warning signals: does a reef 'wobble' more and recover more slowly
as it approaches a tipping point?"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import kendalltau


def wobble(series, smooth_years=10):
    """Remove the slow trend, keep the year-to-year wobble (residuals)."""
    return series - gaussian_filter1d(series, smooth_years)


def rolling_indicators(series, window):
    """Variance and lag-1 autocorrelation of the wobble in a moving window.

    variance        = how WIDE the wobbles are
    autocorrelation = how much this year looks like last year
                      (high = the reef is slow to bounce back)
    Value at index i uses years i-window .. i-1.
    """
    r = wobble(series)
    var = np.full(len(r), np.nan)
    ac = np.full(len(r), np.nan)
    for i in range(window, len(r) + 1):
        w = r[i - window:i]
        var[i - 1] = np.var(w)
        ac[i - 1] = np.corrcoef(w[:-1], w[1:])[0, 1]
    return var, ac


def trend_strength(values):
    """Kendall's tau: +1 = always rising, 0 = no trend, -1 = always falling."""
    ok = ~np.isnan(values)
    return kendalltau(np.arange(ok.sum()), values[ok])[0]
