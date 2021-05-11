import datetime as dt
import numpy as np

def cal_bins(binsize, start_date, end_date):
    """
    Calculates bins of size binsize between start_date and end_date.
    If end_date-start_date is not int multiple of binsize,
        will use a final bin date > end_date
    
    Args:
        binsize: integer bin size in days
        start_date: datetime object indicating the first day
        end_date: datetime object indicating second day

    Returns:
        date_bins: list of datetime objects representing bounds
    """

    if start_date >= end_date:
        raise ValueError("Start date must be before end date.")
    assert type(binsize) == int, "binsize must be integer days"

    date_bins = [start_date]
    while date_bins[-1] < end_date:
        date_bins.append(date_bins[-1] + dt.timedelta(days=binsize))
    return date_bins

def find_bin(d, date_bins):
    """
    Find which bin d belongs to, where bins are right boundaries.
    The first bin is designed to be (-\infty, d0].
    The second bin is (d0, d1], etc.
    The final bin is (dn-1, dn]

    Args:
        d: value of which to find a bin
        date_bins: a list of ordered values for which < is defined
    """
    for j, bd in enumerate(date_bins):
        if d <= bd:
            return j
    raise Exception("Value %s does not fit into any bins." % str(d))

def histogram(dates, date_bins):
    """
    Place dates into date_bins

    Args:
        dates: list of datetime objects
        date_bins: a list of n datetime objects representing 
            right inclusive boundaries of n bins

    Returns:
        counts: np.array of integer counts
    """
    counts = np.zeros(len(date_bins))
    for d in dates:
        j = find_bin(d, date_bins)
        counts[j] += 1
    return counts

