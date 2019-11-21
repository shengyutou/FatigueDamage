# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:32:56 2019

@author: WangPeide
"""

# rainflow.py
"""
Rainflow cycle counting function based on Downing's Method 1 of the 
paper "Simple rainflow counting algorithms" by S.D.Downing and D.F.Socie
International Journal of Fatigue, January 1982.
This method is for limited histories where the data has been obtained 
and stored. Method 2 from the same paper (not implemented here) is 
designed for open-ended histories, typically in a monitoring situation.
Method 1 requires that the signal be first re-arranged so that signal 
starts and finishes with the largest (magnitude) peak or trough (we just
say 'peak' to include both from now on). This is done in the function
_get_peaks
"""
from numpy import array, roll, concatenate, argmax, abs
from numba import jit
from collections import deque, defaultdict

def rainflow_count(data, ndigits=2):
    """
    Return ranges, means of a 1D array using Rainflow cycle counting
    :Parameters:
        raw_data: 1-D ndarray. The data to be analysed. Int or float.
    :Returns:
        ranges, means: 1D ndarrays. The cycle ranges and mean values
    :References:
        Rainflow cycle counting method based on Downing's Method 1 of 
        the paper "Simple rainflow counting algorithms" by S.D.Downing
        and D.F.Socie, International Journal of Fatigue, January 1982.
    """
    peaks = _get_peaks(data)

    if len(peaks) == 0:
        return [], []

    ranges, means = _ranges_means(peaks)
    counts = defaultdict(float)
    for low, high, times in _ranges_means(data):
        delta = round_(abs(high - low) / 2)
        mean = round_((high - low) / 2)
        counts[(delta, mean)] += times
    return sorted(counts.items())

    return array(ranges), array(means)


def _get_peaks(data, ndigits=2):
    """
    Return peaks & troughs of data, starting & finishing at the maximum
    """
    data = _round(data, ndigits)
    # eliminate repeated data points:
    data = data[data != roll(data, 1)]

    if len(data) == 0:
        return []

    # split and rejoin at largest abs value
    max_idx = argmax(abs(data))
    data = concatenate((data[max_idx:], data[:max_idx]))

    # find peaks and troughs
    prv = roll(data, -1)
    nxt = roll(data, 1)

    isPeak = ((data > prv) & (data > nxt)) | ((data < prv) & (data < nxt))

    # Close off the signal with the max (ie first) value, as required.
    return concatenate((data[isPeak], data[:1]))  


@jit(nopython=True)  # gives ~20x speed improvement
def _ranges_means(peaks):
    """
    Return ranges, means of cycles counted using Downing's method 1.
    """
    values = []  # of the current peaks being processed
    ranges = []
    means = []

    for peak in peaks:
        values.append(peak)
        while len(values) > 2:
            X = abs(values[-1] - values[-2])
            Y = abs(values[-2] - values[-3])
            if X < Y:
                break
            ranges.append(Y)
            means.append(0.5*(values[-2] + values[-3]))
            values[-3] = values[-1]
            values.pop()
            values.pop()

    return ranges, means

@jit(nopython=True) 
def _round(x, ndigits=2):
    return array([round(i,ndigits) for i in x])