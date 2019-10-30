#!/usr/bin/env python
# _*_coding:utf-8 _*_
"""
@Time    : 10/30/19 1:35 PM
@Author  : Louis Wang
@Email   : louis1990.wang@outlook.com
@FileName: Fatigue_LouisW.py
@Software: PyCharm
"""


class SnGL(object):
    def __init__(self, mat='GGG', t=120, j=3, gamma_m=1.265, modified_t=False):
        """
        @description:
            Calculation of the basic SN curve parameters according to GL2010.
        @para:
            Rm: {float} Strength of tension. /MPa
            self.Rp: {float}strength of yield. /MPa
            mat: {str} material type. GGG for cast iron.
            t: {float} thickness of the material.
            sign_tc: {bool} whether to consider thickness correction.
            Rz: {float} Surface roughness.
            R: {float} stress ratio.
            j: {float} quality level for component.
            j_0: {float} constant for material and test method.
            S_pu: {float} survival probability.
            gamma_m: {float} partial safety factor for material.
            sigma_b: {float} Modified tensile strength. /MPa
            M: {float} Mean stress sensitivity
        @return:
        """
        self.Rm = 360
        self.Rp = 220
        self.mat = mat
        self.t = t
        if modified_t:
            self.sign_tc = 1
        else:
            self.sign_tc = 0
        self.Rz = 125
        self.R = -1
        self.j = j
        self.j_0 = 1
        self.S_pu = 2 / 3
        self.gamma_M = gamma_m
        self.sigma_b = self.Rm * 1.06
        if self.mat == 'GGG':
            self.M = 0.00035 * self.sigma_b + 0.08
        else:
            self.M = 0.00035 * self.sigma_b + 0.05

    def cal(self, plot=False, save=True):
        """
            Calculation the basic parameters.
        :return: m1, m2, sigma_D, N_d, plot(if plot == True)
        """
        from math import log10, sqrt
        import matplotlib.pyplot as plt
        # Surface roughness factor
        F_o = (1 - 0.22 * (log10(self.Rz)**0.64) * log10(self.sigma_b) + 0.45 *
               (log10(self.Rz)**0.53))
        alpha_k = 1
        n = 1
        belt_k = alpha_k / n
        F_ok = sqrt(belt_k**2 - 1 + 1 / (F_o**2))

        # Fatigue strength of specimen
        if self.mat == 'GGG':
            sigma_w = 0.27 * self.sigma_b + 100
        else:
            sigma_w = 0.27 * self.sigma_b + 85

        # Fatigue strength of component
        sigma_wk = sigma_w / F_ok

        # Slopes of SN curve m1 and m2
        self.m1 = 5.5 / (F_ok**2) + 6
        self.m2 = 2 * self.m1 - 1

        # Factor for influence of mean stress
        u = 1 / (self.M + 1) * sigma_wk / self.sigma_b
        a = (1 + self.R) / (1 - self.R) * sigma_wk / self.sigma_b
        p = (1 / (self.M + 1) - 1 + u**2) / (u**2 - u)
        if a == 0:
            Fm = 1
        else:
            if p <= 1:
                Fm = -1 * (1 + p * a) / (2 * a**2 * (1 - p)) + \
                     sqrt(1 / (1 - p) / a**2 + ((1 + p * a) / 2 / a**2 / (1 - p))**2)
            else:
                Fm = -1 * (1 + p * a) / (2 * a**2 *(1 - p)) - \
                     sqrt(1 / (1 - p) / a**2 + ((1 + p * a) / 2 / a**2 / (1 - p))**2)

        # Stress amplitude at knee of SN curve
        sigma_A = sigma_wk * Fm
        # Number of load cycles at knee of SN curve
        self.N_D = 10**(6.8 - 3.6 * (1 / self.m1))
        # Upgrading factors
        S_d = 0.85**(self.j - self.j_0)
        S_t = (self.t / 25)**((-0.15) * self.sign_tc)
        S = self.S_pu * S_d * S_t

        # Upgraded stress amplitude at knee of SN curve
        self.sigma_D = sigma_A * S / self.gamma_M
        # Upper limit of fatigue life line
        sigma_1 = self.Rp * (1 - self.R) / self.gamma_M
        # Number of load cycles at upper fatigue limit
        N_1 = self.N_D * (2 * self.sigma_D / sigma_1)**self.m1
        N_e = 10**9
        sigma_e = (self.N_D / N_e)**(1 / self.m2) * self.sigma_D
        # Parameters of the haigh diagram
        sigma_p = self.sigma_D / (self.M + 1)
        # Parameters saved to .txt file
        if save:
            with open('SN_Curve.txt', 'w') as f:
                f.write('Summary of the SN curve Parameter\n')
                f.write('*' * 30 + '\n')
                f.write('Material:\n%s\n' % self.mat)
                f.write('Stress ratio R:\n%s\n' % self.R)
                f.write('*' * 30 + '\n')
                f.write(
                    'Calculation of SN curve according to the GL2010 5.B.3.1 and 5.B.3.2.\n'
                )
                f.write('*' * 30 + '\n')
                f.write('Tensile strength:\n')
                f.write('%s\n' % self.Rm)
                f.write('Yield Strength:\n')
                f.write('%s\n' % self.Rp)
                f.write('Alternating stress limit (Amplitude/Range):\n')
                f.write('%.2f/%.2f\n' % (self.sigma_D, 2 * self.sigma_D))
                f.write('Pulsating stress limit (Amplitude/Range):\n')
                f.write('%.2f/%.2f\n' % (sigma_p, 2 * sigma_p))
                f.write(
                    'Slope of the SN curve (Left/Right of the knee point):\n')
                f.write('%.2f/%.2f\n' % (self.m1, self.m2))
                f.write('Number of load cycles at knee of SN curve:\n')
                f.write('%.2e\n' % self.N_D)
                f.write('Mean stress sensitivity:\n')
                f.write('%.3f\n' % self.M)
                f.write('*' * 30 + '\n')
                f.write('Parameter of the Haigh Diagram:\n')
                f.write('*' * 30 + '\n')
                f.write('Mean stress / Amplitude:\n')
                f.write('%11.2f/%11d\n' % ((self.Rp / self.gamma_M), 0))
                if self.mat == 'GGG':
                    f.write('%11.2f/%11.2f\n' %
                            ((self.Rp / self.gamma_M - self.sigma_D) /
                             (1 - self.M), self.M *
                             (self.Rp / self.gamma_M - self.sigma_D) /
                             (self.M - 1) + self.sigma_D))
                else:
                    f.write("")
                f.write('%11.2f/%11.2f\n' % (sigma_p, sigma_p))
                f.write('%11d/%11.2f\n' % (0, self.sigma_D))
                f.write('%11.2f/%11.2f\n' %
                        (self.sigma_D / (self.M - 1), -self.sigma_D /
                         (self.M - 1)))
                f.write('%11.2f/%11d\n' % ((-self.Rp / self.gamma_M), 0))
                f.close()
        # Plot of the basic SN curve according to GL2010
        if plot:
            x = [0, N_1, self.N_D, N_e]
            y = [sigma_1, sigma_1, self.sigma_D, sigma_e]
            plt.loglog(x, y, lw=2, marker='*')
            plt.xlabel('Cycle Numbers')
            plt.ylabel('Stress Amplitude/MPa')
            plt.xlim(10, 10**9)
            plt.yticks([10, 100, 1000])
            plt.annotate(s='(%.2e,%.2f)' % (N_1, sigma_1), xy=(N_1, sigma_1))
            plt.annotate(s='(%.2e,%.2f)' % (self.N_D, self.sigma_D),
                         xy=(self.N_D, self.sigma_D))
            plt.annotate(s='m1=%.2f' % self.m1, xy=(10**3, 142))
            plt.annotate(s='m2=%.2f' % self.m2, xy=(10**7, 40))
            plt.show()
        else:
            pass
        return

from collections import deque, defaultdict
import functools


def _get_round_function(ndigits=None):
    if ndigits is None:
        def func(x):
            return x
    else:
        def func(x):
            return round(x, ndigits)
    return func


def reversals(series, left=False, right=False):
    """Iterate reversal points in the series.
    A reversal point is a point in the series at which the first derivative
    changes sign. Reversal is undefined at the first (last) point because the
    derivative before (after) this point is undefined. The first and the last
    points may be treated as reversals by setting the optional parameters
    `left` and `right` to True.
    Parameters
    ----------
    series : iterable sequence of numbers
    left: bool, optional
        If True, yield the first point in the series (treat it as a reversal).
    right: bool, optional
        If True, yield the last point in the series (treat it as a reversal).
    Yields
    ------
    float
        Reversal points.
    """
    series = iter(series)

    x_last, x = next(series), next(series)
    d_last = (x - x_last)

    if left:
        yield x_last
    for x_next in series:
        if x_next == x:
            continue
        d_next = x_next - x
        if d_last * d_next < 0:
            yield x
        x_last, x = x, x_next
        d_last = d_next
    if right:
        yield x_next


def _sort_lows_and_highs(func):
    "Decorator for extract_cycles"
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for low, high, mult in func(*args, **kwargs):
            if low < high:
                yield low, high, mult
            else:
                yield high, low, mult
    return wrapper


@_sort_lows_and_highs
def extract_cycles(series, left=False, right=False):
    """Iterate cycles in the series.
    Parameters
    ----------
    series : iterable sequence of numbers
    left: bool, optional
        If True, treat the first point in the series as a reversal.
    right: bool, optional
        If True, treat the last point in the series as a reversal.
    Yields
    ------
    cycle : tuple
        Each tuple contains three floats (low, high, mult), where low and high
        define cycle amplitude and mult equals to 1.0 for full cycles and 0.5
        for half cycles.
    """
    points = deque()

    for x in reversals(series, left=left, right=right):
        points.append(x)
        while len(points) >= 3:
            # Form ranges X and Y from the three most recent points
            X = abs(points[-2] - points[-1])
            Y = abs(points[-3] - points[-2])

            if X < Y:
                # Read the next point
                break
            elif len(points) == 3:
                # Y contains the starting point
                # Count Y as one-half cycle and discard the first point
                yield points[0], points[1], 0.5
                points.popleft()
            else:
                # Count Y as one cycle and discard the peak and the valley of Y
                yield points[-3], points[-2], 1.0
                last = points.pop()
                points.pop()
                points.pop()
                points.append(last)
    else:
        # Count the remaining ranges as one-half cycles
        while len(points) > 1:
            yield points[0], points[1], 0.5
            points.popleft()


def count_cycles(series, ndigits=None, left=False, right=False):
    """Count cycles in the series.
    Parameters
    ----------
    series : iterable sequence of numbers
    ndigits : int, optional
        Round cycle magnitudes to the given number of digits before counting.
    left: bool, optional
        If True, treat the first point in the series as a reversal.
    right: bool, optional
        If True, treat the last point in the series as a reversal.
    Returns
    -------
    A sorted list containing pairs of cycle magnitude and count.
    One-half cycles are counted as 0.5, so the returned counts may not be
    whole numbers.
    """
    counts = defaultdict(float)
    round_ = _get_round_function(ndigits)

    for low, high, mult in extract_cycles(series, left=left, right=right):
        delta = round_(abs(high - low))
        counts[delta] += mult
    return sorted(counts.items())

if __name__ == '__main__':
    test = SnGL()
    test.cal(plot=True, save=True)
