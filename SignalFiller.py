#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2019/11/15 13:25
# @Author : Louis Wang
# @Site : https://github.com/pinginf
# @File : SignalFiller.py
# @Software: PyCharm

import numpy as np
from scipy.fftpack import fft
from scipy import signal
import matplotlib.pyplot as plt
from numba import jit

# load read
load = np.loadtxt(r"D:\Code\FatigueDamage\Data\012_002.txt", skiprows=2)[:, -6:-1]
y = load[:, 4]
x = np.linspace(0, 600, len(load))

# Signal filtering
b, a = signal.butter(5, 0.2)
yft = signal.filtfilt(b, a, y)

# FFT
yfft = fft(y).real
yfft_f = fft(yft).real

# plot
fig = plt.figure()
plt.subplot(221)
plt.plot(x, y)
plt.plot(x, yft)
plt.xlim([0, 50])
plt.title('Original/Filtered load')

plt.subplot(222)
plt.plot(x, y - yft)
plt.xlim([0, 50])
plt.title('Relative change')

plt.subplot(223)
plt.plot(x, yfft,'.')
plt.yscale('log')
plt.xlim([-1, 10])
plt.ylim([10**3,10**8])
plt.title('FFT/Original')

plt.subplot(224)
plt.plot(x, yfft_f,'.')
plt.yscale('log')
plt.xlim([-1, 10])
plt.ylim([10**3,10**8])
plt.title('Relative change/Filtered')

fig.tight_layout()
plt.show()


def signal_filte(data,order=8, cutoff=0.3, btype='low', analog=False):
	b, a = signal.butter(order, cutoff,btype,analog)
	return np.array([signal.filtfilt(b, a, data[:,i]) for i in range(data.shape[1])]).T
