# Module for plotting functions
from usp_functions import *
import numpy as np
from scipy.optimize import leastsq
import pylab as plt
import numpy, scipy.optimize
from scipy import interpolate
from numpy import arange
import math


def fft_hann(data, fs_input):
    # Runs a Fast Fourier Transform with a Hann filter for 'smoothing' edges of data window for reducing high freq

    fft_freq, fft_values = signal.welch(data, fs=fs_input, window='hamming', nperseg=len(data), noverlap=None,
                                        nfft=len(data)*5, detrend='constant', return_onesided=True, scaling='density',
                                        axis=-1)

    return fft_freq, fft_values


def get_stimulus_phase(usp_object, window_start, window_end, peak=True, plot_events=False):

    phase_angles = []
    phase_degrees = []

    for i in range(len(usp_object.delta_f[0])):
        time = np.linspace(window_start,window_end,(window_end-window_start)*usp_object.fs)
        y_vals = usp_object.delta_f[0][i][window_start*usp_object.fs:window_end*usp_object.fs]
        res = fit_sine(time, y_vals)
        phase_degrees.append(get_phase_degrees(time, res["fit_func"](time), res["period"], usp_object, plot=False))
        phase_angles.append(get_slope_angle(time, res["fit_func"](time), usp_object.t_pre, plot=False))

        if plot_events:
            print("Amplitude=%(amp)s, Angular freq.=%(omega)s, phase=%(phase)s, offset=%(offset)s, Max. Cov.=%(max_cov)s" % res)
            plt.plot(time, y_vals, "-k", label="y", linewidth=2)
            plt.plot(time, res["fit_func"](time), "r-", label="y fit curve", linewidth=2)
            plt.legend(loc="best")
            plt.show()

    return phase_angles, phase_degrees


def get_phase_degrees(time, vals, period, usp_object, plot=False):

    max_index = np.argmax(vals)
    stim_index = usp_object.t_pre * usp_object.fs
    phase = (max_index - stim_index) / (period*usp_object.fs)  # how many periods away from stim is peak
    if phase < 0:
        phase = phase % -1
        phase = phase * -360

    else:
        phase = phase % 1
        phase = 360 - (phase * 360)

    return phase


def get_slope_angle(x, y, a, plot=False):
    # interpolate the data with a spline

    spl = interpolate.splrep(x, y)
    tangent_time = arange(a-2, a+5)
    f_a = interpolate.splev(a, spl, der=0)     # f(a)
    f_prime = interpolate.splev(a, spl, der=1) # f'(a)
    tan = f_a + f_prime * (tangent_time-a)
    slope = (tan[1] - tan[0] / tangent_time[1]-tangent_time[0])
    slope_angle = math.degrees(math.atan(slope))
    if plot:
        plt.plot(a, f_a, 'om', tangent_time, tan, '--r')

    return slope_angle


def fit_sine(time, y_vals):
    'Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'

    ff = numpy.fft.fftfreq(len(time), (time[1] - time[0]))  # assume uniform spacing
    Fyy = abs(numpy.fft.rfft(y_vals))
    guess_freq = abs(
         ff[numpy.argmax(Fyy[1:]) + 1])  # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(y_vals) * 2. ** 0.5
    guess_offset = numpy.mean(y_vals)
    guess = numpy.array([guess_amp, 2. * numpy.pi * guess_freq, 0., guess_offset])
    popt, pcov = scipy.optimize.curve_fit(sin_func, time, y_vals, p0=guess, maxfev=10000)
    a, w, p, c = popt
    f = w / (2. * numpy.pi)
    fit_func = lambda t: a * numpy.sin(w * t + p) + c

    return {"amp": a, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1. / f, "fit_func": fit_func,
            "max_cov": numpy.max(pcov), "raw_res": (guess, popt, pcov)}


def sin_func(t, a, w, p, c):
    return a * numpy.sin(w * t + p) + c
