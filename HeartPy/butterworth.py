import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter

# For plot_zplane
from matplotlib import patches
from matplotlib.pyplot import axvline, axhline
from collections import defaultdict

# for get_ticks
import math

import os.path as path
import utils as ut
import filter as flt

def plot_zplane(a, b, title='Zeros and Poles', title2=None,
    title_coefficients=True, filename=None):
    '''Plot the complex z-plane given zeros and poles.
    Modified from Christopher Felton, Dec 2011
    '''
    # Get the poles and zeros
    p = np.roots(a)
    z = np.roots(b)
    #print(f"{len(p)=} {len(z)=}")
    
    # get a figure/plot
    plt.figure()
    ax = plt.gca()
    # TODO: should just inherit whatever subplot it's called in?

    # Add unit circle and zero axes    
    unit_circle = patches.Circle((0,0), radius=1, fill=False,
                                 color='black', ls='solid', alpha=0.5)
    ax.add_patch(unit_circle)
    axvline(0, color='0.7')
    axhline(0, color='0.7')
    
    # Plot the poles and set marker properties
    poles = plt.plot(p.real, p.imag, 'x', markersize=9, alpha=0.5)
    
    # Plot the zeros and set marker properties
    zeros = plt.plot(z.real, z.imag,  'o', markersize=9, 
             color='none', alpha=0.5,
             markeredgecolor=poles[0].get_color(), # same color as poles
             )

    # Scale axes to fit
    r = 1.5 * np.amax(np.concatenate((abs(z), abs(p), [1])))
    plt.axis('scaled')
    plt.axis([-r, r, -r, r])
    #ticks = [-1, -.5, .5, 1]
    #plt.xticks(ticks)
    #plt.yticks(ticks)

    #If there are multiple poles or zeros at the same point, put a 
    #superscript next to them.
    #TODO: can this be made to self-update when zoomed?
    
    # Finding duplicates by same pixel coordinates (hacky for now):
    poles_xy = ax.transData.transform(np.vstack(poles[0].get_data()).T)
    zeros_xy = ax.transData.transform(np.vstack(zeros[0].get_data()).T)    

    # dict keys should be ints for matching, but coords should be floats for 
    # keeping location of text accurate while zooming

    # TODO make less hacky, reduce duplication of code
    d = defaultdict(int)
    coords = defaultdict(tuple)
    for xy in poles_xy:
        key = tuple(np.rint(xy).astype('int'))
        d[key] += 1
        coords[key] = xy
    for key, value in d.items():
        if value > 1:
            x, y = ax.transData.inverted().transform(coords[key])
            plt.text(x, y, 
                        r' ${}^{' + str(value) + '}$',
                        fontsize=13,
                        )

    d = defaultdict(int)
    coords = defaultdict(tuple)
    for xy in zeros_xy:
        key = tuple(np.rint(xy).astype('int'))
        d[key] += 1
        coords[key] = xy
    for key, value in d.items():
        if value > 1:
            x, y = ax.transData.inverted().transform(coords[key])
            plt.text(x, y, 
                        r' ${}^{' + str(value) + '}$',
                        fontsize=13,
                        )

    title_used = title
    if title2:
        title_used += ' for ' + title2
    if title_coefficients:
        title_used += f"\n{a=}\n{b=}"
    plt.title(title_used)
    plt.tight_layout()

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
        print('Pole-zero plot saved to ' + str(filename))

def get_ticks(fs):
    ''' Calculate prettier tick labels for semi-log plots
         Uses 10**-1, 10**1, etc if ticks not set
         Assumes values < 1 and > 2k are not needed
    '''
    # Get highest power of 10 of Nyquist frequency = fs/2
    upper = 10 ** math.ceil(math.log10(fs / 2))
    ticks = [1, 10, 20, 50, 100]
    labels = ["1", "10", "20", "50", "100"]
    if upper < 100:
        ticks += [200, 500, 1000]
        labels += ["200", "500", "1k"]
    if upper > 100:
        ticks += [2000, 5000, 10000]
        labels += ["2k", "5k", "10k"]
    if upper > 1000:
        ticks += [20000]
        labels += ["20k"]
    #print(f"{upper=} {labels=}")
    return ticks, labels, upper

def plot_frequency_response(fs, a, b, title='Digital Filter Frequency Response',
        title2=None, title_coefficients=True, use_semilog=True):
    w, h = signal.freqz(b=b, a=a, worN=2000)
    x = w * fs * 1.0 / (2 * np.pi)
    # Avoid divide by zero warning
    np.seterr(divide = 'ignore') 
    y = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.8)
    if use_semilog:
        plt.semilogx(x, y, label='Frequency Response')
        ticks, labels, upper = get_ticks(fs)
        plt.xlim(1, upper)
        plt.xticks(ticks, labels)
    else:
        plt.plot(x, y, label='Group Delay')
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    title_used = title
    if title2:
        title_used += ' for ' + title2
    if title_coefficients:
        title_used += f"\n{fs=}\n{a=}\n{b=}"
    plt.title(title_used)
    plt.grid(which='both', linestyle='-', color='grey')
    plt.show()

def plot_group_delay(fs, a, b, title='Group Delay', title2=None,
        title_coefficients=True, use_semilog=True):
    w, gd = signal.group_delay((b, a))
    x = w * fs * 1.0 / (2 * np.pi)
    y = gd
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.8)
    if use_semilog:
        plt.semilogx(x, y, label='Group Delay')
        ticks, labels, upper = get_ticks(fs)
        plt.xlim(1, upper)
        plt.xticks(ticks, labels)
    else:
        plt.plot(x, y, label='Group Delay')
    plt.xlabel('Frequency, Hz')
    plt.ylabel('Group Delay, samples')
    title_used = title
    if title2:
        title_used += ' for ' + title2
    if title_coefficients:
        title_used += f"\n{fs=}\n{a=}\n{b=}"
    plt.title(title_used)
    plt.grid(which='both', linestyle='-', color='grey')
    plt.legend(loc='best')
    plt.show()

def plot_frequency_response_group_delay(fs, a, b,
    title='Frequency Response and Group Delay', title2=None,
    title_coefficients=True, use_semilog=True):
    # Frequency response
    w, h = signal.freqz(b=b, a=a, worN=2000)
    x1 = w * fs * 1.0 / (2 * np.pi)
    # Avoid divide by zero warning
    np.seterr(divide = 'ignore') 
    y1 = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    # Group delay
    w, gd = signal.group_delay((b, a))
    x2 = w * fs * 1.0 / (2 * np.pi)
    y2 = gd
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.8)
    if use_semilog:
        plt.semilogx(x1, y1, label='Frequency Response')
        plt.semilogx(x2, y2, label='Group Delay')
        ticks, labels, upper = get_ticks(fs)
        plt.xlim(1, upper)
        plt.xticks(ticks, labels)
    else:
        plt.plot(x1, y1, label='Frequency Response')
        plt.plot(x2, y2, label='Group Delay')
    plt.xlabel('Frequency, Hz')
    plt.ylabel('Amplitude, dB / Group Delay, samples')
    title_used = title
    if title2:
        title_used += ' for ' + title2
    if title_coefficients:
        title_used += f"\n{fs=}\n{a=}\n{b=}"
    plt.title(title_used)
    plt.grid(which='both', linestyle='-', color='grey')
    plt.legend(loc='best')
    plt.show()

def plot_frequency_response_2(fs1, a1, b1, fs2, a2, b2, title='Digital Filter Frequency Response'):
    '''Plots the frequency responses for two filters'''
    w, h = signal.freqz(b=b1, a=a1, worN=2000)
    x1 = w * fs1 * 1.0 / (2 * np.pi)
    np.seterr(divide = 'ignore') 
    y1 = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    w, h = signal.freqz(b=b2, a=a2, worN=2000)
    x2 = w * fs2 * 1.0 / (2 * np.pi)
    np.seterr(divide = 'ignore') 
    y2 = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.8)
    plt.semilogx(x1, y1, label=f"{fs1} Hz")
    plt.semilogx(x2, y2, label=f"{fs2} Hz")
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    plt.title(f"{title}\n{a1=}\n{b1=}\n{a2=}\n{b2=}")
    plt.grid(which='both', linestyle='-', color='grey')
    ticks, labels, upper = get_ticks(max(fs1,fs2))
    plt.xlim(1, upper)
    plt.xticks(ticks, labels)
    #plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000], ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def example():
    '''From the web'''
    fs = 44100
    a = 1
    b = [1, -1]
    plot_frequency_response(fs, a, b)

def butter_bandpass(fs, lowcut, highcut, order=5):
    '''Calculates a and b coefficients for a Butterworth filter with these parameters'''
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter_data(data, fs, lowcut, highcut, order=5):
    '''Calculates the filtered values for the input data.
    Uses scipy.signal.lfilter.
    '''
    b, a = butter_bandpass(fs, lowcut, highcut, order=order)
    y = lfilter(b, a, data)
    return y

def filter_data(b, a, data):
    '''Calculates the filtered values for the input data for a filter
    with the given a and b coefficients.
    Uses scipy.signal.lfilter.
    '''
    b, a = butter_bandpass(fs, lowcut, highcut, order=order)
    y = lfilter(b, a, data)
    return y

def plot_butterworth_frequency_response(fs, low_cutoff, high_cutoff, orders, print_coefficients=False):
    '''Plots the frequency response for a Butterworth filter.
    Uses butter_bandpass.
    
        Parameters
        ----------
        fs : float
            The sampling rate.
    
        low_cutoff : float
            The low cutoff.
    
        high_cutoff : float
            The high cutoff.
    
        orders : list of int
            The orders to calculate.
    
        print_coefficients : boolean
            Whether to print the coefficients in Pythion format or not.
            default: False
    
        Returns
        ------
        No return value.
    '''

    # Calculate the frequency response for a few different orders.
    freq_x = []
    freq_y = []
    freq_label = []
    delay_x = []
    delay_y = []
    delay_label = []
    for order in orders:
        b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
        # Generate Python code
        if print_coefficients:
            print(f"\n{order=} {len(a)=} {len(b)=}")
            out = f"A_BUTTERFIELD{order}=["
            for i in range(len(a)):
                out += f"{a[i]}, "
            out += "]"
            print(out)
            out = f"B_BUTTERFIELD{order}=["
            for i in range(len(b)):
                out += f"{b[i]}, "
            out += "]"
            print(out)
        # Frequency response
        w, h = signal.freqz(b, a, worN=2000)
        x = w * fs * 1.0 / (2 * np.pi)
        np.seterr(divide = 'ignore') 
        y = 20 * np.log10(abs(h))
        np.seterr(divide = 'warn')
        freq_x.append(x)
        freq_y.append(y)
        freq_label.append(f"{order=}")
        # Group delay
        w, gd = signal.group_delay((b, a))
        x = w * fs * 1.0 / (2 * np.pi)
        delay_x.append(x)
        delay_y.append(gd)
        delay_label.append(f"{order=}")

    # Plot frequency response
    plt.figure(figsize=(10,6))
    for i in range(len(freq_x)):
        plt.semilogx(freq_x[i], freq_y[i], label=freq_label[i])
    # Mark cutoff
    cutoff_x = [low_cutoff, high_cutoff]
    cutoff_y = [0, 0]
    plt.semilogx(low_cutoff, 0, 'r<', markersize=5)
    plt.semilogx(high_cutoff, 0, 'r>', markersize=5)
    plt.xlabel('Frequency [Hz]')
    plt.title(f"Butterworth Frequency Response ({fs=} Hz {low_cutoff=} {high_cutoff=})")
    plt.grid(which='both', linestyle='-', color='grey')
    plt.legend(loc='best')
    ticks, labels, upper = get_ticks(fs)
    plt.xlim(1, upper)
    plt.xticks(ticks, labels)
    #plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
              #["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.show()

    # Plot group delay
    plt.figure(figsize=(10,6))
    for i in range(len(freq_x)):
        plt.semilogx(delay_x[i], delay_y[i], label=delay_label[i])
    # Mark cutoff
    cutoff_x = [low_cutoff, high_cutoff]
    cutoff_y = [0, 0]
    plt.semilogx(low_cutoff, 0, 'r<', markersize=5)
    plt.semilogx(high_cutoff, 0, 'r>', markersize=5)
    plt.ylabel('Group Delay, samples')
    plt.xlabel('Frequency [Hz]')
    plt.title(f"Butterworth Group Delay ({fs=} Hz {low_cutoff=} {high_cutoff=})")
    plt.grid(which='both', linestyle='-', color='grey')
    plt.legend(loc='best')
    ticks, labels, upper = get_ticks(fs)
    plt.xlim(1, upper)
    plt.xticks(ticks, labels)
    #plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
    #          ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.show()

def get_a0_for_frequency_response(a, b):
    '''Calculates the value of a[0] that puts the maximum of 20log10(abs(h))
    at 0 for a FIR filter with coefficients a, b.
    
        Parameters
        ----------
        a : list of float
            The a coefficients.
    
        b : list of float
            The b coefficients.
    
        Returns
        ------
        a0new : float
            The calculated value of a[0].
    '''
    w, h = signal.freqz(b=b, a=a, worN=2000)
    np.seterr(divide = 'ignore') 
    y = 20 * np.log10(abs(h))
    max_y = np.max(y)
    np.seterr(divide = 'warn')
    a0new = a[0] * 10 ** (max_y / 20)
    return a0new

def test_butterworth(fs, low_cutoff, high_cutoff, order, filename):
    # Filter ecg file
    ecg, _, headers = ut.read_ecg_file(filename)
    x = range(0, len(ecg))
    b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
    filtered = butter_bandpass_filter_data(ecg, fs, low_cutoff, high_cutoff, order=order)
    plt.figure(figsize=(10,6))
    plt.plot(x, ecg, linewidth=1, label='ECG Signal')
    plt.plot(x, filtered, color='#89CFF0', linewidth=1, label=f"Filtered Signal")
    #plt.xlabel('time (seconds)')
    #plt.hlines([-a, a], 0, T, linestyles='--')
    plt.title(f"Butterworth ({fs=} {low_cutoff=} {high_cutoff=} {order=} {len(a)=} {len(b)=}\n{filename}")
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='lower center')
    plt.show()

def test_butterworth_ke(fs, low_cutoff, high_cutoff, order, filename):
    # Filter ecg file
    ecg, _, headers = ut.read_ecg_file(filename)
    x= range(0, len(ecg))
    b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
    filtered = flt.filter_data(a, b, ecg)
    #print(f"{len(ecg)} {ecg[-5:]=}")
    #print(f"{len(filtered)} {filtered[-5:]=}")
    plt.figure(figsize=(10,6))
    plt.plot(x, ecg, linewidth=1, label='ECG Signal')
    plt.plot(x, filtered, color='#89CFF0', linewidth=1, label=f"Filtered Signal")
    #plt.xlabel('time (seconds)')
    #plt.hlines([-a, a], 0, T, linestyles='--')
    plt.title(f"Butterworth Using Filter ({fs=} {low_cutoff=} {high_cutoff=} {order=} {len(a)=} {len(b)=}\n{filename}")
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='lower center')
    plt.show()

def butterworth0():
    '''Butterworth example'''
    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 5000.0
    lowcut = 500.0
    highcut = 1250.0

    # Plot the frequency response for a few different orders.
    plt.figure(1)
    plt.clf()
    for order in [3, 6, 9]:
        b, a = butter_bandpass(fs, lowcut, highcut, order=order)
        w, h = signal.freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')

    # Filter a noisy signal.
    T = 0.05
    nsamples = int(T * fs)
    t = np.linspace(0, T, nsamples, endpoint=False)
    a = 0.02
    f0 = 600.0
    x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
    x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
    x += a * np.cos(2 * np.pi * f0 * t + .11)
    x += 0.03 * np.cos(2 * np.pi * 2000 * t)
    plt.figure(2)
    plt.clf()
    plt.plot(t, x, label='Noisy signal')

    y = butter_bandpass_filter_data(fs, x, lowcut, highcut, order=6)
    plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
    plt.xlabel('time (seconds)')
    plt.hlines([-a, a], 0, T, linestyles='--')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='upper left')

    plt.show()

def ptlow250_130():
    '''Plot Pan-Tompkins Low-Pass at 130 and 250 Hz'''
    fs1 = 250
    a1 = [1, -2, 1]
    b1 = [1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]
    fs2 = 130
    a2 = a1
    b2 = b1
    plot_frequency_response_2(fs1, a1, b1, fs2, a2, b2)

def pan_tompkins_derivative():
    '''Plot Pan-Tompkins Derivative at 130 and 200 Hz'''
    fs1 = 200
    # Empirically choose a[0] to make the gain 0 at the peak
    a1 = [5.4716104436]
    b1 = [2, 1, 0, -1, -2]
    fs2 = 130
    a2 = a1
    b2 = b1
    if True:
        # Calculate what a[0] should be to make y = 0 at the top
        a0new = get_a0_for_frequency_response(a1, b1)
        print(f"{a0new=}")
    plot_frequency_response_2(fs1, a1, b1, fs2, a2, b2,
        title='Pan Tompkins Derivative Frequency Response')

def main():
    print(path.basename(path.normpath(__file__)))
    # Working on computer HR=55
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
    # Walking Falison HR=114
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
    # Walking Blueberry Lake HR=120
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
    # Feb 4 Example New Low Heartrate HR=63
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'

    #ptlow250_130()
    # Pan Tompkins Compare frequency responses for 
    #test_butterworth(130, 5, 15, order=3, filename=filename)
    #test_butterworth(130, 5, 20, order=3, filename=filename)
    #test_butterworth_ke(130, 5, 20, order=3, filename=filename)

    # Pan Tompkins Compare frequency responses for Pan Thompkins Derivative
    if False:
        pan_tompkins_derivative()

    # Get Butterworth frequency responses for a set of orders
    # Optionally print the coefficients in Python code format
    if True:
        plot_butterworth_frequency_response(130, 5, 20, orders=[2,3,4,5,6,9],
                                           print_coefficients=False)

    # Plots for Pan Tompkins derivative
    if False:
        a = flt.A_DERIVATIVE
        b = flt.B_DERIVATIVE
        title2 = 'Pan Tompkins Derivative'
        plot_frequency_response(130, a, b, title='Digital Filter Frequency Response',
                               title2=title2, use_semilog=True)
        #plot_group_delay(130, a, b, title2=title2, use_semilog=False)
        #plot_frequency_response_group_delay(130, a, b, title2=title2, use_semilog=True)
        plot_zplane(a, b, title2=title2)

    # Plots for Butterworth
    if True:
        a = flt.A_BUTTERWORTH3
        b = flt.B_BUTTERWORTH3
        title2= 'Butterworth3 (fs=130 low_cutoff=5 high_cutoff=20)'
        plot_frequency_response(130, a, b, title='Digital Filter Frequency Response',
                               title2=title2, use_semilog=True, title_coefficients=False)
        #plot_group_delay(130, a, b, title2=title2, use_semilog=False, title_coefficients=False)
        #plot_frequency_response_group_delay(130, a, b, title2=title2, use_semilog=True, title_coefficients=False)
        plot_zplane(a, b, title2=title2, title_coefficients=False)

if __name__ == "__main__":
    main()
