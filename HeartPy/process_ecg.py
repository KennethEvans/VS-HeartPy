'''
This module simulates real-time processing of an ECG signal.

run_real_time() processes the input file using moving windows.
run_process_after() processes the input file all at once (not real-time).
'''

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.fft import fft, fftfreq
import os
import time
from collections.abc import Iterable
import utils as ut
import filter as flt

# Sampling rate. These algorithme are based on this particular sampling rate.
FS = 130.0
# Width of plot to show
PLOT_WIDTH_SEC = 3
PLOT_WIDTH_SAMPLES = 3 * FS
PLOT_WIDTH_MARGIN = .1
# Number of mV to show
PLOT_HEIGHT_MV = 2
# Data window size. Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# Moving average window size.
MOV_AVG_WINDOW = 20
# Set whether to be fast (blit=True) with no x ticks or blit=False and slow
FAST = True

def plot_fft(data, fs, title = 'FFT', filename=None):
    if filename:
      title = f"{title}\n{filename}"
    n = len(data)
    t = 1 / fs
    yf = fft(data)
    xf = fftfreq(n, t)[:n // 2]
    plt.plot(xf, 2.0 / n * np.abs(yf[0:n // 2]))
    plt.grid()
    plt.title(title)
    plt.xlabel('frequency')
    plt.ylabel('fft')
    plt.show()

def plot_all(ecg, bandpass, deriv, square, avg, score, title='ECG'):
    plt.figure(figsize=(14,7))
    plt.subplots_adjust(top=0.83)
    plt.suptitle(title)

    plt.subplot(2, 3, 1)
    plt.title('ECG')
    plt.plot(ecg)

    plt.subplot(2, 3, 2)
    plt.title('Bandpass')
    plt.plot(bandpass)

    plt.subplot(2, 3, 3)
    plt.title('Derivative')
    plt.plot(deriv)

    plt.subplot(2, 3, 4)
    plt.title('Square')
    plt.plot(square)

    plt.subplot(2, 3, 5)
    plt.title('Moving Average')
    plt.plot(avg)

    plt.subplot(2, 3, 6)
    plt.title('Score')
    plt.plot(score)

    #plt.subplot(2, 3, 6)
    #plt.title('f5')
    #plt.plot(f5)

    plt.show()

def plot_2_values(ecg, timevals, vals1, vals2=None, label1='1', label2='2',
        title='ECG Filtering Results', use_time_vals=True,
        xlim=None):
    plt.figure(figsize=(12,7))
    if  not timevals: use_time_vals = False
    if use_time_vals == False:
        if ecg:
            plt.plot(ecg, label='ECG')
        if vals1:
            plt.plot(vals1, label = label1)
        if vals2:
            plt.plot(vals2, label = label2)
    else:
        plt.xlabel('time, sec')
        if ecg:
            plt.plot(timevals, ecg, label='ECG')
        if vals1:
            plt.plot(timevals, vals1, label = label1)
        if vals2:
            plt.plot(timevals, vals2, label = label2)
    if(xlim):
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.title(title)
    #plt.xlabel('time')
    #plt.ylabel('mV')
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def run_process_after():
    ecg, headers = ut.read_ecg_file(filename)
    necg = len(ecg)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    if description:
        title = f'Process After\n{filename}\n{description}'
    else:
        title = f'Process After\n{filename}'

    # Set up the plot
    bandpass = []
    deriv = []
    square= []
    avg=[]
    score=[]
    cur_ecg = []
    cur_filt = []
    filter = []
    cur_x = []

    # Butterworth
    a = flt.A_BUTTERWORTH3
    b = flt.B_BUTTERWORTH3
    data = ecg
    bandpass = flt.filter_data(a, b, data)
    print(f"{len(bandpass)=}")

    # Derivative
    a = flt.A_DERIVATIVE
    b = flt.B_DERIVATIVE
    data = bandpass
    deriv = flt.filter_data(a, b, data)
    print(f"{len(deriv)=}")
    print(f"{deriv[:5]=}")

    # Square
    data = deriv
    for i in range(necg):
        new = data[i] * data[i]
        square.append(new)
    print(f"{len(square)=}")

    # Moving average
    a = [1]
    b = [1 / MOV_AVG_WINDOW] * MOV_AVG_WINDOW
    ##DEBUG
    #b = [1 / 5] * 5
    data = square
    avg = flt.filter_data(a, b, data)
    print(f"{len(avg)=}")

    # Score
    threshold = .015
    data = avg
    for i in range(necg):
       new = flt.score(data[i], threshold)
       score.append(new)

    ## Remove low frequency
    #avg = avg - np.mean(avg)

    # FFT
    if False:
        plot_fft(ecg, FS, filename = filename)
        plot_fft(avg, FS, title = 'FFT (Average)', filename = filename)

    #vals1 = bandpass
    #label1 = 'Bandpass'
    #vals1 = deriv
    #label1 = 'Derivative'
    #vals1 = square
    #label1 = 'Square'
    #vals1 = square
    #label1 = 'Score'

    if True:
        vals1 = np.ndarray.tolist(np.array(avg) * 10.)
        label1 = 'Average x 10'
        vals2 = score
        label2 = 'Score'
    else:
        # DEBUG
        vals1 = square
        label1 = 'Square'
        vals2 = np.ndarray.tolist(np.array(avg) * 10.)
        label2 = 'Average x 10'

    #plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
    #    use_time_vals=True, xlim=[0, PLOT_WIDTH_SEC])
    plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
        use_time_vals=True)

    plot_all(ecg, bandpass, deriv, square, avg, score, title=title)

def run_real_time():
    '''
    This version collects cur_bandpass, cur_deriv, etc. separately
    rather than in one cur_filter.
    '''
    ecg, headers = ut.read_ecg_file(filename)
    necg = len(ecg)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    if description:
        title = f'Real Time Processing\n{filename}\n{description}'
    else:
        title = f'Real Time Processing\n{filename}'

    # Set up the plot
    bandpass = []
    deriv = []
    square= []
    avg=[]
    score=[]
    cur_butterworth = []
    cur_deriv = []
    cur_square= []
    cur_avg=[]
    cur_score=[]
    cur_ecg = []
    cur_x = []
    
    keep = DATA_WINDOW # Keep this many items

    # Loop over ECG values
    for i in range(0, necg):
       if len(cur_ecg) == keep:
            cur_ecg.pop(0)
       cur_ecg.append(ecg[i])
       #print(f"{i} {len(cur_ecg)=}")

       # Butterworth
       input = cur_ecg
       if len(cur_butterworth) == keep:
            cur_butterworth.pop(0)
       cur_butterworth.append(0) # Doesn't matter
       new = flt.butterworth3(input, cur_butterworth)
       cur_butterworth[-1] = new
       bandpass.append(new)

       # Derivative
       input = cur_butterworth
       if len(cur_deriv) == keep:
            cur_deriv.pop(0)
       cur_deriv.append(0) # Doesn't matter
       new = flt.derivative(input, cur_deriv)
       cur_deriv[-1] = new
       deriv.append(new)

       # Square
       input = cur_deriv
       if len(cur_square) == keep:
            cur_square.pop(0)
       cur_square.append(0) # Doesn't matter
       new = flt.square(input, cur_square)
       cur_square[-1] = new
       square.append(new)

       # Moving average
       input = cur_square
       if len(cur_avg) == keep:
            cur_avg.pop(0)
       cur_avg.append(0) # Doesn't matter
       new = flt.moving_average1(input, MOV_AVG_WINDOW)
       cur_avg[-1] = new
       avg.append(new)
       #print(f"{i} avg = {avg[0:5]}")

       # Score
       input = cur_avg
       threshold = .01
       new = flt.score(input[-1], threshold)
       score.append(new)

    ## Remove low frequency
    #avg = avg - np.mean(avg)

    # FFT
    if False:
        plot_fft(ecg, FS, filename = filename)
        plot_fft(avg, FS, title = 'FFT (Average)', filename = filename)

    #vals1 = bandpass
    #label1 = 'Bandpass'
    #vals1 = deriv
    #label1 = 'Derivative'
    #vals1 = square
    #label1 = 'Square'
    #vals1 = square
    #label1 = 'Score'
    vals1 = np.ndarray.tolist(np.array(avg) * 10.)
    label1 = 'Average x 10'
    vals2 = score
    label2 = 'Score'
    #plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
    #    use_time_vals=True, xlim=[0, PLOT_WIDTH_SEC])
    plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
        use_time_vals=True)

    plot_all(ecg, bandpass, deriv, square, avg, score, title=title)

def main():
    global filename

    print(os.path.basename(os.path.normpath(__file__)))

    # Working on computer HR=55
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
    # Walking Falison HR=114
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
    # Walking Blueberry Lake HR=120
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
    # Feb 4 Example New Low Heartrate HR=63
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
    # Feb 6 Walking
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-06_13-52.csv'

    #test()
    #test2()

    if True:
        # Processing as we go
        run_process_after()
    if True:
        # Processing as we go
        run_real_time()

if __name__ == "__main__":
    main()
