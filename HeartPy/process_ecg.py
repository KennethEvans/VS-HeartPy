'''
This module simulates real-time processing of an ECG signal.

run_real_time() processes the input file using moving windows.
run_process_after() processes the input file all at once (not real-time).
'''

from pickle import TRUE
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.fft import fft, fftfreq
import sys
import os
import time
from datetime import datetime
from collections.abc import Iterable
import utils as ut
import filter as flt
from filter import Moving_Average

# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0
# Minimum length of score interval
# Determined from mimumum QRS is .08 to .10 sec
MAX_QRS_LENGTH = round(.12 * FS)
# Data window size.  Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# Moving average window size.
MOV_AVG_WINDOW = 20
# Moving average height window size
MOV_AVG_HEIGHT_WINDOW = 5
# Moving average height default
MOV_AVG_HEIGHT_DEFAULT = .025
# Moving average height threshold factor
#     Note: threshold = MOV_AVG_HEIGHT_THRESHOLD_FACTOR * Moving_average.avg()
MOV_AVG_HEIGHT_THRESHOLD_FACTOR = .4

# Width of plot to show in zoomed view
ZOOMED_PLOT_WIDTH_SEC = 6

def plot_peaks(ecg, ecg_x, peak_indices, title='Detected Peaks', filename=None,
               xlim = None):
    peak_vals = [ecg[i] for i in peak_indices]
    peak_x = [ecg_x[i] for i in peak_indices]
    plt.figure(figsize=(10,6))
    #plt.subplots_adjust(top=0.8)
    plt.plot(ecg_x, ecg)
    plt.plot(peak_x, peak_vals, "ro")
    title_used = f"{title} ({len(peak_indices)} Peaks)"
    if filename:
        title_used += f"\n{filename}"
    plt.title(title_used)
    if xlim:
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.xlabel('time, sec')
    plt.tight_layout()
    plt.show()

def calculate_rr(x_ecg, peak_indices, window=5):
    rr = []
    rr_inv = []
    x_rr = []
    for i in range(len(peak_indices)):
        if i == 0:
            rr.append(float("NaN"))
            rr_inv.append(0.)
            x_rr.append(x_ecg[peak_indices[i]])
            continue
        delta = x_ecg[peak_indices[i]] - x_ecg[peak_indices[i - 1]]
        rr.append(delta * 1000.)
        rr_inv.append(60. / delta)
        x_rr.append(x_ecg[peak_indices[i]])
        #print(f'{rr[i]=} {x_ecg[i]=} {x_rr[i]=}')
    hr = flt.moving_average_uniform_filter(rr_inv, window)
    # Mask out the first window values
    for i in range(window):
        hr[i] = float("NaN")
    #print(f'{x_rr=}')
    #print(f'{rr=}')
    return rr, hr, x_rr

def plot_rr(x, rr, hr, title='HR and RR Values', window = None, filename=None):
    lenX = len(x)
    y_fixed_val = 30
    y0 = [y_fixed_val] * lenX
    fig = plt.figure(figsize=(10,6))
    plt.plot(x, rr, "-ob", label='RR', markersize=4)
    plt.ylim(bottom=0)
    plt.ylabel("RR, ms");
    title_used = f"{title}"
    if filename:
        title_used += f"\n{filename}"
    plt.title(title_used)
    plt.xlabel('time, sec')
    # HR on 2nd y axis
    ax2 = plt.gca().twinx()
    if window:
        hr_label = f'HR (avg over {window} peaks)'
    else:
        hr_label = 'HR'
    ax2.plot(x, hr, '-or', label=hr_label, markersize=4)
    ax2.plot(x, y0, "-sg", label='Peaks', markersize=4)
    ax2.set_ylabel('HR, bpm');
    ax2.set_ylim(bottom=0)
    ax = plt.gca()
    fig.legend(loc='lower left', framealpha=0.6, bbox_to_anchor=(0,0),
              bbox_transform=ax.transAxes)
    plt.tight_layout()
    plt.show()

def plot_fft(data, fs, title='FFT', filename=None):
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
        title='ECG Filtering Results', use_time_vals=True, xlim=None):
    '''Plots 2 or optionally 3 sets of data using the same x values.
    
        Parameters
        ----------
        ecg : list of float
            The first set of data.
    
        timevals : list of float
            The x values corresponding to time.
    
        vals1 : list of float
            The second set of data. May be None
    
        vals2 : list of float
            The third set of data. May be None
            default: None
    
        label1 : str
            Label for second set of data. (first is 'ECG'.)
            default: '1'
    
        label2 : str
            Label for third set of data.
            default: '2'
    
        title : str
            Title for the plot.
            default: 'ECG Filtering Results'

        use_time_vals : boolean
            Whether to use sample number or the given time values.
            default: True

        xlim : list of two float values for the start and end of the x axis.
            Sets the x axis limits.
            default: None
     '''
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
    if xlim:
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.title(title)
    #plt.xlabel('time')
    #plt.ylabel('mV')
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def shift_score(score, shift):
    ''' Shifts a list shift places to the right or left depending on the sign
    of shift.
    '''
    shifted = score.copy()
    if shift >= 0:
        for i in range(shift):
            shifted.pop()
            shifted.insert(0, shifted[0])
    else:
        for i in range(-shift):
            shifted.pop(0)
            shifted.append(shifted[-1])
    return shifted

def run_process_after(filename):
    ecg, _, headers = ut.read_ecg_file(filename)
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
    square = []
    avg = []
    score = []
    cur_ecg = []
    cur_filt = []
    filter = []
    cur_x = [i / FS for i in range(necg)]

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
    # TODO Fix this for recalculating moving average height
    # and starting threshold = .01
    #     (MOV_AVG_HEIGHT_THRESHOLD_FACTOR * MOV_AVG_HEIGHT_DEFAULT)
    threshold = .015
    data = avg
    for i in range(necg):
       new = flt.score(data[i], threshold)
       score.append(new)

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

    # Plot the first PLOT_WIDTH_SEC
    if False:
        plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2,
            title=title, use_time_vals=True, xlim=[0, ZOOMED_PLOT_WIDTH_SEC])

    # Plot with score shifted
    if True:
        shift = -17
        shifted = shift_score
        plot_2_values(ecg, cur_x, vals1, shifted, label1=label1,
            label2=f"score shifted by {shift}", title=title, use_time_vals=True)

    # Normal plot
    plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
        use_time_vals=True)

    plot_all(ecg, bandpass, deriv, square, avg, score, title=title)

def score_real_time(filename, show_progress = False):
    '''
    This is the scoring part of run_real_time.
    Returns ecg, x_ecg, peak_indices, headers, avg, score, ecg_ext, cur_x
    '''
    ecg, _, headers = ut.read_ecg_file(filename)
    necg = len(ecg)

    # Set up the plot
    bandpass = []
    deriv = []
    square = []
    avg = []
    score = []
    cur_butterworth = []
    cur_deriv = []
    cur_square = []
    cur_avg = []
    cur_score = []
    cur_ecg = []

    #Scoring
    score_offset = 18 # This is the group delay, used for searching ecg for maxima
    peaks = []
    peak_indices = []
    max_avg = []
    max_avg_indices = []
    scoring = False
    score_start = 0
    score_end = 0
    max_index = min_index = -1

    print_scores = False

   # Create an extended array with zeros in the last score_offset places
    ecg_ext = ecg.copy() + [0.0] * score_offset
    necgext = len(ecg_ext)
    x_ecg = [i / FS for i in range(necg)]
    cur_x = [i / FS for i in range(necgext)]
   
    keep = DATA_WINDOW # Keep this many items

    # Keep a moving average of the moving average heights
    # Initialize it assuming avg peaks are MOV_AVG_HEIGHT_DEFAULT high
    moving_average_height_list = [MOV_AVG_HEIGHT_DEFAULT] * MOV_AVG_HEIGHT_WINDOW
    moving_average_height = Moving_Average(MOV_AVG_HEIGHT_WINDOW,
        moving_average_height_list)
    threshold = MOV_AVG_HEIGHT_THRESHOLD_FACTOR * moving_average_height.avg()
    #print(f'Starting {threshold=} {moving_average_height.avg()=}')

    # Loop over ECG values
    for i in range(0, necgext):
        if len(cur_ecg) == keep:
            cur_ecg.pop(0)
        cur_ecg.append(ecg_ext[i])
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

        # Score
        input = cur_avg
        # Base the threshold on the current average of the moving average heights
        # input[-1] = cur_avg[-1] = avg[i] is the last value of cur_avg
        new = flt.score(input[-1], threshold)
        score.append(new)
        # Process finding the peaks
        if i == 0:
            # At the first point of the [extended] data
            score_start = score_stop = -1
            if new == 1:
                score_start = i
                scoring = True
        elif i == necgext - 1:
            # At last point in the [extended] data so stop
            if scoring:
                score_end = i
                scoring = False
        else:
            if scoring and new == 0:
                score_stop = i
                scoring = False
            if not scoring and new == 1:
                score_start = i
                scoring = True
        if not scoring and score_stop == i:
            # End of interval, process the score
            score_len = score_stop - score_start
            delta_rs = 1000000 # Something large
            if min_index > -1 and max_index > -1:
                delta_rs = min_index - max_index
            # Criterion for using this interval as containing a valid QRS complex
            use_interval = delta_rs >= 0 and delta_rs <= MAX_QRS_LENGTH and\
               score_len > MAX_QRS_LENGTH / 2
            if use_interval:
                peaks.append(ecg_ext[max_index])
                peak_indices.append(max_index)
                # Recalculate the threshold
                moving_average_height.add(max_avg_height)
                threshold = MOV_AVG_HEIGHT_THRESHOLD_FACTOR * moving_average_height.avg()
                #print(f'{i} New {threshold=} {moving_average_height.avg()=}')
            if show_progress:
                print(f'use_interval, min_index, max_index, delta_rs, '
                    f'min_ecg, max_ecg, score_len: '
                    f'{use_interval!s:^5} {min_index:4d} {max_index:4d} '
                    f'{delta_rs:4d} {min_ecg:6.3f} {max_ecg:6.3f} {score_len:4d}')

            max_index = min_index = -1
        if scoring:
            if score_start == i:
                # Start of interval, set up scoring
                if i >= score_offset:
                    max_index = min_index = i - score_offset
                    max_ecg = min_ecg = ecg_ext[i - score_offset]
                else:
                    max_index = min_index = -1
                    max_ecg = -sys.float_info.max
                    min_ecg = sys.float_info.max
                max_avg_height = input[-1]
            else:
                # In interval, accumulate data
                if i >= score_offset:
                    last = ecg_ext[i - score_offset]
                    if last > max_ecg:
                        max_ecg = last
                        max_index = i - score_offset
                    if last < min_ecg:
                        min_ecg = last
                        min_index = i - score_offset
                if input[-1] > max_avg_height:
                    max_avg_height = input[-1]

    return ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score, ecg_ext, cur_x

def run_real_time(filename, show_progress = False, write_csv = False,
            plot_peaks_zoomed = False,
            plot_analysis = True, plot_all_filter_steps = False,
            plot_rr = False):
    '''
    This version collects cur_bandpass, cur_deriv, etc. separately
    rather than in one cur_filter. It handles mostly plotting and printing.
    The algorithm is run in score_real_time.
    '''
    ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score,\
       ecg_ext, cur_x = score_real_time(filename, show_progress = show_progress)

    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    if description:
        title = f'Real Time Processing\n{filename}\n{description}'
    else:
        title = f'Real Time Processing\n{filename}'

    # Plot peaks
    if True:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename
        plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename)

    # Plot peaks zoomed
    if plot_peaks_zoomed:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename
        plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename,
                   xlim=[0, 10])

    # Calculate and plot RR and HR
    if plot_rr:
        window = 10
        rr, hr, x_rr = calculate_rr(x_ecg, peak_indices, window=window)
        plot_rr(x_rr, rr, hr, window=window, filename=peak_filename)

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

    # Debugging
    #
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

    # Plot with specified x axis
    if plot_analysis:
        shift = -18
        shifted = shift_score(score, shift)
        label2 = f"Score shifted by {shift}"
        plot_2_values(ecg_ext, cur_x, vals1, shifted, label1=label1, label2=label2,
            title=title, use_time_vals=True, xlim=[0, ZOOMED_PLOT_WIDTH_SEC])

    # Plot with score shifted
    if plot_analysis:
        shift = -18
        shifted = shift_score(score, shift)
        label2 = f"Score shifted by {shift}"
        plot_2_values(ecg_ext, cur_x, vals1, shifted, label1=label1,
            label2=label2, title=title, use_time_vals=True)

    # Normal plot with score not shifted
    if False and plot_analysis:
        label2 = f"Score not shifted"
        plot_2_values(ecg_ext, cur_x, vals1, vals2, label1=label1, label2=label2,
            title=title, use_time_vals=True)

    # Plot all filter steps
    if plot_all_filter_steps:
        plot_all(ecg_ext, bandpass, deriv, square, avg, score, title=title)

    # Write the result
    if write_csv:
        outfile = r'data\ecg_test_data.csv'
        now = datetime.now()
        timestamp = now.strftime('%Y-%m-%d %H:%M:%S')
        print(f'\nWriting {outfile} at {timestamp}')
        ut.write_ecg_file(ecg_ext, peak_indices, headers, outfile,
            simulation_str = f'HeartPy.process_ecg.run_real_time {timestamp}')

def check_real_time(filename, show_progress = False):
    ''' Runs score_real_time and generates plots to check for how well the
    algorithm did.
    '''
    ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score,\
       ecg_ext, cur_x = score_real_time(filename, show_progress = show_progress)
    print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')

    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    if description:
        title = f'Real Time Processing\n{filename}\n{description}'
    else:
        title = f'Real Time Processing\n{filename}'

    # Print results
    npeaks = ut.find_header_item(headers, 'npeaks')
    if npeaks:
        print(f'\nFound {len(peak_indices)} peaks. File had {npeaks} peaks.')
    else:
        print(f'\nFound {len(peak_indices)} peaks.')

    # Plot peaks zoomed
    if True:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename

        for i in range(3):
            plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename,
                xlim=[i * 10, i * 10 + 10])
        

def main():
    global filename
    print(os.path.basename(os.path.normpath(__file__)))

     # Set prompt to use default filename or prompt with a FileDialog
    prompt = False
    if prompt:
        file_names = ut.prompt_for_files(title='Choose a Polar CVS file with ECG Data',
            multiple=False, type='csv')
        nFiles = len(file_names)
        if nFiles == 0:
            print('Canceled')
            return;
        filename = file_names[0]
    else:
        # 0 Working on computer HR=55
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
        # 1 Walking Falison HR=114
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
        # 2 Walking Blueberry Lake HR=120
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
        # 3 Feb 4 Example New Low Heartrate HR=63
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
        # 4 Feb 6 Walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-06_13-52.csv'
        # After walking. Device, internal HR differ
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_15-08.csv'
        # Before walking 2. Device HR changed from 94 to 54
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-51.csv'
        # Before walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-44.csv'
        # Sitting down. Little out of breath.
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-03_18-43.csv'
        # 2-20-2022 Before walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-44.csv'
        # 2-15-2022-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-15_12-58.csv'
        # 2-18-16 Has extra processed peak
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-18_16-06.csv'
        # 5-13-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-05-13_12-42.csv'
        # 6-15-18 With new algorithm in KE.Net ECG
        filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-06-15_18-30.csv'
    #test()
    #test2()

    if False:
        # Processing entire file
        run_process_after(filename)
    if True:
        # Check the file vs processed peaks
        print('Running check_real_time')
        check_real_time(filename, show_progress = False)
    if True:
        # Processing as we go
        print('Running run_real_time')
        write_csv = False
        run_real_time(filename, show_progress = True, write_csv = write_csv)

if __name__ == "__main__":
    main()
