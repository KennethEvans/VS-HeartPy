'''
This module simulates real-time processing of an ECG signal.

run_real_time() processes the input file using moving windows.
run_process_after() processes the input file all at once (not real-time).
'''

from operator import truediv
from pickle import FALSE, TRUE
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
from datetime import datetime
import utils as ut
import filter as flt

import score0 as sc0
import score1 as sc1
import score1a as sc1a
import score1b as sc1b
import score1f as sc1f
import score1g as sc1g
import score2 as sc2
import wavelet as wav

# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0
# Data window size.  Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# HR 200 Interval
HR_200_INTERVAL = int(60.0 / 200.0 * FS)

# Width of plot to show in custom zoomed view
ZOOMED_PLOT_WIDTH_SEC = 6

class curve_data:
    ''' Class to hold data that can be plotted '''
    def __init__(self, data, peaks, name, color, marker=None, markersize=6):
        self.data = data
        self.peaks = peaks
        self.name = name
        self.color = color
        self.markersize = markersize
        self.marker = marker

def plot_peaks(ecg, ecg_x, peak_indices,
               title='Detected Peaks', filename=None, xlim = None, use_points=False,
               curves=None, do_intervals=False, xlabel='time, sec'):
    '''Plots ecg and the peaks. Optionally adds other curves.
    ecg: The ECG values, assumed the same for all curves
    ecg_x: The x values for the ECG
    xlim: None or xlim[0] = start, xlim[1]=end values for x axis
    use_points: Use markers to show individual points
        markersize = 2
        default is no markers
    do_intervals: Show the HR_200_INTERVAL's
    curves: List of curve_data to use for extra plots
    '''
    peak_vals = [ecg[i] for i in peak_indices]
    peak_x = [ecg_x[i] for i in peak_indices]

    plt.figure(figsize=(10,6))
    #plt.subplots_adjust(top=0.8)
    # Use points instead of line for checking
    if use_points:
        marker = "o"
    else:
        marker = None
    do_legend = False
    if curves:
        do_legend = True
        for i in range(len(curves)):
            curve = curves[i]
            if 'Score' in curve.name:
                marker_used = None
            else:
                marker_used = marker
            if curve.data:
                plt.plot(ecg_x, curve.data, color=curve.color,
                        marker=marker_used, markersize = 2,
                        label=curve.name)
            if curve.peaks:
                y = [ecg[i] for i in curve.peaks]
                x = [ecg_x[i] for i in curve.peaks]
                plt.plot(x, y, color=curve.color,
                        marker=None, markersize = 2,
                        label=f'{curve.name} Peaks')
    plt.plot(ecg_x, ecg, marker=marker, markersize = 2, label='ECG')
    
    if do_intervals:
        intervals = []
        intervals_x = []
        for i in range(len(ecg)):
            if i % HR_200_INTERVAL == 0:
                intervals.append(.25)
                intervals_x.append(ecg_x[i])
        plt.plot(intervals_x, intervals, color='tomato', marker='|')

    plt.plot(peak_x, peak_vals, "ro")

    # Show axis
    plt.axhline(0, color='black', linewidth=.75)

    if do_legend:
        plt.legend(loc='lower right', framealpha=0.6)
    title_used = title
    title_used = f'{title_used} ({len(peak_indices)} Peaks)'
    if filename:
        title_used += f'\n{filename}'
    plt.title(title_used)
    if xlim:
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.xlabel(xlabel)
    plt.tight_layout()
    plt.show()

def plot_comparison(scores, filename=None,
        height_fract = 0, title='Peak Detection Comparison',
        show_segments = False, show_segments_only = False,
        sample_rate = FS):
    #ecg = electrocardiogram()[2000:4000]
    # Determine markersizes
    n_scores = len(scores)
    # Assume all ecg are the same, use the first
    ecg = scores[0].data
    if len(ecg) == 0:
        print('There are no ecg values')
        return
    markersize = 2
    for i in range(n_scores):
        score = scores[i]
        if score:
            markersize = markersize + 2
            score.markersize = markersize

    if not show_segments_only:
        start = 0
        end = len(ecg)
        time = []
        for i in range(start, end):
            time.append(i / sample_rate)
        plt.figure(figsize=(10,6))
        plt.plot(time, ecg)

    # Do in reverse order
    for i in range(n_scores-1, -1, -1):
        score = scores[i]
        if score:
            npeaks = add_peaks(score, 0, len(score.data), sample_rate,
                marker = 'o', color=score.color, markersize = score.markersize,
                label = f'{score.name} ({len(score.peaks)} peaks)')

    if not filename:
        plt.title(f'{title}')
    else:
        plt.title(f'{title}\n{filename}')
    # Show axis
    plt.axhline(0, color='black', linewidth=.75)
    plt.xlabel(f'time, sec (sample_rate={sample_rate})')
    plt.legend(loc='lower right', framealpha=0.6)
    plt.tight_layout
    plt.show()
    
    if not show_segments:
        return
  
  # Plot magnified
    npoints = len(ecg)
    nplots = 3
    npoints1 = (int)(npoints / nplots)
    for nplot in range(0, nplots):
        start = nplot * npoints1
        if nplot == nplots - 1:
            end = npoints
        else:
            end = (nplot+1) * npoints1
        time = []
        for i in range(start, end):
            time.append(i / sample_rate)
        plt.figure(figsize=(10,6))
        plt.plot(time, ecg[start : end])

        # Do in reverse order
        for i in range(n_scores-1, -1, -1):
            score = scores[i]
            if score:
                npeaks = add_peaks(score, start=start, end=end, sample_rate=FS,
                    marker='o', color=score.color, markersize=score.markersize,
                    label=f'{score.name} ({len(score.peaks)} peaks)')

        plt.title(f'{title}\n{filename}')
        plt.xlabel(f'time, sec (sample_rate={sample_rate})')
        plt.legend(loc='lower right', framealpha=0.6)
        plt.tight_layout()
        plt.show()

def add_peaks(score, start=None, end=None, sample_rate=FS,
             color='red', marker=None, markersize=6, label=''):
    ''' Adds a peaks array to plt in the range start to end.'''
    if not score:
        return
    if not start:
        start = 0
    if not end:
        end = len(score)
    peak_x = [x for x in score.peaks if x >= start and x < end]
    peak_vals = [score.data[i] for i in peak_x]
    peak_time = []
    for i in peak_x:
        peak_time.append(i / sample_rate)
    plt.plot(peak_time, peak_vals, color = color, linestyle = '',
             marker = marker, markersize = markersize, label = label)
    return len(peak_time)

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

def run_real_time(filename, score_real_time, show_progress = False, 
            write_csv = False, plot_peaks_zoomed = False,
            plot_analysis = True, plot_all_filter_steps = False,
            plot_rr = False, title=''):
    '''
    This version collects cur_bandpass, cur_deriv, etc. separately
    rather than in one cur_filter. It handles mostly plotting and printing.
    The algorithm is run in score_real_time.
    It always plots the peaks.
         show_progress:         Prints information as it processes
         plot_peaks_zoomed:     Plot the peaks zoomed in (first 10 sec)
         plot_analysis:         Creates plots of the steps, first 10 sec, allz
         plot_all_filter_steps: Creates a 6-figure plot of all steps
         write_csv:             Writes a CSV file
    '''
    ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score,\
       ecg_ext, cur_x = score_real_time(filename, show_progress = show_progress)

    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)

    # Plot peaks
    if True:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename
        plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename,
                   title=title)

    # Plot peaks zoomed
    if plot_peaks_zoomed:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename
        for i in range(3):
            plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename,
                xlim=[i * 10, i * 10 + 10])

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
    #if True:
    #    vals1 = np.ndarray.tolist(np.array(avg) * 10.)
    #    label1 = 'Average x 10'
    #    vals2 = score
    #    label2 = 'Score'
    #else:
    #    # DEBUG
    #    vals1 = square
    #    label1 = 'Square'
    #    vals2 = np.ndarray.tolist(np.array(avg) * 10.)
    #    label2 = 'Average x 10'

    ## Plot with specified x axis
    #if plot_analysis:
    #    shift = -18
    #    shifted = shift_score(score, shift)
    #    label2 = f"Score shifted by {shift}"
    #    plot_2_values(ecg_ext, cur_x, vals1, shifted, label1=label1, label2=label2,
    #        title=title, use_time_vals=True, xlim=[0, ZOOMED_PLOT_WIDTH_SEC])

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

def scale1_stats(score_real_time):
    ''' Collects statistics from Scale1 algorithm for all files in the given
   directiory matching the given pattern and writes them to CSV file in the
   given directory with name of the form: Scale1Stats_YYYY-MM-DD_HH-MM.csv.
   '''
    src_dir = r'C:\Scratch\ECG\Polar ECG\CSV'
    pattern = 'PolarECG-2022*.csv'
    dst_dir = r'C:\Scratch\AAA\Python\HeartPy'

    from pathlib import Path
    date = datetime.now().strftime('%Y-%M-%d_%H-%M')
    dst_file = fr'{dst_dir}\Scale1Stats_{date}.csv'
    with open(dst_file, 'w') as fw:
        fw.write('name,length,score_offset,offset_error,offset_delta_avg,'
                    + f'offset_delta_min,offset_delta_max,'
                    + f'offset_delta_stdev,mean,stddev\n')
        for filename in Path(src_dir).glob(pattern):
                    score_real_time(filename, show_progress = False,
                            statistics_file=fw)
    print(f'Wrote {dst_file}')
    fw.close()

def compare_file_peaks(filename, score_real_time, show_segments = False,
        title='Peak Detection Comparison'):
    ecg, file_peaks, _ = ut.read_ecg_file(filename)
    file_peak_indices = []
    for i in range(len(file_peaks)):
        val = file_peaks[i]
        if val > 0:
            file_peak_indices.append(i)

    ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score,\
       ecg_ext, cur_x = score_real_time(filename)

    #npeaks = ut.find_header_item(headers, 'npeaks')
    #print(f'File had {npeaks} peaks found. score_real_time found {len(peak_indices)} peaks.')

    print(f'{len(file_peak_indices)=} {len(peak_indices)=}')
    file_set  = set(file_peak_indices)
    calc_set  = set(peak_indices)
    diffs_a = []
    for item in file_set:
        if item not in calc_set:
            diffs_a.append(item)
    diffs_a.sort()
    diffs_b = []
    for item in calc_set:
        if item not in file_set:
            diffs_b.append(item)
    diffs_b.sort()

    diffs_a_sec = []
    for i in range(len(diffs_a)):
        diffs_a_sec.append(round(diffs_a[i] / FS, 1))

    diffs_b_sec = []
    for i in range(len(diffs_b)):
        diffs_b_sec.append(round(diffs_b[i] / FS, 1))

    print(f'Not in file: {diffs_b}\n    {diffs_b_sec} sec')
    print(f'Not in calculated: {diffs_a}\n    {diffs_a_sec} sec')

    # Plots
    curves = []
    curves.append(curve_data(ecg, file_peak_indices, 'KE.Net ECG Peaks', 'red'))
    curves.append(curve_data(ecg, peak_indices, 'Python Calculated Peaks', 'blue'))
    plot_comparison(curves, filename=filename,
        height_fract = 0, title=title,
        show_segments = show_segments, show_segments_only = False,
        sample_rate = FS)

def analyze(filename, score_real_time, show_progress = False,
                   plot_peaks_zoomed=False,
                   plot_bandpass=False, plot_deriv=False, plot_square=False,
                   plot_avg=False, plot_score=False,
                   do_intervals=False, use_points=False,
                   use_indices=False,
                   title=''):
    ''' Runs score_real_time and generates plots to check for how well the
    algorithm did.
    '''
    ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score,\
       ecg_ext, cur_x = score_real_time(filename, show_progress = show_progress,\
       print_steps=True)
    print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')

    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)

    # Print results
    npeaks = ut.find_header_item(headers, 'npeaks')
    if npeaks:
        print(f'\nFound {len(peak_indices)} peaks. File had {npeaks} peaks.')
    else:
        print(f'\nFound {len(peak_indices)} peaks.')

    # Plot peaks zoomed
    if plot_peaks_zoomed:
        if description:
            peak_filename = f"{filename}\n{description}"
        else:
            peak_filename = filename

        # Values for the original algorithm are extended
        end = len(ecg)
        for i in range(3):
            curves = []
            if plot_bandpass and bandpass and len(bandpass) > 0:
                curves.append(curve_data(bandpass[0:end], None, 'Bandpass', 'dodgerblue'))
            if plot_deriv and deriv and len(deriv) > 0:
                curves.append(curve_data(deriv[0:end], None, 'Deriv', 'gold'))
            if plot_square and square and len(square) > 0:
                curves.append(curve_data(square[0:end], None, 'Square', 'crimson'))
            if plot_avg and avg and len(avg) > 0:
                curves.append(curve_data(avg[0:end], None, 'Avg', 'mediumorchid'))
            if plot_score and score and len(score) > 0:
                curves.append(curve_data(score[0:end], None, 'Score', 'green'))
            if use_indices:
                end_interval = len(ecg) /3
                x_ecg = list(range(len(x_ecg)))
                xlabel = 'sample'
            else:
                end_interval = 10
                xlabel = None

            xlim=[i * end_interval, (i + 1) * end_interval]
            plot_peaks(ecg, x_ecg, peak_indices, filename=peak_filename,
                xlim=xlim,
                use_points=use_points,
                do_intervals=do_intervals,
                curves=curves,
                xlabel=xlabel)

def multiple_compare():
    ''' Compares all the files in the given directory matching the pattern '''
    src_dir = r'C:\Scratch\ECG\Polar ECG\CSV'
    pattern = 'PolarECG-2022*.csv'

    from pathlib import Path
    src_dir = r'C:\Scratch\ECG\Polar ECG\CSV'
    for filename in Path(src_dir).glob(pattern):
        print(filename)
        scores=[]
        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1.name, 'red'))

        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc0.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc0.name, 'dodgerblue'))

        print(filename)
        plot_comparison(scores, filename,
                        height_fract = 0, title='Peak Detection Comparison',
                        show_segments = False, show_segments_only = False,
                        sample_rate = FS)

def do_run_implementation(filename, algorithm, score_real_time,
            do_compare_file_peaks=False,
            do_comparison = False,
            do_comparison1 = False,
            do_multiple_comparison = False,
            do_analyze = False,
            do_scale1_stats = False,
            do_run_real_time = False):
    ''' Implements what is actually being run for the given flags.
    '''

    score_title = f'QRS Algorithm {algorithm} Detected Peaks'
    score_analysis_title = f'QRS Algorithm {algorithm} Detected Peaks Analysis'

    print(f'QRS Algorithm {algorithm}')

    # Print what is being run
    if do_compare_file_peaks:
        print(f'  do compare_file_peaks')
    if do_comparison:
        print(f'  do comparison')
    if do_comparison1:
        print(f'  do comparison1')
    if do_scale1_stats:
        print(f'  do scale1_stats')
    if do_analyze:
        print(f'  do analyze')
    if do_run_real_time:
        print(f'  do run_real_time')

    # Do a comparison
    if do_comparison:
        print('\nRunning comparison')
        scores=[]
        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1g.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1g.name, 'red'))
        print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')

        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1.name, 'dodgerblue'))
        print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')
        
        #ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
        #   = sc2.score_real_time(filename, show_progress = False)
        #scores.append(curve_data(ecg, peak_indices, sc2.name, 'orange'))
        #print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')
        
        #ecg, x_ecg, peak_indices = wav.run(filename, use_ecg=True, plot=False)
        #scores.append(curve_data(ecg, peak_indices, wav.name, 'green'))
        #print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')

        plot_comparison(scores, filename,
                        height_fract = 0, title='Peak Detection Comparison',
                        show_segments = True, show_segments_only = False,
                        sample_rate = FS)

    # Do a comparison of several versions of Score1
    if do_comparison1:
        print('\nRunning comparison1')
        scores=[]

        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1.name, 'red'))
        print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')
        
        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1b.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1b.name, 'green'))
        print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')
        
        ecg, x_ecg, peak_indices, _, _, _, _, _, _, _, _\
           = sc1a.score_real_time(filename, show_progress = False)
        scores.append(curve_data(ecg, peak_indices, sc1a.name, 'dodgerblue'))
        print(f'{len(ecg)=} {len(x_ecg)=} {len(peak_indices)=}')
        
        plot_comparison(scores, filename,
                        height_fract = 0, title='Peak Detection Comparison',
                        show_segments = True, show_segments_only = False,
                        sample_rate = FS)

    if do_multiple_comparison:
        # All files in a directory
        print('\nRunning do_multiple_comparison')
        multiple_compare()
    
    if False:
        print('\nRunning run_process_after')
        # Processing entire file, uses the specified algorithm
        run_process_after(filename)

    if do_scale1_stats:
       print('\nRunning do_scale1_stats')
       scale1_stats(sc1.score_real_time)

    if do_compare_file_peaks:
        print('\nRunning compare_file_peaks')
        compare_file_peaks(filename, score_real_time,
                           title=f'Peak Detection Comparison (Algorithm {algorithm})',
                          show_segments = True),
 
    if do_analyze:
        # Check the file vs processed peaks, uses the specified algorithm
        print('\nRunning analyze')
        analyze(filename, score_real_time, show_progress = False,
                       plot_peaks_zoomed=True,
                       plot_bandpass = False, plot_deriv=True, plot_square=False,
                       plot_avg=True, plot_score=True,
                       do_intervals=True, use_points=True,
                       use_indices=True,
                       title=score_analysis_title)

    if do_run_real_time:
        # Processing as we go, uses the specified algorithm
        print('\nRunning run_real_time')
        #run_real_time(filename, score_real_time, show_progress = True, write_csv = False)
        run_real_time(filename, score_real_time, show_progress = False, # Prints information as it processes
                     plot_all_filter_steps = False, # Creates a 6-figure plot of all steps
                     plot_peaks_zoomed = True, # Plot the peaks zoomed in
                     plot_analysis = False, # Plots of the steps, one zoomed, one not
                     write_csv = False, title=score_title)

def main():
    #global filename
    ut.print_module_filename(__file__)

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
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-06-15_18-30.csv'
        # 12-17-2022 Walking Hillside
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-17_16-22.csv'
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-17_16-16.csv'
        # 12-31-2022 Working on computer. Testing new QRS algorithm.
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-31_16-52.csv'
        # 01-03-2023 Test ECG = sin (2*pi*f), f=12 Hz
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECGTestSineWave.csv'

        # 02-26-2022 Sitting, HR mostly steady
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-26_16-12.csv'
        # 01-3--2022 Test 3, no data
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-01-30_14-10.csv'
        # 01-24-2022 Testing new HR Calculation
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-01-24_15-34.csv'

        # 02-13-2023 First after pacemaker
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-09_14-16.csv'
        # 02-13-2023 Last before pacemaker
        filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-01_16-01.csv'
        # 02-13-2023 New algorithm
        filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2023-02-13_16-27.csv'
        # 02-13-2023 Treadmill, new algorithm
        #filename =  r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2023-02-13_17-31.csv'
        
    #test()
    #test2()

    # Determine scoring algorithm
    algorithm = 1
    if algorithm == 0:
        score_real_time = sc0.score_real_time
    elif algorithm == 1:
        score_real_time = sc1.score_real_time
    elif algorithm == 2:
        score_real_time = sc2.score_real_time

    do_run_implementation(filename, algorithm, score_real_time,
            do_compare_file_peaks = True,   # Compare to file
            do_comparison = False,          # Compare several versions
            do_comparison1 = False,         # Compare several versions (Different score files)
            do_multiple_comparison = False, # Compare all in a directory
            do_analyze = False,              # Compare file vs processed peaks
            do_scale1_stats = False,        # Scale1 stats for all in a directory
            do_run_real_time = False        # Normal run processing
    )

if __name__ == "__main__":
    main()
