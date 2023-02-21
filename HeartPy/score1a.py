'''
Score1 as used for 12-30-22 changes to KE.Net
'''

import numpy as np
import utils as ut
import filter as flt
import process_ecg as pecg
import math
import sys

name = 'QRS Detection Dec 2022'

# Common variables
# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0
# Data window size.  Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# HR 200 Interval
HR_200_INTERVAL = int(60.0 / 200.0 * FS)
# Expected signal delay (Butterworth ~7, Derivative 4, Square 1)
SIGNAL_DELAY = 12

def score_real_time(filename, show_progress = False, fw=None):
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
    # This is the group delay, used for searching ecg for maxima
    # Only used for scoring
    score_offset = 8
    # Look around score_offset by =/- this much to find true maximum
    offset_error = 5

    # Multiplier for mean + n_stdev * stddev in threshold criteria
    n_stdev = 2
  
    peaks = []
    peak_indices = []
    max_avg = []
    max_avg_indices = []
    scoring = False
    score_start = 0
    score_end = 0
    max_index = min_index = -1

   # Create an extended array with zeros in the last score_offset places
    ecg_ext = ecg.copy() + [0.0] * score_offset
    n_ecgext = len(ecg_ext)
    x_ecg = [i / FS for i in range(necg)]
    cur_x = [i / FS for i in range(n_ecgext)]
   
    keep = DATA_WINDOW # Keep this many items

    ## Keep a moving average of the moving average heights
    ## Initialize it assuming avg peaks are MOV_AVG_HEIGHT_DEFAULT high
    #moving_average_height_list = [MOV_AVG_HEIGHT_DEFAULT] * MOV_AVG_HEIGHT_WINDOW
    #moving_average_height = Moving_Average(MOV_AVG_HEIGHT_WINDOW,
    #    moving_average_height_list)
    #threshold = MOV_AVG_HEIGHT_THRESHOLD_FACTOR * moving_average_height.avg()
    ##print(f'Starting {threshold=} {moving_average_height.avg()=}')

    # Initialize these with observed values
    stat_initial_mean = .00
    stat_initial_stddev = .03
    sumvals = stat_initial_mean * HR_200_INTERVAL
    sumsq = (stat_initial_stddev * stat_initial_stddev -
                  stat_initial_mean * stat_initial_mean) * HR_200_INTERVAL
    n_stat = HR_200_INTERVAL

    # Set up usage statistics
    offset_delta_n = 0
    offset_delta_sum = 0
    offset_delta_sumsq = 0
    offset_delta_max = 0
    offset_delta_min = 0

    peak_index = -1
    max_val = -sys.float_info.max
    # Loop over ECG values
    for i in range(0, n_ecgext):
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

        ## Moving average
        #input = cur_square
        #if len(cur_avg) == keep:
        #     cur_avg.pop(0)
        #cur_avg.append(0) # Doesn't matter
        #new = flt.moving_average1(input, MOV_AVG_WINDOW)
        #cur_avg[-1] = new
        #avg.append(new)

        # Not used
        score = [0.] * n_ecgext # Not using score
        avg = [0.] * n_ecgext # Not using avg

        input = cur_square

        # Process finding the peaks
        if i % HR_200_INTERVAL == 0 or i == n_ecgext - 1:
            # End of interval, process this interval
            if i > 0 and peak_index != -1:

                # Check if this actually at the local maximum
                local_ecg_max = max_ecg_val
                new_peak_index = peak_index
                for i1 in range(peak_index - offset_error, peak_index + offset_error + 1):
                    if ecg_ext[i1] > local_ecg_max:
                        local_ecg_max = ecg_ext[i1]
                        new_peak_index = i1
                # Offset statistics
                offset_delta_n = offset_delta_n + 1
                offset = peak_index - new_peak_index
                offset_delta_sum = offset_delta_sum + offset
                offset_delta_sumsq = offset_delta_sumsq + offset * offset
                if offset > offset_delta_max:
                    offset_delta_max = offset
                if offset < offset_delta_min:
                    offset_delta_min = offset
                # Use the new peak_index
                peak_index = new_peak_index
                # Check if there is a close one in the previous interval
                if len(peak_indices) > 0:
                    last_index = len(peak_indices) - 1 # last index in peak_indices
                    last_peak_index = peak_indices[last_index]
                    if peak_index - last_peak_index < HR_200_INTERVAL:
                        last_max_ecg_val = ecg_ext[last_peak_index]
                        if max_ecg_val >= last_max_ecg_val:
                            peak_indices[last_index] = peak_index
                    else:
                        peak_indices.append(peak_index)
                else:
                    peak_indices.append(peak_index)
                
            # Start a new interval
            peak_index = -1
            max_ecg_val = -sys.float_info.max

        # Accumulate statistics
        val = input[-1]
        ecgval = ecg_ext[i - score_offset]
        n_stat = n_stat + 1
        sumvals = sumvals + val
        sumsq = sumsq + val * val
        mean = sumvals / n_stat
        variance = sumsq/ n_stat + mean*mean
        stddev = math.sqrt(variance)
        if val > mean + n_stdev * stddev:
            if ecgval >= max_ecg_val:
                max_ecg_val = ecg_ext[i - score_offset]
                peak_index = i - score_offset

        # Print mean and stddev values
        if False:
            interval = 500
            if i < 80 or i % interval == 0 or i == n_ecgext-1:
                print(f'{i=} {n_stat=} {mean=} {stddev=}')

        #if i < HR_200_INTERVAL:
        #    if max_val == -sys.float_info.max:
        #        print(f'{i=} {val=:0.3f} {peak_index=} max_val=undef')
        #    else:
        #        print(f'{i=} {val=:0.3f} {peak_index=} {max_val=:.3f}')

    if fw:
        if offset_delta_n:
            offset_delta_avg = offset_delta_sum / offset_delta_n
            offset_delta_stdev = math.sqrt(offset_delta_sumsq / offset_delta_n +\
                offset_delta_avg * offset_delta_avg)
        else:
            offset_delta_avg = offset_delta_stdev = 0
        fw.write(f'{filename},{len(ecg)},{score_offset},{offset_error},'
                    + f'{offset_delta_avg:.3f},{offset_delta_min},'
                    + f'{offset_delta_max},{offset_delta_stdev:.3f},'
                    + f'{mean:.3f},{stddev:.3f}\n')
    if True:
        print(f'{n_stat=} {mean=} {stddev=}')

    return ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score, ecg_ext, cur_x

