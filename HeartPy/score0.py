import numpy as np
import utils as ut
import filter as flt
import process_ecg as pecg
import math
import sys
from filter import Moving_Average

# Common variables
# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0
# Data window size.  Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# HR 200 Interval
HR_200_INTERVAL = int(60.0 / 200.0 * FS)

# Minimum length of score interval
# Determined from mimumum QRS is .08 to .10 sec
MAX_QRS_LENGTH = round(.12 * FS)
# Moving average window size.
MOV_AVG_WINDOW = 20
# Moving average height window size
MOV_AVG_HEIGHT_WINDOW = 5
# Moving average height default
MOV_AVG_HEIGHT_DEFAULT = .025
# Moving average height threshold factor
#     Note: threshold = MOV_AVG_HEIGHT_THRESHOLD_FACTOR * Moving_average.avg()
MOV_AVG_HEIGHT_THRESHOLD_FACTOR = .4

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
    # This is the group delay, used for searching ecg for maxima
    # Only used for scoring
    score_offset = 18
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

