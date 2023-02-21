import utils as ut
import filter as flt
import process_ecg as pecg
import math
import sys

name = 'QRS Detection Jan 2023 [F]'

# Common variables
# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0
# Data window size.  Must be large enough for maximum number of coefficients.
DATA_WINDOW = 20
# HR 200 Interval
HR_200_INTERVAL = int(60.0 / 200.0 * FS)
# Expected signal delay (Butterworth ~7, Derivative 4, Square 1)
SIGNAL_DELAY = 12

def score_real_time(filename, show_progress = False, statistics_file=None,
                   print_steps=False):
    '''
    This is the scoring part of run_real_time.
    Returns ecg, x_ecg, peak_indices, headers, avg, score, ecg_ext, cur_x
    show_progress: Printout the steps taken at the end of each interval
    statistics_file: A file object for writing CSV statistics about the offsets
    used at each interval.
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
    score_offset = 0
    # Look around score_offset by =/- this much to find true maximum
    offset_error = 2
    # Peak is past point of max derivitive
    # Extend interval for looking for peak by this amount
    search_extend = 2

    # Multiplier for mean + n_stdev * stddev in threshold criteria
    n_sigma = 2.5
  
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

    # Initialize these with observed values
    mean = stat_initial_mean = .00
    stddev = stat_initial_stddev = .1
    threshold = mean + n_sigma * stddev
    sumvals = stat_initial_mean * HR_200_INTERVAL
    sumsq = (stat_initial_stddev * stat_initial_stddev -
                  stat_initial_mean * stat_initial_mean) * HR_200_INTERVAL
    n_stat = HR_200_INTERVAL

    # Set up usage statistics
    if statistics_file:
        offset_delta_n = 0
        offset_delta_sum = 0
        offset_delta_sumsq = 0
        offset_delta_max = 0
        offset_delta_min = 0

    peak_index = -1
    ecg_max_val = -sys.float_info.max
    # Loop over ECG values
    for i in range(0, n_ecgext):
        if len(cur_ecg) == keep:
            cur_ecg.pop(0)
        cur_ecg.append(ecg_ext[i])
        #print(f"{i} {len(cur_ecg)=}")

        # Butterworth (!!!! not using)
        input = cur_ecg
        if len(cur_butterworth) == keep:
             cur_butterworth.pop(0)
        cur_butterworth.append(0) # Doesn't matter
        new = flt.butterworth3(input, cur_butterworth)
        cur_butterworth[-1] = new
        bandpass.append(new)

        # Derivative and Square together
        input = cur_ecg
        if len(cur_deriv) == keep:
             cur_deriv.pop(0)
        cur_deriv.append(0) # Doesn't matter
        if len(cur_square) == keep:
             cur_square.pop(0)
        cur_square.append(0) # Doesn't matter
        new = flt.derivative5(input, cur_deriv)
        cur_deriv[-1] = new
        deriv.append(new)
        newsquare = max(new, 0)
        cur_square[-1] = newsquare
        square.append(newsquare)
        n_stat = n_stat + 1
        sumvals = sumvals + newsquare
        sumsq = sumsq + newsquare * newsquare
        mean = sumvals / n_stat
        variance = sumsq/ n_stat + mean*mean
        stddev = math.sqrt(variance)
        threshold = mean + n_sigma * stddev
        score.append(threshold)
    
        ## Square
        #input = cur_deriv
        #if len(cur_square) == keep:
        #     cur_square.pop(0)
        #cur_square.append(0) # Doesn't matter
        #new = flt.square(input, cur_square)
        #cur_square[-1] = new
        #square.append(new)

        # Not used
        #score = None # Not using score
        avg = None # Not using avg

        input = cur_square

        # Process finding the peaks
        if i % HR_200_INTERVAL == 0 or i == n_ecgext - 1:
            # End of interval, process this interval
            if i > 0 and peak_index != -1:
                # Look between min_peak_index and max_peak_index for maximum ecg value
                ecg_max_val = -sys.float_info.max
                start_search = max(i - HR_200_INTERVAL, min_peak_index)
                if start_search < 0: start_search = 0;
                end_search = min(i, max_peak_index + search_extend)
                if end_search > len(ecg_ext): end_search = len(ecg_ext)
                for i1 in range(start_search, end_search + 1):
                    if ecg_ext[i1] > ecg_max_val:
                        ecg_max_val = ecg_ext[i1]
                        peak_index = i1
                if print_steps:
                    print(f'Peak  at {peak_index},'
                        + f' searched {start_search} to {end_search},'
                        + f' found ecg_max_val at {peak_index}')
                # Check if there is a close one in the previous interval
                if len(peak_indices) > 0:
                    last_index = len(peak_indices) - 1 # last index in peak_indices
                    last_peak_index = peak_indices[last_index]
                    if peak_index - last_peak_index < HR_200_INTERVAL:
                        last_max_ecg_val = ecg_ext[last_peak_index]
                        if ecg_max_val >= last_max_ecg_val:
                            # Replace the old one
                            peak_indices[last_index] = peak_index
                            if print_steps:
                                print(f'Near a previous one,'
                                    + f' replaced ecg_max_val at {last_peak_index}'
                                    + f' with new one at {peak_index}')
                    else:
                        # Is not near a previous one, add it
                        peak_indices.append(peak_index)
                        if print_steps:
                            print(f'Not near a previous one,'
                                + f' added ecg_max_val at {peak_index}')
                else:
                    # First peak
                    peak_indices.append(peak_index)
                    if print_steps:
                        print(f'First peak, added ecg_max_val at {peak_index}')
                
            # Start a new interval
            peak_index = -1
            min_peak_index = -1
            max_peak_index = -1
            ecg_max_val = -sys.float_info.max

        # Accumulate statistics
        val = input[-1]
        scoreval = score[-1]
        # Note: This gets end values if i < score_offset, these are 0
        if i >= score_offset:
            ecg_val = ecg_ext[i - score_offset]
            #n_stat = n_stat + 1
            #sumvals = sumvals + val
            #sumsq = sumsq + val * val
            #mean = sumvals / n_stat
            #variance = sumsq/ n_stat + mean*mean
            #stddev = math.sqrt(variance)
            #threshold = mean + n_stdev * stddev
            if val > scoreval:
                peak_index = i - score_offset
                if peak_index > max_peak_index:
                    max_peak_index = peak_index
                if peak_index < min_peak_index:
                    min_peak_index = peak_index

        # Print running mean and stddev values
        if False:
            interval = 500
            if i < 80 or i % interval == 0 or i == n_ecgext-1:
                print(f'{i=} {n_stat=} {mean=} {stddev=}')

        #if i < HR_200_INTERVAL:
        #    if ecg_max_val == -sys.float_info.max:
        #        print(f'{i=} {val=:0.3f} {peak_index=} ecg_max_val=undef')
        #    else:
        #        print(f'{i=} {val=:0.3f} {peak_index=} {ecg_max_val=:.3f}')

    if statistics_file:
        if offset_delta_n:
            offset_delta_avg = offset_delta_sum / offset_delta_n
            offset_delta_stdev = math.sqrt(offset_delta_sumsq / offset_delta_n +\
                offset_delta_avg * offset_delta_avg)
        else:
            offset_delta_avg = offset_delta_stdev = 0
        statistics_file.write(f'{filename},{len(ecg)},{score_offset},{offset_error},'
                    + f'{offset_delta_avg:.3f},{offset_delta_min},'
                    + f'{offset_delta_max},{offset_delta_stdev:.3f},'
                    + f'{mean:.3f},{stddev:.3f}\n')
    if True:
        print(f'{n_stat=} {mean=:.3f} {stddev=:.3f} {threshold=:.3f}')

    return ecg, x_ecg, peak_indices, headers, bandpass, deriv, square, avg, score, ecg_ext, cur_x

