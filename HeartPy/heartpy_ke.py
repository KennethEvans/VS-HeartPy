'''
Has functions from HeartPy that have been modified to function differently
    process_rr
    mean_rr_cor (not modified)
'''

from datetime import datetime
import time
import os
import sys

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.signal import butter, filtfilt, welch, periodogram, resample_poly, resample

#from . import exceptions
#from .datautils import get_data, get_samplerate_mstimer, get_samplerate_datetime,\
#                       rolling_mean, outliers_iqr_method, outliers_modified_z, \
#                       load_exampledata
#from .preprocessing import scale_data, scale_sections, interpolate_clipping, \
#                           flip_signal, enhance_peaks, enhance_ecg_peaks
#from .filtering import filter_signal, hampel_filter, hampel_correcter, \
#                       remove_baseline_wander, smooth_signal
#from .peakdetection import make_windows, append_dict, fit_peaks, check_peaks, \
#                           check_binary_quality, interpolate_peaks
#from .visualizeutils import plotter, segment_plotter, plot_poincare, plot_breathing
from heartpy.analysis import calc_rr, calc_rr_segment, clean_rr_intervals, calc_ts_measures, \
                      calc_fd_measures, calc_breathing, calc_poincare

#from . import config
#config.init() #initialize global conf vars

def mean_rr_cor(rr_list, rr_list_mask):
    sum = 0
    n_cor = 0;
    count = len(rr_list)
    for i in range(0, count):
        sum = sum + 1
        rr = rr_list[i]
        if(rr_list_mask[i] == 0):
            n_cor = n_cor +1
            sum = sum + rr
    if(n_cor > 0):
       mean = sum/n_cor
    else:
       mean = 0
    return mean

def process_rr(rr_list, threshold_rr=False, clean_rr=False,
               clean_rr_method='quotient-filter', calc_freq=False,
               freq_method='welch', welch_wsize=240, square_spectrum=True,
               breathing_method='welch', measures={}, working_data={}):
    '''process rr-list

    Function that takes and processes a list of peak-peak intervals (tachogram).
    Computes all measures as computed by the regular process() function, and
    sets up all dicts required for plotting poincare plots.

    Note: This method assumes ms-based tachogram

    Several filtering methods are available as well.

    Parameters
    ----------
    rr_list : 1d array or list
        list or array containing peak-peak intervals (in ms).

    threshold_rr : bool
        if true, the peak-peak intervals are cleaned using a threshold filter, which
        rejects all intervals that differ 30% from the mean peak-peak interval, with
        a minimum of 300ms.
        default : false

    clean_rr : bool
        if true, the RR_list is further cleaned with an outlier rejection pass. This pass
        is performed after threshold_rr, if that is specified.
        default : false

    clean_rr_method: str
        how to find and reject outliers. Available methods are 'quotient-filter',
        'iqr' (interquartile range), and 'z-score'.
        default : 'quotient-filter'

    calc_freq : bool
        whether to compute time-series measurements
        default : False

    freq_method : str
        method used to extract the frequency spectrum. Available: 'fft' (Fourier Analysis),
        'periodogram', and 'welch' (Welch's method).
        default : 'welch'

    welch_wsize : float
        Size of window (in sec) when welch method used to compute the spectrogram.
        This choice is based on a trade-off btw temporal res and freq res of the resulting spectrum
        60 sec may seem reasonable, but this would greatly limit frequency resolution!
          1/60 s = 0.017 Hz, so you would only have 2 points in the VLF band
        Therefore, the default is 4 min (9 points in the VLF band)
        default : 240

    square_spectrum : bool
        whether to square the power spectrum returned.
        default : true

    measures : dict
        dictionary object used by heartpy to store computed measures. Will be created
        if not passed to function.

    working_data : dict
        dictionary object that contains all heartpy's working data (temp) objects.
        will be created if not passed to function

    Returns
    -------
    working_data : dict
        dictionary object used to store temporary values.

    measures : dict
        dictionary object used by heartpy to store computed measures.

    Examples
    --------
    Let's generate an RR-list first.

    >>> import heartpy as hp
    >>> data, timer = hp.load_exampledata(2)
    >>> sample_rate = hp.get_samplerate_datetime(timer, timeformat = '%Y-%m-%d %H:%M:%S.%f')
    >>> wd, m = hp.process(data, sample_rate)
    >>> rr_list = wd['RR_list']

    Using only the RR-list (in ms!) we can now call this function, and let's put the results
    into a differently named container so we're sure all measures are unique:
    >>> wd2, m2 = process_rr(rr_list, threshold_rr = True, clean_rr = True, calc_freq = True)
    >>> '%.3f' %m2['rmssd']
    '45.641'

    If you want to, you can turn off all filters and rejection features:
    >>> wd2, m2 = process_rr(rr_list, threshold_rr = False, clean_rr = False)
    >>> '%.3f' %m2['rmssd']
    '162.645'

    In this case it seems the filtering was necessary: without the RMSSD lies outside the
    range expected in healthy humans.
    '''

    working_data['RR_list'] = rr_list

    if threshold_rr:
        #do thresholding pass
        mean_rr = np.mean(rr_list)
        #upper_threshold = mean_rr + 300 if (0.3 * mean_rr) <= 300 else mean_rr + (0.3 * mean_rr)
        #lower_threshold = mean_rr - 300 if (0.3 * mean_rr) <= 300 else mean_rr - (0.3 * mean_rr)
        upper_threshold = 3333
        lower_threshold = 300
        rr_list_cor = [x for x in rr_list if x > lower_threshold and x < upper_threshold]
        rr_mask = [1 if x <= lower_threshold or x >= upper_threshold else 0 for x in rr_list]
        working_data['RR_list_cor'] = rr_list_cor
        working_data['RR_masklist'] = rr_mask

    if clean_rr:
        #do clean_rr pass
        working_data = clean_rr_intervals(working_data = working_data, method = clean_rr_method)

    if not threshold_rr and not clean_rr:
        working_data['RR_list_cor'] = rr_list
        working_data['RR_masklist'] = [0 for i in range(len(rr_list))]
        rr_diff = np.abs(np.diff(rr_list))
        rr_sqdiff = np.power(rr_diff, 2)
    else:
        rr_diff = np.abs(np.diff(working_data['RR_list_cor']))
        rr_sqdiff = np.power(rr_diff, 2)

    #compute ts measures
    working_data, measures = calc_ts_measures(rr_list = working_data['RR_list_cor'], rr_diff = rr_diff,
                                              rr_sqdiff = rr_sqdiff, measures = measures,
                                              working_data = working_data)

    measures = calc_poincare(rr_list = working_data['RR_list'], rr_mask = working_data['RR_masklist'],
                             measures = measures, working_data = working_data)
    if calc_freq:
        #compute freq measures
        working_data, measures = calc_fd_measures(method=freq_method, welch_wsize=240, square_spectrum=square_spectrum,
                                                  measures=measures, working_data=working_data)

    #compute breathing
    try:
        measures, working_data = calc_breathing(working_data['RR_list_cor'], method = breathing_method,
                                                measures=measures, working_data=working_data)
    except:
        measures['breathingrate'] = np.nan

    return working_data, measures

