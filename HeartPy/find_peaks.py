'''
Finds peaks using scipy.signal.find_peaks
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
import utils as ut
import process_ecg as pr
import os.path as path
import glob
from pathlib import Path


def plot_peaks(ecg, peaks = None, height_fract = 0,title='Peak Detection',
        show_segments = True, show_segments_only = False,
        sample_rate = 130, do_scipy_peaks = True, peaks_proc = None):
    #ecg = electrocardiogram()[2000:4000]
    # Determine markersizes
    markersize = 6
    peaks_sep = proc_sep = scipy_sep = ''
    if(peaks):
        if(peaks_proc):
            markersize = markersize + 2
            peaks_sep = ', '
            markersize_proc = markersize
            if do_scipy_peaks:
                proc_sep = ', '
                markersize = markersize + 2
                markersize_scipy = markersize
        else:
            if do_scipy_peaks:
                markersize_scipy = markersize + 2
                peaks_sep = ', '
    else:
        if(peaks_proc):
            markersize_proc = markersize
            if do_scipy_peaks:
                proc_sep = ', '
                markersize = markersize + 2
                markersize_scipy = markersize
        else:
            markersize_scipy = markersize

    if not show_segments_only:
        start = 0
        end = len(ecg)
        time = []
        for i in range(start, end):
            time.append(i / sample_rate)
        plt.figure(figsize=(10,6))
        plt.plot(time, ecg)

       # Scipy peaks
        title_scipy = ""
        if do_scipy_peaks:
            maxval = np.max(ecg)
            height = height_fract * maxval
            peaks_scipy, _ = find_peaks(ecg, height=height)
            npeaks = add_peaks(ecg, peaks_scipy, 0, len(ecg), sample_rate,
                color = 'orange', marker = 'o', markersize = markersize_scipy,
                label = 'SciPy Peaks')
            title_scipy = (f'{npeaks} SciPy peaks for height_fract={height_fract}'
                          + scipy_sep)
            print(f'\nfind_peaks: Found {len(peaks_scipy)}' +
                f' for height_fract={height_fract}, height={height}')
       # Process peaks
        title_proc = ""
        if peaks_proc:
            npeaks = add_peaks(ecg, peaks_proc, 0, len(ecg), sample_rate,
                    color = 'b', marker = 'o', markersize = markersize_proc,
                    label = 'Process Peaks')
            title_proc = f'{npeaks} Process peaks' + proc_sep
        # File peaks
        title_peaks = ""
        if peaks:
            peak_vals = [ecg[i] for i in peaks]
            peak_time = []
            for i in peaks:
                peak_time.append(i / sample_rate)
            plt.plot(peak_time, peak_vals, "ro", label = 'Peaks')
            title_peaks = f'{len(peaks)} peaks' + peaks_sep 
        #plt.plot(np.zeros_like(x), "--", color="gray")
        plt.title(title + f"\n{title_peaks}{title_proc}{title_scipy}")
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
        # SciPy peaks
        title_scipy = ""
        if do_scipy_peaks:
            npeaks = add_peaks(ecg, peaks_scipy, start, end, sample_rate,
                color = 'orange', marker = 'o',
                markersize = markersize_scipy, label = 'SciPy Peaks')
            title_scipy = f'{npeaks} of {len(peaks_scipy)} SciPy peaks'\
                f'for height_fract={height_fract}' + scipy_sep
        # Proc peaks
        title_proc = ""
        if peaks_proc:
            npeaks = add_peaks(ecg, peaks_proc, start, end, sample_rate,
               color = 'b', marker = 'o', markersize = markersize_proc,
               label = 'Process Peaks')
            title_proc = f'{npeaks} of {len(peaks_proc)} Process peaks'\
               + proc_sep
        # File peaks
        title_peaks = ""
        if peaks:
            npeaks = add_peaks(ecg, peaks, start, end, sample_rate, color = 'r',
                     marker = 'o', markersize = 6, label = 'Peaks')
            title_peaks = f'{npeaks} of {len(peaks)} peaks' + peaks_sep
        plt.title(title + f"\n{title_peaks}{title_proc}{title_scipy}")
        plt.xlabel(f'time, sec (sample_rate={sample_rate})')
        plt.legend(loc='lower right', framealpha=0.6)
        plt.tight_layout()
        plt.show()

def add_peaks(ecg, peaks, start, end, sample_rate, color, marker, markersize,
        label):
    ''' Adds a peaks array to plt in the range start to end.'''
    peak_x = [x for x in peaks if x >= start and x < end]
    peak_vals = [ecg[i] for i in peak_x]
    peak_time = []
    for i in peak_x:
        peak_time.append(i / sample_rate)
    plt.plot(peak_time, peak_vals, color = color, linestyle = '',
             marker = marker, markersize = markersize, label = label)
    return len(peak_time)

def check_peak_diffs(search_path = r'C:\Scratch\ECG\Polar ECG\CSV',
        csv_name = r'data\peaks-processed-peaks-differences.csv',
        pattern = 'PolarECG-2022', in_both_only = True):
    '''
    Finds all the files in search_path that match pattern
    If write_csv then writes the results to a file.
    If in_both_only only writes the ones that are in both and differ.
   '''
    glob_path = f'{search_path}\{pattern}*.csv'
    print(glob_path)
    files = glob.glob(glob_path)
    #print(f'Found {len(files)} files');
    diffs = []
    for file in files:
        ecg, is_peak, headers = ut.read_ecg_file(file)
        if is_peak:
            peaks = ut.get_peak_values(ecg, is_peak)
        else:
            peaks = None
        _, _, peaks_proc, _, _, _, _, _, _, _, _, _, _, _ = pr.score_real_time(file)
        
        npeaks = 0
        if peaks:
            npeaks = len(peaks)
        nproc = 0
        if peaks_proc:
            nproc = len(peaks_proc)
        not_proc = []
        not_peaks = []
        if peaks and peaks_proc:
            for i in peaks:
                found = False
                for j in peaks_proc:
                    if i == j:
                        found = True
                        break
                if not found:
                    not_proc.append(i)
            for i in peaks_proc:
                found = False
                for j in peaks:
                    if i == j:
                        found = True
                        break
                if not found:
                    not_peaks.append(i)
            if in_both_only:
                if(npeaks > 0 and nproc > 0 and (len(not_proc) > 0 or len(not_peaks) > 0)):
                    diffs.append(f'{Path(file).name}\t{npeaks}\t{nproc}\t{not_proc}\t{not_peaks}')
            else:
                diffs.append(f'{Path(file).name}\t{npeaks}\t{nproc}\t{not_proc}\t{not_peaks}')
        else:
            if not in_both_only:
                diffs.append(f'{Path(file).name}\t{npeaks}\t{nproc}\t{not_proc}\t{not_peaks}')
    ## Print the results
    #for diff in diffs:
    #    print(diff)

    # Write the results to a CSV file
    if csv_name:
        with open(csv_name, 'w') as f:
            timestamp = ut.timestamp('%m/%d/%Y at %H:%M:%S')
            f.write(f'Files with pattern "{pattern}*.csv" that differ in peaks and processed peaks\n')
            f.write(f'Processed on {timestamp}\n')
            for diff in diffs:
                f.write(diff + '\n')
        print(f'Wrote {len(diffs)} lines to {csv_name}')

def check_peaks():
    '''
    Finds all the files in search_path that match pattern
    If write_csv then writes the results to a file.
    If in_both_only only writes the ones that are in both and differ.
   '''
    search_path = r'C:\Scratch\ECG\Polar ECG\CSV'
    filenames = [
        'PolarECG-2022-02-10_13-42.csv',
        'PolarECG-2022-02-10_13-50.csv',
        'PolarECG-2022-02-12_14-44.csv',
        'PolarECG-2022-02-15_12-58.csv',
        'PolarECG-2022-02-15_13-31.csv',
        'PolarECG-2022-02-15_13-57.csv',
        'PolarECG-2022-02-15_14-09.csv',
        'PolarECG-2022-02-18_16-06.csv',
        'PolarECG-2022-02-20_13-44.csv',
        'PolarECG-2022-02-20_15-08.csv',
        'PolarECG-2022-02-26_16-12.csv',
        'PolarECG-2022-03-04_15-39.csv',
        'PolarECG-2022-03-09_12-17.csv',
        'PolarECG-2022-03-14_17-06.csv',
        'PolarECG-2022-03-17_14-23.csv',
        'PolarECG-2022-03-20_17-02.csv',
        'PolarECG-2022-03-20_17-14.csv',
        'PolarECG-2022-03-20_17-25.csv',
        'PolarECG-2022-03-27_13-57.csv',
        'PolarECG-2022-04-02_13-18.csv',
        'PolarECG-2022-04-07_14-59.csv',
        'PolarECG-2022-04-12_13-36.csv',
        'PolarECG-2022-04-20_17-00.csv',
        'PolarECG-2022-05-13_12-42.csv',
        'PolarECG-2022-05-15_11-17.csv'
    ]

    for name in filenames:
        filename = f'{search_path}\\{name}'
        print(f'Processing {filename}')
        # Get peaks from file
        ecg, is_peak, headers = ut.read_ecg_file(filename)
        if is_peak:
            peaks = ut.get_peak_values(ecg, is_peak)
        else:
            peaks = None
            continue

        # Get peaks from process.py
        _, _, peaks_proc, _, _, _, _, _, _, _, _ = pr.score_real_time(filename)

        print(f'{len(peaks)} peaks, {len(peaks_proc)} processed peaks')

        sample_rate = 130.0
        res = ut.find_header_item(headers, 'samplingrate')
        if not res:
            continue
        sample_rate = float(res)
        description = ut.find_header_description(headers)
        if description:
            title = f'Find Peaks\n{name}\n{description}'
        else:
            title = f'Find Peaks\n{name}'

        # Do segments only
        plot_peaks(ecg, peaks = peaks, peaks_proc = peaks_proc, height_fract=.4,
            title= title, sample_rate = sample_rate, do_scipy_peaks = False,
            show_segments = True, show_segments_only = True)

def run():
    print(path.basename(path.normpath(__file__)))

    # Find files where peaks and processed peaks vary and return
    find_peak_differences = False
    if find_peak_differences:
        print('\nFiles with peak vs. processed peaks differences')
        check_peak_diffs(in_both_only = True)
        return

    do_checkPeaks = False
    if do_checkPeaks:
        check_peaks()
        return

    # Set prompt to use default filename or prompt with a FileDialog
    prompt = False
    if prompt:
        file_names = ut.prompt_for_files(\
            title='Choose ECG (Not QRSHR or DeviceHR) file',
            multiple=False, type='csv')
        filename = file_names[0]
    else:
        # 10-31-2021 Working on computer HR=55
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
        # 10-14-21 Walking Falison HR=114
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
        # 10-19-21 Walking Blueberry Lake HR=120
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
        # 2-4-2022 Example New Low Heartrate HR=63
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
        # 2-6-2022 Walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-06_13-52.csv'
        # 2-20-2022 After walking. Device, internal HR differ
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_15-08.csv'
        # 2-20-2022 Before walking 2. Device HR changed from 94 to 54
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-51.csv'
        # 2-20-2022 Before walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-44.csv'
        # 2-3-2022 Sitting down. Little out of breath.
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-03_18-43.csv'
        # 2-15-2022-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-15_12-58.csv'
        # 5-13-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-05-13_12-42.csv'
        # 6-15-18 With new algorithm in KE.Net ECG
        filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-06-15_18-30.csv'

    print(filename)
    ecg, is_peak, headers = ut.read_ecg_file(filename)
    #print(f'{len(ecg)=}')
    #print(f'{len(is_peak)=}')
    # List of ECG values which are peaks
    if is_peak:
        peaks = ut.get_peak_values(ecg, is_peak)
        print(f'Number of peaks: {len(peaks)}')
    else:
        peaks = None
        print(f'Number of peaks: 0')

    # Get peaks from process.py
    _, _, peaks_proc, _, _, _, _, _, _ , _, _ = pr.score_real_time(filename)
    print(f'Number of processed peaks: {len(peaks_proc)}')

    ## Print different indices
    #print('\nPeak vs. processed peaks differences')
    #if not peaks:
    #    print('No file peaks')
    #elif not peaks_proc:
    #    print('No file peaks')
    #else:
    #    print(f'Peaks {len(peaks)} Processed peaks {len(peaks_proc)}')
    #    for i in peaks:
    #        found = False
    #        for j in peaks_proc:
    #            if i == j:
    #                found = True
    #        if not found:
    #            print(f'Not found in processed peaks: {i} ({i / pr.FS:.3f} sec)')
    #    for i in peaks_proc:
    #        found = False
    #        for j in peaks:
    #            if i == j:
    #                found = True
    #        if not found:
    #            print(f'Not found in peaks: {i} ({i / pr.FS:.3f} sec)')

    # Print header
    print()
    sample_rate = 130.0
    for header in headers:
        print(header)
        # Look for the sampling rate
        if header.startswith("samplingrate"):
            try:
                index = header.index("=")
                sample_rate = float(header[index+1: len(header)])
            except:
                None

    description = ut.find_header_description(headers)
    if description:
        title = f'Find Peaks\n{filename}\n{description}'
    else:
        title = f'Find Peaks\n{filename}'

    plot_peaks(ecg, peaks = peaks, peaks_proc = peaks_proc, height_fract=.4,
              title= title, sample_rate = sample_rate)

if __name__ == "__main__":
    run()
