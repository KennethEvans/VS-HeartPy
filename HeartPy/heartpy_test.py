'''
Investigates HeartPy
'''

import heartpy as hp
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import resample
import os.path as path

import heartpy_ke as kehp
import utils as ut

sample_rate = 130

#import sys
#for p in sys.path:
#    print(p)

def plot_data_rr(data, peaklist, rr_masklist, sample_rate, title, do_rejects=False):
    count = len(data)
    peak_count = len(peaklist)
    if peak_count > len(rr_masklist): peak_count = len(rr_masklist)
    x_peak = []
    y_peak = []
    if do_rejects:
        x_peak_reject = []
        y_peak_reject = []
    for i in range(0, peak_count):
        if rr_masklist[i] == 0:
            x_peak.append(peaklist[i] / sample_rate)
            y_peak.append(data[peaklist[i]])
        else:
            if do_rejects:
                x_peak_reject.append(peaklist[i] / sample_rate)
                y_peak_reject.append(data[peaklist[i]])
    x_data = []
    for i in range(0, count):
        x_data.append(i / sample_rate)
    #print()
    #print(len(x_peak))
    #print(len(y_peak))
    #print(len(x_data))
    #print(len(data))
    plt.figure(figsize=(12,4))
    plt.xlabel('Time (s)')
    plt.plot(x_data, data, label='ECG Data');
    plt.plot(x_peak, y_peak, 'go', label='R Peaks')
    if do_rejects:
        plt.plot(x_peak_reject, y_peak_reject, 'ro', label='Rejected Peaks', markersize = 4)
    plt.title(title)
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def plot_data_rr_2(data, peaklist, rr_masklist, sample_rate, title):
    count = len(data)
    peak_count = len(peaklist)
    if peak_count > len(rr_masklist): peak_count = len(rr_masklist)
    x_peak = []
    y_peak = []
    x_peak_reject = []
    y_peak_reject = []
    for i in range(0, peak_count):
        if rr_masklist[i] == 0:
            x_peak.append(peaklist[i] / sample_rate)
            y_peak.append(data[peaklist[i]])
        else:
            x_peak_reject.append(peaklist[i] / sample_rate)
            y_peak_reject.append(data[peaklist[i]])
    x_data = []
    for i in range(0, count):
        x_data.append(i / sample_rate)
    #print()
    #print(len(x_peak))
    #print(len(y_peak))
    #print(len(x_data))
    #print(len(data))
    plt.figure(figsize=(12,4))
    plt.xlabel('Time (s)')
    plt.plot(x_data, data, label='ECG Data');
    plt.plot(x_peak, y_peak, 'go', label='R Peaks')
    plt.plot(x_peak_reject, y_peak_reject, 'ro', label='Rejected Peaks', markersize = 4)
    plt.title(title)
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def fix_rr(rr_list, rr_masklist, peaklist):
    #print(peaklist)
    count = len(rr_list)
    peak_count = len(peaklist)
    new_rr_list = []
    new_rr_masklist = []
    new_peaklist = []
    rejected = False
    for i in range(0, count):
        if rr_masklist[i] == 0:
            new_rr_list.append(rr_list[i])
            if rejected:
                new_peak_count = len(new_peaklist)
                if new_peak_count > 0 and peak_count > i:
                    new_peaklist[new_peak_count-1] = peaklist[i]
            else:
                new_peaklist.append(peaklist[i])
            rejected = False
        else:
            if rejected:
                new_peak_count = len(new_peaklist)
                if new_peak_count > 0 and peak_count > i:
                    new_peaklist[new_peak_count-1] = peaklist[i]
            else :
                new_peaklist.append(peaklist[i])
            rejected = True
    # The last point
    if len(new_peaklist) == count:
        new_peaklist.append(count -1)
    # new_masklist = all 0's
    for i in range(0, len(new_rr_list)): new_rr_masklist.append(0)
    return new_rr_list, new_rr_masklist, new_peaklist

def print_rr(wd, measures, title=None, do_peak_list=True):
    if title: print('\n' + title)
    rr_list = wd['RR_list']
    rr_list_cor = wd['RR_list_cor']
    rr_list_mask = wd['RR_masklist']
    if not 'peaklist' in wd:
        do_peak_list = False
    elif do_peak_list:
        peaklist = wd['peaklist']
    print('rr_list size=%d' % (len(rr_list)))
    print('rr_list_cor size=%d' % (len(rr_list_cor)))
    print('rr_list_mask size=%d' % (len(rr_list_mask)))
    print('peaklist size=%d' % (len(peaklist)))
    count = len(rr_list)
    print('%d RR values of which %d were rejected and %d kept:'
         % (count, count - len(rr_list_cor), len(rr_list_cor)))
    sample_rate = wd['sample_rate']
    sum = np.sum(rr_list) / 1000 # (in sec)
    mean_rr = np.mean(rr_list)
    hr = (60000 / mean_rr) if mean_rr != 0 else 0
    print('  sample_rate=%d' % (sample_rate))
    print('  mean_rrr=%.2f sum(sec)=%.2f hr=%.2f' % (mean_rr, sum, hr))
    mean_rr_cor = kehp.mean_rr_cor(rr_list, rr_list_mask)
    sum_cor = np.sum(rr_list_cor) / 1000  #(in sec)
    hr_cor = (60000 / mean_rr_cor) if mean_rr_cor != 0 else 0
    print('  mean_rr_cor=%.2f sum(sec)=%.2f hr=%.2f' % (mean_rr_cor, sum_cor, hr_cor))
    print('  measures hr=%.2f' % (measures['bpm']))
    for i in range(0, count):
        rr = rr_list[i]
        if do_peak_list:
            start = int(peaklist[i])
            end = int(peaklist[i+1])
            delta = end - start
            if(rr_list_mask[i] == 0):
                print('%2d %8.2f %d to %d (%d)' % (i + 1, rr, start, end, delta))
            else:
                print('%2d %8.2f %d to %d (%d) %s' % (i + 1, rr, start, end, delta, 'Rejected'))
        else:
            if(rr_list_mask[i] == 0):
              print('%2d %8.2f' % (i + 1, rr))
            else:
             print('%2d %8.2f %s' % (i + 1, rr, 'Rejected'))

def test1():
    #first let's load the clean PPG signal
    data, timer = hp.load_exampledata(0)
    plt.figure(figsize=(12,4))
    plt.plot(data)
    plt.show()

def test2():
    #first let's load the clean PPG signal
    data, timer = hp.load_exampledata(0)
    #run the analysis
    wd, m = hp.process(data, sample_rate = 100.0)

    #display measures computed
    for measure in m.keys():
        print('%s: %f' %(measure, m[measure]))

    #call plotter
    hp.plotter(wd, m, figsize=(12,4))
    input('Press Enter to continue...')   

def test3():
    # Working on computer HR=55
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
    # Walking Falison HR=114
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
    # Walking Blueberry Lake HR=120
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
    # Feb 4 Example New Low Heartrate HR=63
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
    data, _, headers = ut.read_ecg_file(filename)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)

    # Note first two arguments are required
    wd, m = hp.process(data, sample_rate, hampel_correct=False)

    #Debugging
    if True:
        print(f"best={wd['best']}")
        print(f"rrsd={wd['rrsd']}")

    # Print the keys
    if False:
        print('\nwd keys:')
        for key in wd.keys():
            print(key)
        print('\nm keys:')
        for key in m.keys():
            print(key)

    # Scale
    if False:
        wd, m = hp.process(hp.scale_data(data), sample_rate)

    # Resample
    if False:
        #resample the data. Usually 2, 4, or 6 times is enough depending on original sampling rate
        resample_factor = 2
        resampled_data = resample(data, len(data) * resample_factor)

        #And run the analysis again. Don't forget to up the sample rate as well!
        wd, m = hp.process(hp.scale_data(resampled_data), sample_rate * 2)
        #wd, m = hp.process(data, sample_rate * resample_factor)

    # Display measures computed
    if True:
        print('\nMeasure Information:')
        for measure in m.keys():
            print('%s: %f' %(measure, m[measure]))

    # Process RR separately using KE version
    if True:
        rr_list_1 = wd['RR_list']
        print('........')
        wd2, m2 = kehp.process_rr(rr_list_1, threshold_rr=True, clean_rr=True, working_data=wd)
        print('........')
        # Print the keys
        if False:
            print('\nwd keys:')
            for key in wd2.keys():
                print(key)
        print_rr(wd2, m2, title='RR from kehp.process_rr')
        if True:
            # Plot
            plot_data_rr(data, wd2['peaklist'], wd2['RR_masklist'], sample_rate, 
                           'KE Modified ECG Data and R Peak Data\n' + filename, do_rejects=False)

    # Fix RR
    if True:
        new_rr_list, new_rr_masklist, new_peaklist = fix_rr(wd['RR_list'], wd['RR_masklist'] ,wd['peaklist'])
        wd3 = {}
        wd3['RR_list'] = new_rr_list
        wd3['RR_masklist'] = new_rr_masklist
        wd3['RR_list_cor'] = new_rr_masklist
        # TODO
        new_peaklist.append(0)
        wd3['peaklist'] = new_peaklist
        wd3['sample_rate'] = sample_rate
        print_rr(wd3, m2, title='RR from kehp.process_rr after fix')
        #print()
        #print(len(data))
        #print(len(new_peaklist))
        #print(len(new_rr_list))

        # Plot
        plot_data_rr(data, new_peaklist, new_rr_masklist, sample_rate, 
                       'Original ECG Data and R Peak Data\n' + filename, do_rejects=False)
  
    # Display RR values
    if True:
        print_rr(wd, m, title='RR from hp.process')

    # Plot data and all peaks
    plot_data_rr(data, wd['peaklist'], wd['RR_masklist'], sample_rate, 
                   'Original ECG Data and Accepted and Rejected R Peak Data\n' + filename, do_rejects=True)

    # Call plotter
    if False:
        hp.plotter(wd, m, figsize=(12,4),
                  title='Heart Rate Signal Peak Detection\n' + filename,
                  moving_average=True)

    input('\nPress Enter to continue...')   

def main():
    print(path.basename(path.normpath(__file__)))
    test3()

if __name__ == "__main__":
    main()
