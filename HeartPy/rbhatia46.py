'''
Implementation of https://github.com/rbhatia46/Pan-Tompkins-Algorithm-Python
without using Pandas.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import csv
import utils as ut

def plot_all(ecg, f1, f2, f3, f4, f5, title='Pan-Tompkins'):
    plt.figure(figsize=(12,7))
    plt.suptitle(title)

    plt.subplot(2, 3, 1)
    plt.title('ECG')
    plt.plot(ecg)

    plt.subplot(2, 3, 2)
    plt.title('f1')
    plt.plot(f1)

    plt.subplot(2, 3, 3)
    plt.title('f2')
    plt.plot(f2)

    plt.subplot(2, 3, 4)
    plt.title('f3')
    plt.plot(f3)

    plt.subplot(2, 3, 5)
    plt.title('f4')
    plt.plot(f4)

    plt.subplot(2, 3, 6)
    plt.title('f5')
    plt.plot(f5)

    plt.show()

def plot(ecg, timevals, label='Data', title='Pan-Tompkins Algorithm', use_time_vals=True):
    plt.figure(figsize=(12,7))
    if  not timevals: use_time_vals = False
    if use_time_vals == False:
        plt.plot(ecg, label=label)
    else:
        plt.plot(timevals, ecg, label=label)
    plt.title(title)
    #plt.xlabel('time')
    #plt.ylabel('mV')
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def read_csv_file(filename):
    data = []
    timevals = []
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                timevals.append(float(row[0]))
                data.append(float(row[1]))
                line_count += 1
    return timevals, data

#low-pass filter
def lpf(x):
    y = x.copy()
    count = len(x)
    
    for n in range(0, count):
        if(n < 12):
            continue
        y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[n-6] + x[n-12] 
    return y


#high-pass filter
def hpf(x):
    y = x.copy()
    count = len(x)
    
    for n in range(0, count):
        if(n < 32):
            continue
        y[n] = y[n-1] - x[n]/32 + x[n-16] - x[n-17] + x[n-32]/32
    return y

#derivative of signal
def deriv(x):
    y = x.copy()
    count = len(x)

    for n in range(0, count):
        if(n < 4):
            continue
        y[n] = (2*x[n] + x[n-1] - x[n-3] - 2*x[n-4])/4
    return y

#squaring the signal
def squaring(x):
    y = x.copy()
    count = len(x)

    for n in range(0, count):
        y[n] = x[n]**2
    return y

#integral of the signal for a moving window of ws size.
def win_sum(x, ws):
    y = x.copy()
    l = int(ws/2)
    count = len(x)
    
    for n in range(0, count):
        tmp_sum = 0
        
        if(n > 933-l):
            break

        if(n < l):
            continue
        for j in range(n-l,n+l+1):
            tmp_sum += x[j]
        y[n] = tmp_sum/(l+1)		
    return y

# Not used
def detection(x):
    y = x.copy()

def run():
    print('pan-tompkins.py\n');

    if True:
        # Feb 4 Example New Low Heartrate HR=63
        filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
        ecg, _, headers = ut.read_ecg_file(filename)

    else:
        filename = r'data\ecg_data.csv'
        timevals, ecg = read_csv_file(filename)

    if True:
        try:
            print('pan-tompkins.py')
            print(f"type(ecg)={type(ecg)}")
            ecglen = len(ecg);
            print(f"len(ecg)={len(ecg)}")
            delta_time=timevals[ecglen - 1] - timevals[0]
            print(f"delta_time={delta_time}")
            print(f"time_step={delta_time / ecglen}")
            print(f"sampling_rate assuming 5 sec={ecglen / 5}")
            print(f"time assuming sampling_rate is 200 Hx={ecglen / 200}")
        except:
            pass

    #Application of lpf
    f1 = lpf(ecg)
    #Application of hpf
    f2 = hpf(f1)
    #Application of the derivative
    f3 = deriv(f2)
    #squaring signal
    f4 = squaring(f3)

    window_size = 22 
    f5 = win_sum(f4, window_size)

    #plot(ecg, timevals, use_time_vals=False)
    plot(f5, None, title='Pan-Tompkins\n' + filename, use_time_vals=False)

    plot_all(ecg, f1, f2, f3, f4, f5, title='Pan-Tompkins\n' + filename)

def main():
    run()

if __name__ == "__main__":
    main()

