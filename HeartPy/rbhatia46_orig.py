import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as sig

def plot_all(timevals, ecg, f1, f2, f3, f4, f5):
    plt.figure(figsize=(12,7))

    plt.subplot(2, 3, 1)
    plt.title('ECG')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.subplot(2, 3, 2)
    plt.title('f1')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.subplot(2, 3, 3)
    plt.title('f2')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.subplot(2, 3, 4)
    plt.title('f3')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.subplot(2, 3, 5)
    plt.title(timevals, 'f4')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.subplot(2, 3, 6)
    plt.title('f5')
    plt.plot(f5.iloc[:,0], f5.iloc[:,1])
    plt.title('Pan-Tompkins algorithm')
    plt.xlabel('time')
    plt.ylabel('mV')

    plt.suptitle('Pan-Tompkins')
    plt.show()

def plot(ecg, timevals, label='Data', title='Pan-Tompkins Algorithm', use_time_vals=True):
    if use_time_vals == 0:
        plt.plot(ecg, label=label)
    else:
        plt.plot(timevals, ecg, label=label)
    plt.title(title)
    #plt.xlabel('time')
    #plt.ylabel('mV')
    plt.legend(loc=4, framealpha=0.6)
    plt.show()



#low-pass filter
def lpf(x):
    y = x.copy()
    
    for n in x.index:
        if(n < 12):
            continue
        y.iloc[n,1] = 2*y.iloc[n-1,1] - y.iloc[n-2,1] + x.iloc[n,1] - 2*x.iloc[n-6,1] + x.iloc[n-12,1] 
    return y


#high-pass filter
def hpf(x):
    y = x.copy()
    
    for n in x.index:
        if(n < 32):
            continue
        y.iloc[n,1] = y.iloc[n-1,1] - x.iloc[n,1]/32 + x.iloc[n-16,1] - x.iloc[n-17,1] + x.iloc[n-32,1]/32
    return y

#defivative of signal
def deriv(x):
    y = x.copy()

    for n in x.index:
        if(n < 4):
            continue
        y.iloc[n, 1] = (2*x.iloc[n,1] + x.iloc[n-1,1] - x.iloc[n-3,1] - 2*x.iloc[n-4,1])/4
        y[1] =         (2*x[n]        + x[n-1]        - x[n-3]        - 2*x[n-4])/4
    return y

#squarring the signal
def squaring(x):
    y = x.copy()

    for n in x.index:
        y.iloc[n,1] = x.iloc[n,1]**2
    return y

#integral of the signal for a moving window of ws size.
def win_sum(x, ws):
    y = x.copy()
    l = int(ws/2)
    
    for n in x.index:
        tmp_sum = 0
        
        if(n > 933-l):
            break

        if(n < l):
            continue
        for j in range(n-l,n+l+1):
            tmp_sum += x.iloc[j,1]
        y.iloc[n,1] = tmp_sum/(l+1)		
    return y

def detection(x):
    y = x.copy()

ecg = pd.read_csv(r'data\ecg_data.csv')

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

plt.figure(figsize=(12,7))
plt.suptitle('Pan-Tompkins Algorithm')

plt.subplot(2, 3, 1)
plt.plot(ecg.iloc[:,0], ecg.iloc[:,1])
plt.title('ecg')
#plt.xlabel('time')
#plt.ylabel('mV')

plt.subplot(2, 3, 2)
plt.plot(f1.iloc[:,0], f1.iloc[:,1])
plt.title('f1')
#plt.xlabel('time')
#plt.ylabel('mV')

plt.subplot(2, 3, 3)
plt.plot(f2.iloc[:,0], f2.iloc[:,1])
plt.title('f2')
#plt.xlabel('time')
#plt.ylabel('mV')

plt.subplot(2, 3, 4)
plt.plot(f3.iloc[:,0], f3.iloc[:,1])
plt.title('f3')
#plt.xlabel('time')
#plt.ylabel('mV')

plt.subplot(2, 3, 5)
plt.plot(f4.iloc[:,0], f4.iloc[:,1])
plt.title('f4')
#plt.xlabel('time')
#plt.ylabel('mV')

plt.subplot(2, 3, 6)
plt.plot(f5.iloc[:,0], f5.iloc[:,1])
plt.title('f5')
#plt.xlabel('time')
#plt.ylabel('mV')
plt.show()

#plot_all(f5.iloc[:,0], ecg, f1.iloc[:,1], f2.iloc[:,1], f3.iloc[:,1], f4.iloc[:,1], f5.iloc[:,1])
