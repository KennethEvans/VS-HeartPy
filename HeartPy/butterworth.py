import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter
import os.path as path
import utils as ut
import filter as flt

def plot_frequency_response(sample_rate, a, b, title='Digital Filter Frequency Response'):
    w, h = signal.freqz(b=b, a=a)
    x = w * sample_rate * 1.0 / (2 * np.pi)
    # Avoid divide by zero warning
    np.seterr(divide = 'ignore') 
    y = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    plt.figure(figsize=(10,6))
    plt.semilogx(x, y)
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    plt.title(f"{title}" + "\n{a=}\n{b=}")
    plt.grid(which='both', linestyle='-', color='grey')
    plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000], ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.show()

def plot_frequency_response2(sample_rate1, a1, b1, sample_rate2, a2, b2, title='Digital Filter Frequency Response'):
    w, h = signal.freqz(b=b1, a=a1)
    x1 = w * sample_rate1 * 1.0 / (2 * np.pi)
    np.seterr(divide = 'ignore') 
    y1 = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    w, h = signal.freqz(b=b2, a=a2)
    x2 = w * sample_rate2 * 1.0 / (2 * np.pi)
    np.seterr(divide = 'ignore') 
    y2 = 20 * np.log10(abs(h))
    np.seterr(divide = 'warn') 
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.8)
    plt.semilogx(x1, y1, label=f"{sample_rate1} Hz")
    plt.semilogx(x2, y2, label=f"{sample_rate2} Hz")
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    plt.title(f"{title}\n{a1=}\n{b1=}\n{a2=}\n{b2=}")
    plt.grid(which='both', linestyle='-', color='grey')
    plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000], ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.legend(loc=4, framealpha=0.6)
    plt.show()

def example():
    sample_rate = 44100
    a = 1
    b = [1, -1]
    plot_frequency_response(sample_rate, a, b)

def butter_bandpass(fs, lowcut, highcut, order=5):
    '''Calculates a and b coefficients for a Butterworth filter with these parameters'''
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter_data(data, fs, lowcut, highcut, order=5):
    '''Calculates the filtered values for the input data.
    Uses scipy.signal.lfilter.
    '''
    b, a = butter_bandpass(fs, lowcut, highcut, order=order)
    y = lfilter(b, a, data)
    return y

def filter_data(b, a, data):
    '''Calculates the filtered values for the input data for a filter
    with the given a and b coefficients.
    Uses scipy.signal.lfilter.
    '''
    b, a = butter_bandpass(fs, lowcut, highcut, order=order)
    y = lfilter(b, a, data)
    return y

def plot_butterworth_frequency_response(fs, low_cutoff, high_cutoff, orders, print_coefficients=False):
    '''Plots the frequency response for a Butterworth filter.
    Uses butter_bandpass.
    
        Parameters
        ----------
        fs : float
            The sampling rate.
    
        low_cutoff : float
            The low cutoff.
    
        high_cutoff : float
            The high cutoff.
    
        orders : list of int
            The orders to calculate.
    
        print_coefficients : boolean
            Whether to print the coefficients in Pythion format or not.
            default: False
    
        Returns
        ------
        No return value.
    '''

    # Plot the frequency response for a few different orders.
    #plt.figure(1)
    #plt.clf()
    plt.figure(figsize=(10,6))
    for order in orders:
        b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
        print(f"\n{order=} {len(a)=} {len(b)=}")
        if print_coefficients:
            # Generate Python code
            out = f"A_BUTTERFIELD{order}=["
            for i in range(len(a)):
                out += f"{a[i]}, "
            out += "]"
            print(out)
            out = f"B_BUTTERFIELD{order}=["
            for i in range(len(b)):
                out += f"{b[i]}, "
            out += "]"
            print(out)

        w, h = signal.freqz(b, a, worN=2000)
        x = w * fs * 1.0 / (2 * np.pi)
        y = 20 * np.log10(abs(h))
        plt.semilogx(x, y, label="order = %d" % order)
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    plt.title(f"Butterworth ({fs=} Hz l{low_cutoff=} {high_cutoff=})")
    plt.grid(which='both', linestyle='-', color='grey')
    plt.legend(loc='best')
    plt.xticks([20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000],
              ["20", "50", "100", "200", "500", "1K", "2K", "5K", "10K", "20K"])
    plt.show()

def get_a0_for_frequency_response(a, b):
    '''Calculates the value of a[0] that puts the maximum of 20log10(abs(h))
    at 0 for a FIR filter with coefficients a, b.
    
        Parameters
        ----------
        a : list of float
            The a coefficients.
    
        b : list of float
            The b coefficients.
    
        Returns
        ------
        a0new : float
            The calculated value of a[0].
    '''
    w, h = signal.freqz(b=b, a=a)
    np.seterr(divide = 'ignore') 
    y = 20 * np.log10(abs(h))
    max_y = np.max(y)
    np.seterr(divide = 'warn')
    a0new = a[0] * 10 ** (max_y / 20)
    return a0new

def test_butterworth(fs, low_cutoff, high_cutoff, order, filename):
    # Filter ecg file
    ecg, headers = ut.read_ecg_file(filename)
    x= range(0, len(ecg))
    b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
    filtered = butter_bandpass_filter_data(ecg, fs, low_cutoff, high_cutoff, order=order)
    plt.figure(figsize=(10,6))
    plt.plot(x, ecg, linewidth=1, label='ECG Signal')
    plt.plot(x, filtered, color='#89CFF0', linewidth=1, label=f"Filtered Signal")
    #plt.xlabel('time (seconds)')
    #plt.hlines([-a, a], 0, T, linestyles='--')
    plt.title(f"Butterworth ({fs=} {low_cutoff=} {high_cutoff=} {order=} {len(a)=} {len(b)=}\n{filename}")
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='lower center')
    plt.show()

def test_butterworth_ke(fs, low_cutoff, high_cutoff, order, filename):
    # Filter ecg file
    ecg, headers = ut.read_ecg_file(filename)
    x= range(0, len(ecg))
    b, a = butter_bandpass(fs, low_cutoff, high_cutoff, order=order)
    filtered = flt.filter_data(a, b, ecg)
    #print(f"{len(ecg)} {ecg[-5:]=}")
    #print(f"{len(filtered)} {filtered[-5:]=}")
    plt.figure(figsize=(10,6))
    plt.plot(x, ecg, linewidth=1, label='ECG Signal')
    plt.plot(x, filtered, color='#89CFF0', linewidth=1, label=f"Filtered Signal")
    #plt.xlabel('time (seconds)')
    #plt.hlines([-a, a], 0, T, linestyles='--')
    plt.title(f"Butterworth Using Filter ({fs=} {low_cutoff=} {high_cutoff=} {order=} {len(a)=} {len(b)=}\n{filename}")
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='lower center')
    plt.show()

def butterworth0():
    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 5000.0
    lowcut = 500.0
    highcut = 1250.0

    # Plot the frequency response for a few different orders.
    plt.figure(1)
    plt.clf()
    for order in [3, 6, 9]:
        b, a = butter_bandpass(fs, lowcut, highcut, order=order)
        w, h = signal.freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')

    # Filter a noisy signal.
    T = 0.05
    nsamples = int(T * fs)
    t = np.linspace(0, T, nsamples, endpoint=False)
    a = 0.02
    f0 = 600.0
    x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
    x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
    x += a * np.cos(2 * np.pi * f0 * t + .11)
    x += 0.03 * np.cos(2 * np.pi * 2000 * t)
    plt.figure(2)
    plt.clf()
    plt.plot(t, x, label='Noisy signal')

    y = butter_bandpass_filter_data(fs, x, lowcut, highcut, order=6)
    plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
    plt.xlabel('time (seconds)')
    plt.hlines([-a, a], 0, T, linestyles='--')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='upper left')

    plt.show()


def ptLow250():
    sample_rate = 250
    a = [1, -2, 1]
    b = [1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]
    plot_frequency_response(sample_rate, a, b)

def ptLow130():
    sample_rate = 130
    a = [1, -2, 1]
    b = [1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]

    plot_frequency_response(sample_rate, a, b)

def ptlow250_130():
    sample_rate1 = 250
    a1 = [1, -2, 1]
    b1 = [1, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]
    sample_rate2 = 130
    a2 = a1
    b2 = b1
    plot_frequency_response2(sample_rate1, a1, b1, sample_rate2, a2, b2)

def pan_tompkins_derivative():
    sample_rate1 = 200
    # Empirically choose a[0] to make the gain 0 at the peak
    a1 = [5.4716104436]
    b1 = [2, 1, 0, -1, -2]
    sample_rate2 = 130
    a2 = a1
    b2 = b1
    if True:
        # Calculate what a[0] should be to make y = 0 at the top
        a0new = get_a0_for_frequency_response(a1, b1)
        print(f"{a0new=}")
    plot_frequency_response2(sample_rate1, a1, b1, sample_rate2, a2, b2,
        title='Pan Tompkins Derivative Frequency Response')

def main():
    print(path.basename(path.normpath(__file__)))
    # Working on computer HR=55
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
    # Walking Falison HR=114
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
    # Walking Blueberry Lake HR=120
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
    # Feb 4 Example New Low Heartrate HR=63
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'

    if True:
        #example()
        #ptlow250_130()
        #pan_tompkins_derivative()
        #plot_butterworth_frequency_response(130, 5, 20, orders=[2,3,4,5,6,9], print_coefficients=True)
        #test_butterworth(130, 5, 15, order=3, filename=filename)
        test_butterworth(130, 5, 20, order=3, filename=filename)
        test_butterworth_ke(130, 5, 20, order=3, filename=filename)
    else:
        ptLow250()
        ptLow130()

if __name__ == "__main__":
    main()
