'''
Finds peaks using scipy.signal.find_peaks
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
import utils as ut
import os.path as path

def test(x, height_fract = 0, title='Peak Detection Using Scipy.signal'):
    #x = electrocardiogram()[2000:4000]
    max = np.max(x)
    height = height_fract * max
    peaks, _ = find_peaks(x, height=height)
    #print(f"data_points={len(x)} peaks={len(peaks)}")
    peak_vals = [x[i] for i in peaks]
    x_vals = range(0, len(x))
    #print(peaks)
    #print(peak_vals)
    plt.figure(figsize=(10,6))
    #plt.subplots_adjust(top=0.8)
    plt.plot(x)
    plt.plot(peaks, peak_vals, "ro")
    #plt.plot(np.zeros_like(x), "--", color="gray")
    plt.title(title + f"\n{len(peaks)} peaks: height_fract={height_fract} height={height}")
    plt.show()

def run():
    print(path.basename(path.normpath(__file__)))
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
    print(filename)
    ecg, headers = ut.read_ecg_file(filename)

    test(ecg, height_fract=.4, title='Peak Detection Using Scipy.signal\n' + filename)

if __name__ == "__main__":
    run()
