''' 
This module reads and plot ECG Session files.
'''

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
from datetime import datetime
import utils as ut
from csv import reader

def plot(time1, hr1, rr1, time2, hr2, rr2, subtitle1=None, subtitle2=None):
    fig = plt.figure(figsize=(10,6))

    minTime = min(time1[0], time2[0])
    maxTime = max(time1[-1], time2[-1])
    duration = (maxTime - minTime).total_seconds() / 60


    plt.plot(time1, hr1, 'r', label='HR1')
    plt.plot(time2, hr2, color='hotpink', label='HR2')
    plt.xlabel('time, sec')
    plt.xlabel('HR, bpm')
    plt.ylim(bottom=0)
    title_used = f'Hr/RR Comparison (Duration={duration:.1f} min)'
    if subtitle1:
        title_used += f"\n{subtitle1}"
    if subtitle2:
        title_used += f"\n{subtitle2}"
    plt.title(title_used)

    ax2 = plt.gca().twinx()
    ax2.plot(time1, rr1, color='b', label='RR1')
    ax2.plot(time2, rr2, color='cornflowerblue', label='RR2')
    ax2.set_ylabel('RR, ms')
    ax2.set_ylim(bottom=0)
    fig.legend(loc='lower right', framealpha=0.6, bbox_to_anchor=(1,0),
              bbox_transform=ax2.transAxes)
    plt.tight_layout()
    plt.show()

def main():
    print(os.path.basename(os.path.normpath(__file__)))

    # Set prompt to use default filename or prompt with a FileDialog
    prompt = True
    if prompt:
        nFiles = -1
        while not nFiles == 2 and not nFiles == 0:
            file_names = ut.prompt_for_files(title='Choose 2 session files (not ECG files)',
                multiple=True, type='csv')
            nFiles = len(file_names)
            if not nFiles == 2:
                print('You must choose exactly two files')
        if nFiles == 0:
            print('Canceled')
            return;
        filename1 = file_names[0]
        filename2 = file_names[1]
    else:
        # First file
        filename1 = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-DeviceHR-2022-01-07_14-28.csv'
        # Second file
        filename2 = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-QRSHR-2022-01-07_14-28.csv'

    time1, hr1, rr1 = ut.read_session_file(filename1)
    time2, hr2, rr2 = ut.read_session_file(filename2)

    plot(time1, hr1, rr1, time2, hr2, rr2,
        subtitle1=filename1, subtitle2=filename2) 

if __name__ == "__main__":
    main()

