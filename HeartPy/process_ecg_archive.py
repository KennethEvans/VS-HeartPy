'''
This module contains functions that were formerly in process_ecg. But which
are no longer so much of interest. These include the animations, which were
never fully developed.

run_animate_ion() animates this using plt.ion().
run_animate() animates this using FuncAnimation.
    It is necessary to use global variables to do this.
    It is faster than using plt.ion with blit=True, but the axes labels are not redrawn.
    It seems to be slightly faster than using plt.ion with blit=False.
'''

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.fft import fft, fftfreq
import os
import time
from collections.abc import Iterable
import utils as ut
import filter as flt

import process_ecg as proc

# Sampling rate. These algorithme are based on this particular sampling rate.
FS = 130.0
# Width of plot to show
PLOT_WIDTH_SEC = 3
PLOT_WIDTH_SAMPLES = 3 * FS
PLOT_WIDTH_MARGIN = .1
# Number of mV to show
PLOT_HEIGHT_MV = 2
# Window size.  Must be large enough for maximum number of coefficients
WINSIZE = 20
# Set whether to be fast (blit=True) with no x ticks or blit=False and slow
FAST = True

def test():
    global fig, ax, redDot
    TWOPI = 2*np.pi

    fig, ax = plt.subplots()

    t = np.arange(0.0, TWOPI, 0.001)
    s = np.sin(t)
    l = plt.plot(t, s)

    ax = plt.axis([0,TWOPI,-1,1])

    redDot, = plt.plot([0], [np.sin(0)], 'ro')

    # create animation using the animate() function
    myAnimation = animation.FuncAnimation(fig, animate1, frames=np.arange(0.0, TWOPI, 0.1), \
                                          interval=10, blit=True, repeat=False)

    plt.show()

def animate1(i):
    print(f"animate1 {i}")
    redDot.set_data(i, np.sin(i))
    return redDot,

def init():
    global line1
    #print('init start')
    ax.set_xlim(0, PLOT_WIDTH_SEC)
    ax.set_ylim(-PLOT_HEIGHT_MV, PLOT_HEIGHT_MV)
    #if isinstance(line1 , Iterable):
    #    print(f"line1 is iterable") 
    #else:
    #   print(f"line1 is not iterable")
    #temp = (line1)
    #if isinstance(temp , Iterable):
    #    print(f"(line1) is iterable") 
    #else:
    #   print(f"(line1) is not iterable")
    #temp = (line1, ax)
    #if isinstance(temp , Iterable):
    #    print(f"(line1, ax) is iterable") 
    #else:
    #   print(f"(line1, ax) is not iterable")
    print(f"ax={hash(ax):#x}")
    #print(f"line1={hash(line1):#x}")
    #print('init end')
    # Note the comma
    return line1,
    #return (line1, ax)

def animate(i):
    #print(f"animate {i}")
    # Set up the window
    if i <= WINSIZE:
        window = ecg.copy();
    else:
        window = ecg[-WINSIZE:].copy()
    nwin = len(window)

    cur_ecg.append(ecg[i])
    cur_x.append(i / FS)
    line1.set_xdata(cur_x)
    line1.set_ydata(cur_ecg)
    if i > (PLOT_WIDTH_SEC - PLOT_WIDTH_MARGIN) * FS:
        ax.set_xlim(i / FS - PLOT_WIDTH_SEC + PLOT_WIDTH_MARGIN,
            i / FS + PLOT_WIDTH_MARGIN)
    #print(ax.get_xlim())
    #print(f"ax={hash(ax):#x}")
    #print(f"line1={hash(line1):#x}")

    # Note the comma
    return line1,
    #return (line1, ax)

def run_animate_ion(filename):
    '''Reads the file and loops over ECG values.
    
        Parameters
        ----------
        filename : string
            The file name to read.
    
        Returns
        ------
        None
    '''

    ecg, headers = ut.read_ecg_file(filename)
    description = ut.find_header_description(headers)
    necg = len(ecg)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)

    # Set up the plot
    plot_width = FS * 5
    cur_ecg = []
    cur_x = []
    plt.ion()
    fig, ax = plt.subplots(figsize=(8, 6))
    line1, = ax.plot(cur_x, cur_ecg)
    ax.set_xlim(0, PLOT_WIDTH_SEC)
    ax.set_ylim(-PLOT_HEIGHT_MV, PLOT_HEIGHT_MV)
    if description:
        title = f'Real Time Processing Using Pyplot.ion() {filename}\n{description}'
    else:
        title = f'Real Time Processing Using Pyplot.ion()\n{filename}'
    plt.title(title)
    plt.xlabel('time, sec')
    plt.ylabel('ECG, mV')
 
    # Loop
    #for i in range(0, necg):
    for i in range(0, 500):
        # Set up the window
        if i <= WINSIZE:
            window = ecg.copy();
        else:
            window = ecg[-WINSIZE:].copy()
        nwin = len(window)

        cur_ecg.append(ecg[i])
        cur_x.append(i / FS)
        line1.set_xdata(cur_x)
        line1.set_ydata(cur_ecg)
        if i > (PLOT_WIDTH_SEC - PLOT_WIDTH_MARGIN) * FS:
            ax.set_xlim(i / FS - PLOT_WIDTH_SEC + PLOT_WIDTH_MARGIN,
                i / FS + PLOT_WIDTH_MARGIN)
        fig.canvas.draw()
        if i < necg - 2:
            fig.canvas.flush_events()
        # Cannot make if faster
        #time.sleep(.000001)
    plt.ioff()
    # This is necessary to keep the plot up
    plt.show()

def run_animate():
    global ecg, necg, fig, ax, line1, cur_ecg, cur_x, window, nwin 

    ecg, headers = ut.read_ecg_file(filename)
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    necg = len(ecg)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)

    # Set up the plot
    cur_ecg = []
    cur_x = []
    fig, ax = plt.subplots(figsize=(8, 6))
    line1, = ax.plot(cur_x, cur_ecg)
    if description:
        title = f'Simulated Real Time Processing\n{filename}\n{description}'
    else:
        title = f'Simulated Real Time Processing\n{filename}'
    plt.title(f'ECG\n{filename}')
    plt.xlabel('time, sec')
    plt.ylabel('ECG, mV')
    if FAST:
        blit = True
        # Turn off tick labels
        ax.set_xticklabels([])
    else:
      blit = False

    #npoints = 500
    npoints = necg
    ani = animation.FuncAnimation(fig, animate, np.arange(npoints), interval=1,
        init_func=init, blit=blit, repeat=False)
    plt.show()
 
def run_real_time_1():
    '''
    Implementation using cascading filters with one cur_ecg and cur_filter.
    '''
    ecg, headers = ut.read_ecg_file(filename)
    necg = len(ecg)
    print(filename)
    print('\nHeader Information:')
    for header in headers:
        print(header)
    description = ut.find_header_description(headers)
    if description:
        title = f'Real Time Processing\n{filename}\n{description}'
    else:
        title = f'Real Time Processing\n{filename}'

    # Set up the plot
    bandpass = []
    deriv = []
    square= []
    avg=[]
    score=[]
    cur_ecg = []
    cur_filt = []
    cur_x = []
    
    keep = WINSIZE # Keep this many items

    # Loop over ECG values
    for i in range(0, necg):
       if len(cur_ecg) == keep:
            cur_ecg.pop(0)
            cur_filt.pop(0)
       cur_ecg.append(ecg[i])
       cur_filt.append(0) # Doesn't matter
       cur_x.append(i / FS)
       #print(f"{i} {len(cur_ecg)=}")

       # Butterworth
       new = flt.butterworth3(cur_ecg, cur_filt)
       cur_filt[-1] = new
       bandpass.append(new)

       # Derivative
       new = flt.derivative(cur_ecg, cur_filt)
       cur_filt[-1] = new
       deriv.append(new)

       # Square
       new = flt.square(cur_ecg, cur_filt)
       cur_filt[-1] = new
       square.append(new)

       # Moving average
       new = flt.moving_average1(cur_filt, WINSIZE)
       cur_filt[-1] = new
       avg.append(new)
       #print(f"{i} avg = {avg[0:5]}")

       ## Score
       threshold = .045
       new = flt.score(cur_filt[-1], threshold)
       score.append(new)

    ## Remove low frequency
    #avg = avg - np.mean(avg)

    # FFT
    if False:
        plot_fft(ecg, FS, filename = filename)
        plot_fft(avg, FS, title = 'FFT (Average)', filename = filename)

    #vals1 = bandpass
    #label1 = 'Bandpass'
    #vals1 = deriv
    #label1 = 'Derivative'
    #vals1 = square
    #label1 = 'Square'
    #vals1 = square
    #label1 = 'Score'
    vals1 = np.ndarray.tolist(np.array(avg) * 10.)
    label1 = 'Average x 10'
    vals2 = score
    label2 = 'Score'
    plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
        use_time_vals=True, xlim=[0, PLOT_WIDTH_SEC])
    plot_2_values(ecg, cur_x, vals1, vals2, label1=label1, label2=label2, title=title,
        use_time_vals=True)

    plot_all(ecg, bandpass, deriv, square, avg, score, title=title)

def main():
    global filename

    print(os.path.basename(os.path.normpath(__file__)))

    # Working on computer HR=55
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
    # Walking Falison HR=114
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
    # Walking Blueberry Lake HR=120
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
    # Feb 4 Example New Low Heartrate HR=63
    #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
    # Feb 6 Walking
    filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-06_13-52.csv'

    #test()
    #test2()

    if False:
        # Animate with ion
        # Needs to be in a try/except to avoid error when terminated manually
        try:
            run_animate_ion(filename)
        except Exception as ex:
            name = f"{type(ex)}"
            if 'TclError' in name:
                print('TclError: Process probably terminated manually')
            else:
                print(ex)
                traceback.print_exc()
    if True:
        # Animate with FuncAnimation
        run_animate()

if __name__ == "__main__":
    main()
