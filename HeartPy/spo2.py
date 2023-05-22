''' This module handles EMAY Oximeter SpO2 files.
    Also Event Timer and Halo Rise files.
'''

#from this import d # Gives Zen of Python    
from pickle import NONE
from string import printable
import matplotlib.pyplot as plt
import os
import numpy as np
import os
import utils as ut
from datetime import datetime, timedelta

# This is the glocal Halo file that has all sessions
HALO_FILE = r'C:\Users\evans\Documents\Personal\Calc\Halo_Rise_Output_Sleep_Sessions.csv'

def raw(path):
    '''Converts a path with / to one with \\'''
    new_path = path.replace('/', '\\');
    return new_path

def plot_session(time1, hr1, time2, hr2, filename1=None, filename2=None):
    fig = plt.figure(figsize=(10,6))

    interval1 = (time1[-1] - time1[0]).total_seconds()
    interval2 = (time2[-1] - time2[0]).total_seconds()
    duration = ut.format_duration(interval1) + ', ' + ut.format_duration(interval2)

    if filename1 and os.path.exists(filename1):
        label1 = os.path.basename(filename1).split('.')[0]
    else:
        label1 = 'HR1'
    if filename2 and os.path.exists(filename2):
        label2 = os.path.basename(filename2).split('.')[0]
    else:
        label2 = 'HR2'

    plt.plot(time1, hr1, 'r', label=label1)
    plt.plot(time2, hr2, color='hotpink', label=label2)
    plt.xlabel('time, sec')
    plt.ylabel('HR, bpm')
    #plt.ylim(bottom=0)
    title_used = f'HR Comparison (Durations={duration} min)'
    if filename1:
        title_used += f"\n{filename1}"
    if filename2:
        title_used += f"\n{filename2}"
    plt.title(title_used)
    plt.legend(loc=4, framealpha=0.6)
    plt.tight_layout()
    plt.show()

def plot_sp02_hr(time, spo2, hr, time2=None, halo_session=None,
                notes=None, subtitle1=None, miny=90):
    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.83)
    interval = (time[-1] - time[0]).total_seconds()
    duration = ut.format_duration(interval)
    title_used = f'EMAY Oximeter Results (Duration={duration})'
    if subtitle1:
        title_used += f"\n{subtitle1}"
        if halo_session:
            title_used += f'\n{HALO_FILE} ({halo_session[0][0].strftime("%Y-%m-%d")})'
    plt.suptitle(title_used)
    ax1 = plt.subplot(2, 1, 1)
    ymin = 80
    ax1.set_ylim([ymin, 100])
    plt.plot(time, spo2, 'cornflowerblue', label='SpO2')
    if halo_session:
        halo_x = []
        halo_y = []
        for item in halo_session:
            halo_x.append(item[0])
            halo_y.append(ymin + 3 + 2*item[1])
        plt.plot(halo_x, halo_y, 'darkblue', label='Halo Rise')
    # Add annotation
    if time2  and notes:
        for i in range(len(time2)):
            # Scatter the values around 99. Higher index is later time.
            pos = 99 - i
            plt.text(time2[i], pos, '  ' + notes[i], verticalalignment="center")
            plt.plot(time2[i], pos, marker="o", color="black")

    plt.xlabel('time')
    plt.ylabel('SpO2, %')
    #plt.yticks(np.arange(90, 100, 2))

    plt.subplot(2, 1, 2)
    plt.plot(time, hr, color='red', label='HR')
    plt.ylabel('HR, bpm')
    plt.xlabel('time')
    if halo_session:
        halo_x = []
        halo_y = []
        for item in halo_session:
            halo_x.append(item[0])
            halo_y.append(30 + 7 * item[1])
        plt.plot(halo_x, halo_y, 'darkblue', label='Halo Rise')
    # Add annotation
    if time2 and notes:
        # Can't use sum(hr) and len(hr) if some values are None
        sum = 0
        nvals = 0
        for i in range(len(hr)):
            if hr[i]:
                sum = sum + hr[i]
                nvals = nvals +1
        if nvals > 0:
            avg = float(sum) / nvals
        else:
            avg =60
        # Scatter the values around avg. Higher index is later time.
        for i in range(len(time2)):
            pos = avg - -4 - 2 * i
            plt.text(time2[i], pos, '  ' + notes[i], verticalalignment="center")
            plt.plot(time2[i], pos, marker="o", color="black")

    plt.tight_layout()
    plt.show()

def plot_halo_intervals(start_times, disturbances, intervals,
                        title=f'Halo Rise Disturbances and Average Time Between Disturbances\n{HALO_FILE}'):
    start_date = start_times[-1].strftime("%Y-%m-%d")
    end_date = start_times[0].strftime("%Y-%m-%d")
    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.83)
    plt.suptitle(f'{title}\nStart Date={start_date} End Date={end_date}')
    ax1 = plt.subplot(2, 1, 1)
    ymin = 80
    plt.plot(start_times, disturbances, 'cornflowerblue', label='Disturbances')
    plt.xlabel('time')
    plt.ylabel('Disturbances')

    plt.subplot(2, 1, 2)
    plt.plot(start_times, intervals, color='red', label='Avg. Interval')
    plt.xlabel('intervals')
    plt.ylabel('Average Time Between, hr')
 
    plt.tight_layout()
    plt.show()

def plot_halo_rise(time, halo_session):
    utc_offset = ut.get_local_utc_offset()
    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.83)
    title_used = f'\n{HALO_FILE} ({halo_session[0][0].strftime("%Y-%m-%d")})'
    title_used += f'\nFrom the bottom: Deep, Light, REM, Disturbance'
    plt.suptitle(title_used)
    plt.xlabel('time')
    ax = plt.gca()
    ax.axes.yaxis.set_visible(False)
    ymin = 0
    halo_x = []
    halo_y = []
    for item in halo_session:
        halo_x.append(item[0] + utc_offset)
        halo_y.append(ymin + item[1] / 4)
    plt.plot(halo_x, halo_y, 'darkblue', label='Halo Rise')
    plt.tight_layout()
    plt.show()
    return

def read_event_timer(prompt = False, do_halo_session=False, adjust=False):
    print(os.path.basename(os.path.normpath(__file__)))

    # Set prompt to use default filename or prompt with a FileDialog
    if prompt:
        file_names = ut.prompt_for_files(title='Choose EMAY Oximeter file',
            multiple=False, type='csv')
        if file_names and len(file_names) > 0:
            filename = raw(file_names[0])
        else:
            print('Canceled')
            return
    else:
        #filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220908-001532.csv'
        #filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220910-140911.csv'
        filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20221130-231117.csv'

    # Get the event timer file
    if prompt:
        file_names = ut.prompt_for_files(title='Choose Event Timer file',
            multiple=False, type='csv')
        if file_names and len(file_names) > 0:
            filename2 = file_names[0]
        else:
            print('Canceled Event Timer file')
            time2 = None
            notes = None
            filename2 = None
    else:
        filename2 = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EventTimer-Sleep_Test_Nov_30.csv'

    time, spo2, hr = ut.read_emay_spo2_file(filename)
    if filename2:
        time2, notes, headers = ut.read_event_timer_file(filename2)
    #print(f'{len(time)=} {len(time2)=}')
    #print(f'{type(time[0])=} {type(time2[0])=}')
    #for i in range(len(time2)):
    #    print(time2[i])

    # get the Halo session
    halo_session = None
    if do_halo_session:
        halo_sessions = ut.read_halo_rise_file(HALO_FILE)
        halo_session = get_halo_session_for_date(time[0], halo_sessions, adjust)
        if False:
            # Set a specific start time
            datetime_str = '04/24/23 12:00:00'
            time_test = datetime.strptime(datetime_str, '%m/%d/%y %H:%M:%S')
            halo_session = get_halo_session_for_date(time_test, halo_sessions, adjust)

    plot_sp02_hr(time, spo2, hr, time2=time2, halo_session=halo_session,
                notes=notes, subtitle1=filename) 

def read_halo_rise(datetime_str, adjust=False):
    print(os.path.basename(os.path.normpath(__file__)))
    print(f'Processing Halo Rise for date {datetime_str}\n  {HALO_FILE}')

    halo_sessions = ut.read_halo_rise_file(HALO_FILE)
    time = datetime.strptime(datetime_str, '%m/%d/%y %H:%M:%S')
    halo_session = get_halo_session_for_date(time, halo_sessions, adjust)

    plot_halo_rise(time, halo_session=halo_session)

def read_spo2(prompt = False):
    print(os.path.basename(os.path.normpath(__file__)))

    # Set prompt to use default filename or prompt with a FileDialog
    if prompt:
        file_names = ut.prompt_for_files(title='Choose EMAY Oximeter file',
            multiple=True, type='csv')
        if file_names and len(file_names) > 0:
            filename = file_names[0]
        else:
            print('Canceled')
            return
    else:
        #filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220908-001532.csv'
        #filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220910-140911.csv'
        filename = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220911-095717.csv'

    time, spo2, hr = ut.read_emay_spo2_file(filename)
    #print(f'{len(time)=}')
    plot_sp02_hr(time, spo2, hr, subtitle1=filename) 

def read_spo2_session(prompt = False):
    ut.print_module_filename(__file__)

    # Set prompt to use default filename or prompt with a FileDialog
    if prompt:
            file_names = ut.prompt_for_files(title='Choose 1st EMAY Oximeter or Session file',
                multiple=True, type='csv')
            if file_names and len(file_names) > 0:
                filename1 = file_names[0]
            else:
                print('Canceled')
                return

            file_names = ut.prompt_for_files(title='Choose 2nd EMAY Oximeter or Session file',
                multiple=True, type='csv')
            if file_names and len(file_names) > 0:
                filename2 = file_names[0]
            else:
                print('Canceled')
                return
    else:
        filename1 = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220909-092752.csv'
        filename2 = r'C:\Scratch\ECG\Android\S22\Current\SpO2\Kenneth_Evans_2022-09-09_09-27-57_Walking_Proud_Lake.session.csv'

        #filename1 = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220910-140911.csv'
        #filename2 = r'C:\Scratch\ECG\Android\S22\Current\BLE Cardiac Monitor\BCM-2022-09-10-14-09-00-Combined.csv'

        #filename1 = r'C:\Scratch\ECG\Android\S22\Current\SpO2\EMAY SpO2-20220911-095717.csv'
        #filename2 = r'C:\Scratch\ECG\Android\S22\Current\BLE Cardiac Monitor\BCM-2022-09-11-09-57-10.csv'


    if 'EMAY' in filename1:
        time1, spo21, hr1 = ut.read_emay_spo2_file(filename1)
    else:
        time1, hr1, rr1 = ut.read_session_file(filename1)
    if 'EMAY' in filename2:
        time2, spo21, hr2 = ut.read_emay_spo2_file(filename2)
    else:
        time2, hr2, rr2 = ut.read_session_file(filename2)

    #print(f'{len(hr1)=} {len(time1)=} {len(hr2)=} {len(time2)=}')
    #lentime = len(time1)
    #lenhr = len(hr1)
    #minlen = min(lentime, lenhr)
    #for i in range(minlen):
    #    print(f'{i} {time1[i]=} {hr1[i]=}')
    #if lentime > minlen:
    #    print(f'{lentime-1} {time1[lentime-1]=} hr1[lentime-1]=Undefined')
    #if lenhr > minlen:
    #    print(f'{lenhr-1} time1[lenhr-1=Undefined {hr1[lenhr-1]=}')

    plot_session(time1, hr1, time2, hr2,
        filename1=filename1, filename2=filename2)

def get_halo_session_for_date(time, sessions, adjust=False):
    '''Searches the array of sessions for the one which has the start time
       closest to the given time'''
    delta_min = timedelta.max.total_seconds()
    session_min = None
    for halo_session in sessions:
        if halo_session and len(halo_session[0]) == 2:
            halo_start = halo_session[0][0]
            delta = (halo_start - time).total_seconds()
            if abs(delta) < abs(delta_min):
                delta_min  = delta
                session_min = halo_session
    # Adjust time values to start at time
    if adjust:
        for item in session_min:
            item[0] = item[0] - timedelta(seconds=delta_min)
    return session_min

def get_halo_intervals(plot=True, print=False):
    utc_offset = ut.get_local_utc_offset()
    start_times = []
    intervals = []
    disturbances = []
    halo_sessions = ut.read_halo_rise_file(HALO_FILE)
    n_sessions = len(halo_sessions);
    for n in range(len(halo_sessions)):
        halo_session = halo_sessions[n]
        n_items = len(halo_session)
        if n_items == 0:
            continue
        start_time = None
        end_time = None
        disturbances_session = 0
        total_duration = 0
        disturbance_duration = 0
        interval = 0
        for m in range(0, n_items, 2):
            start_item = halo_session[m]
            end_item = halo_session[m+1]
            start = start_item[0] + utc_offset
            if not start_time:
                start_time = start
            end = end_item[0] + utc_offset
            end_time = end_item[0]
            duration = (end - start).total_seconds()
            total_duration = total_duration + duration
            # Is a disturbance if it lasts 5 min = 300 sec
            # Start and end are not included
            height = start_item[1]
            if height == 4 and m > 0 and m < n_items - 3 and duration >= 300:
                 disturbances_session = disturbances_session + 1
                 disturbance_duration = disturbance_duration + duration
                 #print(f'     {disturbances_session}\t{start.strftime("%H:%M")}\t{end.strftime("%H:%M")}\t{duration / 60}')
            None
        if disturbances_session > 1:
            interval = (total_duration - disturbance_duration) / (disturbances_session - 1) / 3600
            start_times.append(start_time)
            disturbances.append(disturbances_session)
            intervals.append(interval)
        #print(f'{start_time}      {disturbances}\t{interval}\t{total_duration}  \t{disturbances_duration}')
    if print:
        for s in range(len(start_times)):
            print(f'{start_times[s]}      {disturbances[s]}\t{intervals[s]}')
    if plot:
        plot_halo_intervals(start_times, disturbances, intervals)
    return start_times, disturbances, intervals


def run():
    ## Reads multiple Emay session files
    #read_spo2(True)

    ## Reads two session files, Emay or HR
    #read_spo2_session(True)

    ## Reads just the Halo Rise file for a given specific start time
    #datetime_str = '05/05/23 12:00:00'
    #read_halo_rise(datetime_str)

    # Get disturbances and intervals
    get_halo_intervals()
    
    ## Reads one Emay and one Event Timer file (if not canceled)
    #read_event_timer(prompt=True, do_halo_session=True, adjust=False)

if __name__ == "__main__":
    run()
