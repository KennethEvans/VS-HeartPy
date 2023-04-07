''' This module handles EMAY Oximeter SpO2 files.
'''

import matplotlib.pyplot as plt
import os
import numpy as np
import os
import utils as ut
from datetime import timedelta

# This is the glocal Halo file that has all sessions
HALO_FILE = r'C:\Scratch\AAA\Amazon Data\Halo_Rise_Output_Sleep_Sessions.csv'

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

def plot_sp02_hr(time, spo2, hr, time2=None, halo_session=None, notes=None,
                subtitle1=None, miny=90):
    interval = (time[-1] - time[0]).total_seconds()
    duration = ut.format_duration(interval)

    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.83)
    title_used = f'EMAY Oximeter Results (Duration={duration})'
    if subtitle1:
        title_used += f"\n{subtitle1}"
        if halo_session:
            title_used += f'\n{HALO_FILE}'
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

    plot_sp02_hr(time, spo2, hr, time2=time2, halo_session=halo_session,
                notes=notes, subtitle1=filename) 

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
    for item in session_min:
        item[0] = item[0] - timedelta(seconds=delta_min)
    return session_min

def run():
    # Reads multiple Emay session files
    #read_spo2(True)

    # Reads two session files, Emay or HR
    #read_spo2_session(True)
    
    # Reads one Emay and one Event Timer file (if not canceled)
    read_event_timer(prompt=True, do_halo_session=True, adjust=False)

if __name__ == "__main__":
    run()
