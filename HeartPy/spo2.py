''' This module handles EMAY Oximeter SpO2 files.
'''

import matplotlib.pyplot as plt
import os
import utils as ut

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

def plot_sp02_hr(time, spo2, hr, subtitle1=None, miny=90):
    interval = (time[-1] - time[0]).total_seconds()
    duration = ut.format_duration(interval)

    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(top=0.83)
    title_used = f'EMAY Oximeter Results (Duration={duration})'
    if subtitle1:
        title_used += f"\n{subtitle1}"
    plt.suptitle(title_used)

    ax1 = plt.subplot(2, 1, 1)
    ax1.set_ylim([90, 100])
    plt.plot(time, spo2, 'cornflowerblue', label='SpO2')

    plt.xlabel('time')
    plt.ylabel('SpO2, %')
    #plt.yticks(np.arange(90, 100, 2))

    plt.subplot(2, 1, 2)
    plt.plot(time, hr, color='red', label='HR')
    plt.ylabel('HR, bpm')
    plt.xlabel('time')

    plt.tight_layout()
    plt.show()

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
    print(os.path.basename(os.path.normpath(__file__)))

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

def run():
    #read_spo2(True)
    read_spo2_session(False)

if __name__ == "__main__":
    run()
