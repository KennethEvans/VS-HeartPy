'''
Commonly used utility functions for ECG processaing.
'''

import os
from csv import reader
from datetime import datetime, timedelta

def timestamp(format = '%Y-%m-%d %H:%M:%S.%f'):
    '''Generates a timestamp for the current time with the given format
    '''
    now = datetime.now()
    return  now.strftime(format)

def local_utc_offset():
    '''Gets the UTC offset of the local time'''
    millis = 1288483950000
    ts = millis * 1e-3
    # local time == (utc time + utc offset)
    utc_offset = datetime.fromtimestamp(ts) - datetime.utcfromtimestamp(ts)
    return utc_offset

def prompt_for_files(title='Choose Files', multiple=False, type='csv'):
    '''Brings up a system dialog to prompt for files of the given type and
    an 'All files (*.*)" files option. You can pick the dialog title, and
    whether multiple files are allowed.
    '''
    # Prompt for the file name
    import tkinter as tk
    from tkinter import filedialog
    import ctypes
    ctypes.windll.shcore.SetProcessDpiAwareness(1)
    root = tk.Tk()
    root.withdraw() 
    file_names = []
    # Prompt to open a GPX file
    filetypes = ((f'{type.upper()} files', f'*.{type.lower()}'),
        ('All files', '*.*'))
    if multiple:
        file_names_ret = filedialog.askopenfilenames(title=title,
            multiple= True, filetypes=filetypes)
        for file_name in file_names_ret:
            file_names.append(file_name)
    else:
        file_names_ret = filedialog.askopenfilename(title=title,
            filetypes=filetypes)
        if file_names_ret and len(file_names_ret) > 0:
            file_names.append(file_names_ret)
    root.destroy()
    return file_names

def read_ecg_file(fileName):
    '''Reads an ECG file with one column of ECG values and potentially one
column of peak values, treating lines that are not a float as part of the header.

    Parameters
    ----------
    filename : str
        Name of the file to read

    Returns
    -------
    ecg : list
        The ecg values.
    peaks : list
        The peak values.
    headers : list
        The headers
    '''
    # This gives an array of lines w/o CR or LF
    with open(fileName) as fd:
        lines = fd.read().splitlines()
    ecg = []
    peaks = []
    headers = []
    # Just use lines that are floats
    for line in lines:
        if len(line) == 0:
            continue
        try:
            line.replace('\n', '')
            tokens = line.split(",");
            ecgval = float(tokens[0])
            ecg.append(ecgval)
            if len(tokens) > 1:
                if int(tokens[1]) == 0:
                    peaks.append(False)
                else:
                   peaks.append(True)
        except:
            headers.append(line)
    return ecg, peaks, headers

def read_event_timer_file(fileName):
    '''Reads an Event Timer CSV file.

    Parameters
    ----------
    filename : str
        Name of the file to read

    Returns
    -------
    time : list of datetiem.datetime
        The time values.
    note : list of string
        The note values.
    headers : list of string
        The headers
    '''
    # This gives an array of lines w/o CR or LF
    with open(fileName) as fd:
        lines = fd.read().splitlines()
        times = []
        notes = []
        headers = []
        try:
            # Just use lines that are floats
            data_found = False
            for line in lines:
                if (not data_found):
                    if line.startswith('time'):
                        # This is the column names line
                        data_found = True
                        continue
                if (not data_found):
                    # Header
                    headers.append(line)
                else:
                    # Data
                    tokens = line.split("\t");
                    if len(tokens) != 3:
                        continue
                    note = tokens[1]
                    timestamp = float(tokens[2])/1000
                    timeval = datetime.fromtimestamp(timestamp)
                    times.append(timeval)
                    notes.append(note)
        except Exception as e:
            print(e)
    return times, notes, headers

def read_emay_spo2_file(fileName):
    '''Reads an EMAY CSV file with columns for date, time, Spo2, and PR.

    Parameters
    ----------
    filename : str
        Name of the file to read

    Returns
    -------
    time : list of datetime.datetime
        The time values.
    hr : list of int
        The hr values.
    spo2 : list of int
        The SpO2 values.
    headers : list
        The headers
    '''
    # This gives an array of lines w/o CR or LF
    with open(fileName) as fd:
        lines = fd.read().splitlines()
    time = []
    hr = []
    spo2 = []
    # Just use lines that are floats
    nlines = 0;
    for line in lines:
        # Skip the header
        nlines = nlines + 1
        if nlines == 1:
            continue
        try:
            line.replace('\n', '')
            tokens = line.split(",");
            if len(tokens) != 4:
                continue
            # Windows does not handle single-digit month and day
            mdy = tokens[0].split('/');
            if len(mdy[0]) == 1:
                month = f'0{mdy[0]}'
            else:
                month = f'{mdy[0]}'
            if len(mdy[1]) == 1:
                day = f'0{mdy[1]}'
            else:
                day = f'{mdy[1]}'
            year = f'{mdy[2]}';
            combinedtime = f'{month}/{day}/{year} {tokens[1]}'
            timeval = datetime.strptime(combinedtime, '%m/%d/%Y %I:%M:%S %p')
            time.append(timeval)
            if len(tokens[2]) == 0:
                spo2.append(None)
            else:
                spo2.append(int(tokens[2]))
            if len(tokens[3]) == 0:
                hr.append(None)
            else:
                hr.append(int(tokens[3]))
        except Exception as e:
            print(e)
            if nlines < 2:
                continue
            else:
                break
    return time, spo2, hr

def get_peak_values(ecg, is_peak):
    '''Finds the ECG values that are peaks.

    Parameters
    ----------
    ecg : list
        List of ECG values.
    isPeak : list
        List of booleans indicating if it is a peak or not.

    Returns
    -------
    peaks : list
        The indices corresponding to ECG peaks.
    '''
    peaks = [];
    for i in range(len(ecg)):
        if is_peak[i]:
            peaks.append(i)
    return peaks

def find_header_item(headers, item):
    '''Finds the value of the line starting with <item>= in the headers.
    Returns None if not found. Will only check new header format.
    '''
    use_new = False;
    if len(headers) > 0:
        if headers[0].startswith("application="):
            use_new = True;
    if not use_new:
       return None
    for line in headers:
            # Of the form note=<note>
            if line.startswith(item + '='):
                return line[len(item) + 1:]
    return None

def find_header_description(headers):
    use_new = False;
    if len(headers) > 0:
        if headers[0].startswith("application="):
            use_new = True;
    for item in headers:
        if use_new:
            # Of the form note=<note>
            if item.startswith("note="):
                return item[5:]
        else:
            # Line that does not have : or start with HR or Polar
            if not (':' in item) and not item.startswith("HR") and not item.startswith('Polar'):
                return item
    return None

def read_session_file(file_name):
    '''Reads a Session file with 3 columns, time, hr, rr.

    Parameters
    ----------
    flename : str
        Name of the file to read

    Returns
    -------
    time : list
        The time values.
    hr : list
        The HR values.
    rr : list
        The RR values.
    '''
    time = []
    hr = []
    rr = []
    with open(file_name, 'r') as read_obj:
        # pass the file object to reader() to get the reader object
        csv_reader = reader(read_obj)
        # Iterate over each row in the csv using reader object
        for row in csv_reader:
            # row variable is a list that represents a row in csv
            len_row = len(row)
            if not row or len_row == 0:
                continue
            dt = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S.%f')
            time.append(dt)
            if len_row < 2:
                hr.append(None)
                rr.append(None)
                continue
            hr_val = int(row[1])
            hr.append(hr_val)
            if len_row < 3:
                rr.append(None)
                continue
            rr_vals = row[2].split()
            if len(rr_vals) == 0:
                rr.append(None)
            else:
                # TODO make this right
                if 'Invalid' in rr_vals:
                    rr.append(None)
                else:
                    first = True
                    for rr_val in rr_vals:
                        rr1 = int(rr_val)
                        if first:
                            rr.append(rr1)
                            first = False
                        else:
                            # Add rr1 as milliseconds (Actually it is 1/1024 sec)
                            time.append(dt + timedelta(milliseconds=rr1))
                            hr.append(hr_val)
                            rr.append(rr1)

    return time, hr, rr

def read_halo_rise_file(file_name):
    '''Reads a Halo Rise file with 3 columns, time, hr, rr.

    Parameters
    ----------
    flename : str
        Name of the file to read

    Returns
    -------
    halo_sessions : list of list of time, height
        values for each segment each session
    '''
    if not os.path.isfile(file_name):
        print(f'Does not exist: {file_name}')
        return None
    with open(file_name, 'r') as read_obj:
        # pass the file object to reader() to get the reader object
        csv_reader = reader(read_obj)
        #nrows = len(read_obj.readlines())
        #print(f'{csv_reader=} {nrows=}')
        # Iterate over each row in the csv using reader object
        rownum = -1
        analysis_col = None
        halo_sessions = []
        for row in csv_reader:
            # row variable is a list that represents a row in csv
            len_row = len(row)
            rownum += 1
            if not row or len_row == 0:
                continue
            if rownum == 0:
                colnum = -1
                for col in row:
                    colnum = colnum + 1
                    if "Sleep Analysis" in col:
                        analysis_col = colnum
                        if not analysis_col:
                            print('f{row=}')
                            print('read_halo_rise_file: Sleep Analysis column not found')
                            return None
                        #else:
                        #    print(f'{analysis_col=}')
                continue
            tokens = row[analysis_col].split(';')
            halo_session = []
            for item in tokens:
                #print(f'{item=}')
                items = item.split('|')
                if len(items) != 3:
                    continue
                if items[0] == 'WAKE':
                    height = 4
                elif items[0] == 'REM':
                    height = 3
                elif items[0] == 'LIGHT':
                    height = 2
                elif items[0] == 'DEEP':
                    height = 1
                else:
                    # Shouldn't happen
                    height = 0
                start = datetime.strptime(items[1], '%Y-%m-%dT%H:%M:%S.%fz')
                end = datetime.strptime(items[2], '%Y-%m-%dT%H:%M:%S.%fz')
                #print(f'{height=} {start=} {end=}')
                halo_session.append([start, height])
                halo_session.append([end, height])
            halo_sessions.append(halo_session)
            #print(f'{rownum} {len(tokens)=}')
            #print(f'{rownum=} {row[analysis_col]=}\n')
            #dt = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S.%f')
        return halo_sessions

def write_ecg_file(ecg, peaks, header, filename = r'data\ecg_test_data.csv',
                   simulation_str='Heart.py'):
    '''Writes a CSV file with 2 columns, ecg, is_peak. There is a header at the top.

    Parameters
    ----------
    ecg : list
        The ecg values.
    peaks : list
        The peak values.
    headers : list
        The headers
    filename : str
    '''
    with open(filename, 'w') as f:
        if simulation_str:
            f.write(f'simulation={simulation_str}\n')
        for line in header:
            f.write(line + '\n')
        for i in range(len(ecg)):
            if not peaks is None:
                try:
                    if peaks.index(i):
                        is_peak = 1
                except:
                    is_peak = 0
            f.write(f'{ecg[i]:.3f},{is_peak}\n')

def format_duration(seconds: int):
    if seconds is not None:
        seconds = int(seconds)
        d = seconds // (3600 * 24)
        h = seconds // 3600 % 24
        m = seconds % 3600 // 60
        s = seconds % 3600 % 60
        if d > 0:
            return '{:02d}D {:02d}H {:02d}m {:02d}s'.format(d, h, m, s)
        elif h > 0:
            return '{:02d}H {:02d}m {:02d}s'.format(h, m, s)
        elif m > 0:
            return '{:02d}m {:02d}s'.format(m, s)
        elif s > 0:
            return '{:02d}s'.format(s)
    return '-'

def pip_list(print_stdout=True):
    '''Runs pip list and optionally prints the output'''
    import subprocess, sys
    print('pip_list:')
    args = [sys.executable, '-m', 'pip', 'list']
    p = subprocess.run(args, check=True, capture_output=True)
    if print_stdout:
        print(p.stdout.decode())
    return p.stdout.decode()

def pip_show(package, print_stdout=True ):
    '''Runs pip show for the package and optionally prints the output'''
    import subprocess, sys
    args = [sys.executable, "-m", "pip", "show", package]
    p = subprocess.run(args, check=True, capture_output=True)
    if print_stdout:
        print(p.stdout.decode())
    return p.stdout.decode()

def print_module_filename(file):
    ''' Prints the name of the module. Use __file__ for the desired module.
        Example: print_module_name(__file__)
    '''
    print(os.path.basename(os.path.normpath(file)))