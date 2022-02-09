'''
Commonly used utility functions for ECG processaing.
'''

import os
from csv import reader
from datetime import datetime

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
        if file_names_ret:
            file_names.append(file_names_ret)
    root.destroy()
    return file_names

def read_ecg_file(fileName):
    '''Reads an ECG file with one column of ECG values and potentially one
column of peak values, treating lines that are not a float as part of the header.

    Parameters
    ----------
    flename : str
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
                if tokens[1] == 0:
                    peaks.append(False)
                else:
                   peaks.append(True)
        except:
            headers.append(line)
    return ecg, peaks, headers

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
            dt = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S.%f')
            rr_vals = row[2].split()
            hr_val = int(row[1])
            if len(rr_vals) == 0:
                time.append(dt)
                hr.append(hr_val)
                rr.append(None)
            else:
                # TODO make this right
                for rr_val in rr_vals:
                    rr1 = int(rr_val)
                    time.append(dt)
                    hr.append(hr_val)
                    rr.append(rr1)
    return time, hr, rr
