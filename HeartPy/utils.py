'''
Commonly used utility functions for ECG processaing.
'''

import os

def read_ecg_file(fileName):
    '''Reads an ECG file with one column of ECG values, treating lines that
    are not a float as part of the header.

    Parameters
    ----------
    flename : str
        Name of the file to read

    Returns
    -------
    data : list
        The ecg values.
    headers : list
        The headers
    '''
    # This gives an array of lines w/o CR or LF
    with open(fileName) as fd:
        lines = fd.read().splitlines()
    data = []
    headers = []
    # Just use lines that are floats
    for line in lines:
        try:
            line.replace('\n', '')
            val = float(line)
            data.append(val)
        except:
            headers.append(line)
    return data, headers

def find_header_description(headers):
    for item in headers:
        if not (':' in item) and not item.startswith("HR") and not item.startswith('Polar'):
            return item
    return None