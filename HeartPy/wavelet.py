import matplotlib.pyplot as plt
import numpy as np
import pywt
import utils as ut

name = 'Wavelet QRS Detection'

# Sampling rate.  These algorithms are based on this particular sampling rate.
FS = 130.0

def plot_peaks(ecg, ecg_x, peak_indices=None, title='Wavelet Detected Peaks',
              filename=None, xlim = None):
    plt.figure(figsize=(10,6))
    #plt.subplots_adjust(top=0.8)
    plt.plot(ecg_x, ecg)
    if peak_indices:
        peak_vals = [ecg[i] for i in peak_indices]
        peak_x = [ecg_x[i] for i in peak_indices]
        plt.plot(peak_x, peak_vals, "ro")
        title_used = f"{title} ({len(peak_indices)} Peaks)"
    else:
        title_used = title;
    if filename:
        title_used += f"\n{filename}"
    plt.title(title_used)
    if xlim:
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.xlabel('time, sec')
    plt.tight_layout()
    plt.show()

def plot_modulus(sig_thresh, ecg, avg, std, title='Detected High Values',
              filename=None, xlim = None):
    sig_thresh_1 = [i if i > 0 else None for i in sig_thresh]
    sig_thresh_1 = sig_thresh
    std_1 = [std for i in range(len(sig_thresh))]
    avg_1 = [avg for i in range(len(sig_thresh))]

    plt.figure(figsize=(10,6))
    #plt.subplots_adjust(top=0.8)
    plt.plot(ecg, color='darkgray')
    plt.plot(sig_thresh_1)
    plt.plot(std_1, color='red')
    plt.plot(avg_1, color='gray')
    title_used = f"{title} ({avg=:0.2f} {std=:0.2f})\n{filename})"
    plt.title(title_used)
    if xlim:
        plt.gca().set_xlim(xlim[0], xlim[1])
    plt.tight_layout()
    plt.show()

def qrs_detection(filename, ecg, sample_rate=FS, max_bpm=200, use_padding=True,
                 use_level=2, use_ecg=False, plot=True):
    '''
    Based on
    https://medium.com/@andrewtan_36013/electrocardiograms-qrs-detection-using-wavelet-analysis-a1070505efee

    ecg:         The input ecg values  
    max_bpm:     Determines the window for finding maximum modulus
                      size is (60.0 / max_bpm) * sample_rate
    use_padding: Pad the ecg values to have length a power of 2
    use_level:   (1 - max allowable for length of [padded] ecg)
                     Use detail for this level for finding peaks
    use_ecg:     Use the [padded] ecg directly, do not use detail coefficients
    '''
    ## Test odd size array
    #ecg = ecg + [0]

    max_level = pywt.swt_max_level(len(ecg))
    print(f'{len(ecg)=} {max_level=} {use_level=}')

    # Pad the array
    if use_padding:
        ecg_p = pad_array_max(ecg)
    else:
        ecg_p = ecg
    max_level_p = pywt.swt_max_level(len(ecg_p))
    if use_level > max_level_p:
        print(f'Cannot continue: use_level > max_level_p {max_level_p=} {use_level=}')
        return;

    if use_ecg:
        detail_used = ecg_p
    else:
        ## Stationary Wavelet Transform
        #n_levels = 2 # 2 is the max for these data (length is a power of 2) (could pad)
        n_levels = pywt.swt_max_level(len(ecg_p))
        # This is the index of the level to use
        # (The coeffs start at a high level and decrease)
        use_index = n_levels - use_level
        print(f'{len(ecg_p)=} {max_level_p=} {use_level=} {use_index=}')
        coeffs = pywt.swt(ecg_p, wavelet = "haar", level=n_levels, start_level=0, axis=-1)
        detail_used = coeffs[use_index][1] ##2nd level detail coefficients -> No the lowest level, highest frequency

    # Plot the coefficients
    if False:
        coeffs = pywt.swt(ecg_p, wavelet = "haar", level=n_levels, start_level=0, axis=-1)
        for i in range(n_levels):
            dd = coeffs[i][1] 
            nRows = len(coeffs)
            nCols = len(coeffs[0])
            x_ecg = [i / FS for i in range(len(dd))]
            if plot:
                plot_peaks(dd, x_ecg, title=f'Level {n_levels-i} Coefficients', filename=filename)
    
    ## Threhold the detail coefficients
    avg = np.mean(detail_used)
    std = np.std(detail_used)
    ## Notice it uses abs(i) not just i
    if use_ecg and False:
        sig_thresh = [i if i > 2.0 * std else 0 for i in detail_used - avg]
    else:
        sig_thresh = [abs(i) if abs(i) > 2.0 * std else 0 for i in detail_used - avg]

    non_zero = np.count_nonzero(sig_thresh)
    print(f'{avg=:.3f} {std=:.3f} {non_zero=} {len(sig_thresh)=}')
    title = 'Detected High Values'
    if use_padding:
        title= f'{title} (Use Padding)'
    if use_ecg:
        title = f'{title} (Use ECG)'
    else:
        title = f'{title} (Use Level={use_level})'
    if plot:
        plot_modulus(sig_thresh, detail_used, avg=avg, std=std, filename=filename, title=title)
    
    ## Find the maximum modulus in each window
    window = int((60.0 / max_bpm) * sample_rate)
    ecg_p_len = len(ecg_p)
    n_windows = int(ecg_p_len / window)
    modulus, qrs = [], []
    
    ##Loop through windows and find max modulus
    for i in range(n_windows):
        start = i * window
        end = min([(i+1) * window, ecg_p_len])
        mx = max(sig_thresh[start:end])
        if mx > 0:
            modulus.append( (start + np.argmax(sig_thresh[start:end]), mx))
    
    ## Merge if within max bpm
    merge_width = int((0.2) * sample_rate)
    i=0
    while i < len(modulus) - 1:
        ann = modulus[i][0]
        if modulus[i+1][0] - modulus[i][0] < merge_width:
            if modulus[i+1][1] > modulus[i][1]: # Take larger modulus
                ann = modulus[i+1][0]
            i += 1
        qrs.append(ann)
        i += 1 

    ## Pin point exact qrs peak
    window_check = int(sample_rate / 6)
    #ecg_normed = np.absolute((ecg_p-np.mean(ecg_p))/(max(ecg_p)-min(ecg_p)))
    r_peaks = [0] * len(qrs)
    
    for i, loc in enumerate(qrs):
        start = max(0, loc - window_check)
        end = min(ecg_p_len, loc + window_check)
        wdw = np.absolute(ecg_p[start:end] - np.mean(ecg_p[start:end]))
        pk = np.argmax(wdw)
        r_peaks[i] = start + pk

    print(f'{window=} {merge_width=} {window_check=}')
    print(f'{len(modulus)=} {len(qrs)=} {len(r_peaks)=}')
        
    return ecg_p, r_peaks

def next_power_of_2(length, max_power=32):
    i = 0;
    while length >= 2**i:
        i= i+1
    return 2**i

def pad_array_max(array):
    length = next_power_of_2(len(array))
    first = int((length - len(array)) / 2);
    last = first
    delta = length - len(array) - first - last
    last = last + delta
    print(f'{len(array)=} {length=} {first=} {last=}')
    return pywt.pad(array, (first, last), 'zero')

'''
#def QRS_detection_orig(signal,sample_rate,max_bpm):
#    ## Stationary Wavelet Transform
#    coeffs = pywt.swt(signal, wavelet = "haar", level=2, start_level=0, axis=-1)
#    d2 = coeffs[1][1] ##2nd level detail coefficients
    
#    ## Threhold the detail coefficients
#    avg = np.mean(d2)
#    std = np.std(d2)
#    sig_thresh = [abs(i) if abs(i)>2.0*std else 0 for i in d2-avg]
    
#    ## Find the maximum modulus in each window
#    window = int((60.0/max_bpm)*sample_rate)
#    sig_len = len(signal)
#    n_windows = int(sig_len/window)
#    modulus,qrs = [],[]
    
#    ##Loop through windows and find max modulus
#    for i in range(n_windows):
#        start = i*window
#        end = min([(i+1)*window,sig_len])
#        mx = max(sig_thresh[start:end])
#        if mx>0:
#            modulus.append( (start + np.argmax(sig_thresh[start:end]),mx))
    
#    ## Merge if within max bpm
#    merge_width = int((0.2)*sample_rate)
#    i=0
#    while i < len(modulus)-1:
#        ann = modulus[i][0]
#        if modulus[i+1][0]-modulus[i][0] < merge_width:
#            if modulus[i+1][1]>modulus[i][1]: # Take larger modulus
#                ann = modulus[i+1][0]
#            i+=1
                
#        qrs.append(ann)
#        i+=1 
#    ## Pin point exact qrs peak
#    window_check = int(sample_rate/6)
#    #signal_normed = np.absolute((signal-np.mean(signal))/(max(signal)-min(signal)))
#    r_peaks = [0]*len(qrs)
    
#    for i,loc in enumerate(qrs):
#        start = max(0,loc-window_check)
#        end = min(sig_len,loc+window_check)
#        wdw = np.absolute(signal[start:end] - np.mean(signal[start:end]))
#        pk = np.argmax(wdw)
#        r_peaks[i] = start+pk
        
#    return r_peaks
'''

def run(filename, sample_rate=FS, max_bpm=200, use_padding=True, use_level=2,
       use_ecg=False, plot=True):
    ecg, is_peak, headers = ut.read_ecg_file(filename)
    ecg_p, peak_indices = qrs_detection(filename, ecg, sample_rate=sample_rate,
                                        plot=plot,
                                       max_bpm=max_bpm,
                                       use_padding=use_padding,
                                       use_level=use_level,
                                       use_ecg=use_ecg)
    #print(f'{len(ecg)=} {len(ecg_p)=} {len(peak_indices)=}')
    #print(peak_indices)
    x_ecg = [i / FS for i in range(len(ecg_p))]
    title = 'Wavelet Detected Peaks'
    if use_padding:
        title= f'{title} (Use Padding)'
    if use_ecg:
        title = f'{title} (Use ECG)'
    else:
        title = f'{title} (Use Level={use_level})'
    if plot:
        plot_peaks(ecg_p, x_ecg, peak_indices, filename=filename, title=title)
    return ecg, x_ecg, peak_indices

def main():
    global filename
    ut.print_module_filename(__file__)

     # Set prompt to use default filename or prompt with a FileDialog
    prompt = False
    if prompt:
        file_names = ut.prompt_for_files(title='Choose a Polar CVS file with ECG Data',
            multiple=False, type='csv')
        nFiles = len(file_names)
        if nFiles == 0:
            print('Canceled')
            return;
        filename = file_names[0]
    else:
        # 0 Working on computer HR=55
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-31_15-27.csv'
        # 1 Walking Falison HR=114
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-14_16-22.csv'
        # 2 Walking Blueberry Lake HR=120
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-10-19_15-30.csv'
        # 3 Feb 4 Example New Low Heartrate HR=63
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-04_11-08.csv'
        # 4 Feb 6 Walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2021-02-06_13-52.csv'
        # After walking. Device, internal HR differ
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_15-08.csv'
        # Before walking 2. Device HR changed from 94 to 54
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-51.csv'
        # Before walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-44.csv'
        # Sitting down. Little out of breath.
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-03_18-43.csv'
        # 2-20-2022 Before walking
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-20_13-44.csv'
        # 2-15-2022-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-15_12-58.csv'
        # 2-18-16 Has extra processed peak
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-02-18_16-06.csv'
        # 5-13-2022
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-05-13_12-42.csv'
        # 6-15-18 With new algorithm in KE.Net ECG
        filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-06-15_18-30.csv'
        # 12-17-2022 Walking Hillside
        #filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-17_16-22.csv'
        filename = r'C:\Scratch\ECG\Polar ECG\CSV\PolarECG-2022-12-17_16-16.csv'

    if True:
        #ut.pip_show('PyWavelets')
        #ut.pip_list()
        run(filename, sample_rate=FS, max_bpm=200, use_padding=True,
            use_level=2, use_ecg=True)

if __name__ == "__main__":
    main()

