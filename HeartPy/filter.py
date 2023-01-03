'''
General digital filter routines for
    y[n] = (sumb - suma) / a[0]
    suma = sum from 1 to Q of a[i] * y[n-i]
    sumb = sum from 0 to P of b[i] * x[n-i]
    where Q = len(a) and P = len(x)
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
import os.path as path
import math

A_BUTTERWORTH2 = [1.0, -2.657395687639593, 2.8635274589099367,
               -1.5335175290358518, 0.36156443847608805,]
B_BUTTERWORTH2 = [0.08560470548729655, 0.0, -0.1712094109745931,
                0.0, 0.08560470548729655,]

A_BUTTERWORTH3 = [1.0, -4.026234474291334, 7.118704187414651,
                -7.142612123715484, 4.314550872956459,
                -1.4837877480823038, 0.2259301306922936,]
B_BUTTERWORTH3 = [0.025966345753506013, 0.0, -0.07789903726051804,
                0.0, 0.07789903726051804, 0.0, -0.025966345753506013,]

A_BUTTERWORTH4 = [1.0, -5.3887532415518145, 13.209147933667055,
                -19.31830613213619, 18.48092634163037, 
                -11.843119984451862, 4.960480657624009,
               -1.2420934569719146, 0.1429645156472169,]
B_BUTTERWORTH4 = [0.007820300889836627, 0.0, -0.03128120355934651,
                0.0, 0.04692180533901976, 0.0, -0.03128120355934651, 
                0.0, 0.007820300889836627,]

A_BUTTERWORTH5 = [1.0, -6.748184881012428, 21.133788637805445,
                -40.54370801140382, 52.837492602564595, -48.89699858837488, 
                32.53617491879561, -15.368572932198212, 4.933321311760557,
                -0.9726695875227551, 0.08959307936615561,]
B_BUTTERWORTH5 = [0.0023484342471058297, 0.0, -0.011742171235529147,
                0.0, 0.023484342471058295, 0.0, -0.023484342471058295,
                0.0, 0.011742171235529147, 0.0, -0.0023484342471058297,]

A_BUTTERWORTH6 = [1.0, -8.105939566353365, 30.892312936899007,
                -73.30136239159113, 120.72105657212964, -145.43788068444553, 
                131.4348366537356, -89.7692461839638, 45.98717984400798,
                -17.234871226828773, 4.4871579281068295, -0.7292139058407499,
                0.05601484457087554,]
B_BUTTERWORTH6 = [0.0007042113377843282, 0.0, -0.00422526802670597, 
                0.0, 0.010563170066764924, 0.0, -0.014084226755686564,
                0.0, 0.010563170066764924, 0.0, -0.00422526802670597,
                0.0, 0.0007042113377843282,]

A_BUTTERWORTH9 = [1.0, -12.174492360450209, 71.17028661809752,
                -265.5948905445492, 709.3254684475859, -1440.08805643756,
                2303.6996352214533, -2969.2334712044417, 3125.9411413326643,
                -2708.114145562456, 1935.2100339233334, -1137.9827622847183,
                546.5999314454891, -211.5613050734211, 64.54153745080335,
                -14.978374948972618, 2.4908765009296667, -0.26499323309018463,
                0.013581013010116836,]
B_BUTTERWORTH9 = [1.8914555515522274e-05, 0.0, -0.00017023099963970047, 
                0.0, 0.0006809239985588019, 0.0, -0.001588822663303871, 
                0.0, 0.0023832339949558063, 0.0, -0.0023832339949558063, 
                0.0, 0.001588822663303871, 0.0, -0.0006809239985588019, 
                0.0, 0.00017023099963970047, 0.0, -1.8914555515522274e-05,]
A_DERIVATIVE = [1]
B_DERIVATIVE = [0.5, 0.25, -0.25, -0.5]

def filter_alt(a, b, x, y):
    '''Calculates a result for a generalized filter. Only includes terms
   for which the x and y arrays are long enough to provide values. (Effectively
   uses x or y = 0 for these values.)

    y[n] = 1 / a[0] * (suma - sumb)
    suma = sum from 1 to q of a[j] * y[n-j], q = len(a)
    sumb = sum from 0 to p of b[j] * x[n-j], p = len(b)
    Uses the values at the ends of x and y.
        
        Parameters
        ----------
        a : list of float
            The a filter coefficients.
    
        b : list of float
            The b filter coefficients.
    
        x : list of float same length as y
            The x values.
    
        y : list of float, same length as x
            The y values. May be None if a only has a[0].
    
        Returns
        ------
        a[0] : float
            The calculated a[0].
    '''
    na = lena = len(a)
    nb = lenb = len(b)
    lenx = len(x)
    if(lenx < lenb):
        nb = lenx
    suma = 0
    if(y and lena > 1):
        leny = len(y)
        if(leny < lena):
           na = leny
        for i in range(1, na):
            suma += a[i] * y[leny - i - 1]
    sumb = 0
    for i in range(nb):
        sumb += b[i] * x[lenx - i - 1]
    return (sumb - suma) / a[0]


def filter(a, b, x, y):
    '''Calculates a result for a generalized filter. Returns 0 if x and y are 
    not long enough to provide sufficient values for the sums over the
    coefficients.

    y[n] = 1 / a[0] * (suma - sumb)
    suma = sum from 1 to q of a[j] * y[n-j], q = len(a)
    sumb = sum from 0 to p of b[j] * x[n-j], p = len(b)
    Uses the values at the ends of x and y.
        
        Parameters
        ----------
        a : list of float
            The a filter coefficients.
    
        b : list of float
            The b filter coefficients.
    
        x : list of float same length as y
            The x values.
    
        y : list of float, same length as x
            The y values. May be None if a only has a[0].
    
        Returns
        ------
        a[0] : float
            The calculated a[0] or 0 if the arrays are not long enough.
    '''

    # TODO Consider handling lenx < lenb and leny < lena differently than exit
    lena = len(a)
    lenb = len(b)
    lenx = len(x)
    if(lenx < lenb):
        return 0
    suma = 0
    if(y):
        leny = len(y)
        if(leny < lena):
           return 0
        for i in range(1, lena):
            suma += a[i] * y[leny - i - 1]
    sumb = 0
    for i in range(lenb):
        sumb += b[i] * x[lenx - i - 1]
    return (sumb - suma) / a[0]

def filter_data(a, b, data):
    '''Calculates the filtered values for the input data for a filter
    with the given a and b coefficients. Intended for processing the
    whole ecg array not for real-time. There are two choices for handling
    values for the first items, where the length of the data is too short
    to use all the coefficients.
    (a) Initialize filtered with zeros.
    (b) Initialize filtered with data.
    We are using (a). Note that the results are actually invalid for these cases.
    
        Parameters
        ----------
        a : list of float
            The a filter coefficients.
    
        b : list of float
            The b filter coefficients.

        data : list of float
            The input data.
    
        Returns
        ------
        filtered : list of float
            The filtered data.
    '''    
    count = len(data)
    filtered = [0] * count
    #filtered = data.copy()
    #print(f"{len(data)} {data[:5]=}")
    #print(f"{n=} {count=}")
    #print(f"{data[:0]=}")
    #print(f"{data[:1]=}")

    for n in range(1, count):
    #for n in range(1, 5):
        x = data[:n + 1]
        y = filtered[:n + 1]
        filtered[n] = filter(a, b, x, y)
    #print(f"{len(y)} {y[:5]=}")
    #print(f"{len(data)} {data[:5]=}")
    #print(f"{len(filtered)} {filtered[:5]=}")
    return filtered

def butterworth(a, b, x, y):
    '''Generalized Butterworth with coefficients given. Just runs filter
    '''
    return filter(a, b, x, y)
    
def butterworth2(x, y):
    a = A_BUTTERWORTH2
    b = B_BUTTERWORTH2
    return filter(a, b, x, y)
    
def butterworth3(x, y):
    a = A_BUTTERWORTH3
    b = B_BUTTERWORTH3
    return filter(a, b, x, y)
    
def butterworth4(x, y):
    a = A_BUTTERWORTH4
    b = B_BUTTERWORTH4
    return filter(a, b, x, y)
    
def butterworth5(x, y):
    a = A_BUTTERWORTH5
    b = B_BUTTERWORTH5
    return filter(a, b, x, y)
    
def butterworth6(x, y):
    a = A_BUTTERWORTH6
    b = B_BUTTERWORTH6
    return filter(a, b, x, y)
    
def derivative(x, y):
    a = A_DERIVATIVE
    b = B_DERIVATIVE
    return filter(a, b, x, y)
    
def square(x, y):
    #print(f"square: type(x)={type(x)}")
    return x[-1] * x[-1]

def moving_average1(x, winsize):
    '''Calculates a moving average.
    Implemented as a filter.

    Parameters
    ----------
    x: array of x values
    winsize: The size of the moving window.
    '''
    n = winsize if len(x) > winsize else len(x)
    #b = (np.ones(n) / n).tolist()
    b = [1 / n] * n
    a = [1]
    val = filter(a, b, x[-n:], None)
    #print(f"moving_average1: len(x)={len(x)} len(x[-n:])={len(x[-n:])} n={n}
    #val={val} mean={np.mean(x[-n:])}")
    #print(f"b={b}")
    #print(f"x={x}")
    #print(f"x[-n:]={x[-n:]}")
    return val
    
def moving_average2(x, winsize):
    '''Calculates a moving average.
    Uses numpy.mean.

    Parameters
    ----------
    x: array of x values
    winsize: The size of the moving window.
    '''
    n = winsize if len(x) > winsize else len(x)
    val = np.mean(x[-n:])
    #print(f"moving_average2: len(x)={len(x)} len(x[-n:])={len(x[-n:])} n={n}
    #val={val} mean={np.mean(x[-n:])}")
    return val

def moving_average_uniform_filter(x, winsize):
    '''Calculates a moving average for an array.
    Uses scipy.ndimage.filters.uniform_filter1d. Not for real-time filtering.

    Parameters
    ----------
    x: array of x values
    winsize: The size of the moving window.
    '''
    n = winsize if len(x) > winsize else len(x)
    return uniform_filter1d(x, size=n)

class Moving_Average(object):
    '''Class to handle a moving average'''
    def __init__(self, window, items=[]):
        self.window = window
        if type(items) is list:
            if len(items) <= window:
                self.items = items
            else:
                # Use last window item
                self.ietms = items[-window]
        else:
            self.items = []
            self.items.append(item);

    def avg(self):
        return sum(self.items) / len(self.items)

    def add(self, val):
        self.items.append(val)
        if len(self.items) > self.window:
            self.items.pop(0)

def score(value, threshold):
    '''Scores a value to 0 if it is below a threshold or 1 otherwise.
    
        Parameters
        ----------
        value : float
            The value to check.

        threshold : float
            The threshold to use.
    
        Returns
        ------
        return : boolean
            0 or 1 depending if below or at or above threshold.
    '''      
    if value < threshold:
        return 0
    else:
        return 1

def test_moving_average():
    n = 100
    #test = np.arange(n)
    #print(f"test[0]={test[0]}")
    #print(f"test[99]={test[99]}")
    #print(f"test[-1]={test[-1]}")
    #print(f"test[-2]={test[-2]}")
    #print(f"test-[0:]={test[-0:]}")
    #print(f"test[-1:]={test[-1:]}")
    #print(f"test[-2:]={test[-2:]}")

    #test1 = [2., 8., 0., 4., 1., 9., 9., 0.]
    #print(f"test1={test1}")
    #print(f"uniform_filter1d={uniform_filter1d(test1, size=3)}")

    # Need to use dtype or will get integers and roundoff in the results for
    # moving_average1
    #pure = np.linspace(-1, 1, n, dtype=float)
    pure = np.arange(0, n, dtype=float)
    # Make it reproducible
    np.random.seed(10)
    noise = np.random.normal(0, 1, n) * 20
    signal = pure + noise
    #signal = pure
    filtered1 = signal.copy()
    filtered2 = signal.copy()
    filtered3 = signal.copy()
    filtered4 = signal.copy()
    winsize = 3
    keep = 10 # Keep this many items
    mov_avg = Moving_Average(winsize)
    for i in range(n):
        if(i < keep):
            filtered1[i] = moving_average1(signal[:i + 1], winsize)
            filtered2[i] = moving_average2(signal[:i + 1], winsize)
        else:
            filtered1[i] = moving_average1(signal[i - keep:i + 1], winsize)
            filtered2[i] = moving_average2(signal[i - keep:i + 1], winsize)
        mov_avg.add(signal[i])
        filtered4[i] = mov_avg.avg()
    filtered3 = moving_average_uniform_filter(signal, winsize)
    plt.figure(figsize=(10, 6))
    #plt.xlabel('Time (s)')
    #plt.plot(x_data, data, label='ECG Data');
    plt.plot(signal, label='signal')
    plt.plot(filtered1, label='moving_average1', linestyle='solid')
    plt.plot(filtered2, label='moving_average2', linestyle='dashed')
    plt.plot(filtered4, label='moving_average_from_class', linestyle= (0, (5, 10)))
    plt.plot(filtered3, label=f'moving_average_uniform_filter (no delay)', linestyle='solid')
    #plt.plot(pure - .5 * winsize, filtered1, label='Filtered 1 Adjusted')
    #plt.plot(filtered1 - filtered2, label='1 - 2')
    #plt.plot(filtered1 - filtered3, label='1 - uniform')
    plt.title('Moving Average Calculation Comparison')
    plt.legend(loc='upper left', framealpha=0.6)
    plt.show()

def test_derivative():
    n = 101
    keep = 10   # Could be as low as 4
    signal = [math.sin(i * 4 * math.pi / (n - 1)) for i in range(n)]
    cur_ecg = []
    cur_deriv=[]
    cur_square1=[]
    cur_square2=[]

    deriv = []
    square1 = []
    square2 = []
    for i in range(n):
        if len(cur_ecg) == keep:
            cur_ecg.pop(0)
        cur_ecg.append(signal[i])

        # Derivative
        input = cur_ecg
        if len(cur_deriv) == keep:
             cur_deriv.pop(0)
        cur_deriv.append(0) # Doesn't matter
        new = derivative(input, cur_deriv)
        cur_deriv[-1] = new
        cur_deriv.append(new)
        deriv.append(new)

        # Square without derivative
        input = cur_ecg
        if len(cur_square1) == keep:
             cur_square1.pop(0)
        cur_square1.append(0) # Doesn't matter
        new = square(input, cur_square1)
        cur_square1[-1] = new
        square1.append(new)

        # Square with derivtive
        input = cur_deriv
        if len(cur_square2) == keep:
             cur_square2.pop(0)
        cur_square2.append(0) # Doesn't matter
        new = square(input, cur_square2)
        cur_square2[-1] = new
        square2.append(new)

    plt.figure(figsize=(10, 6))
    #plt.xlabel('Time (s)')
    #plt.plot(x_data, data, label='ECG Data');
    plt.plot(signal, label='Signal', marker='o', markersize=2)
    plt.title('Derivative Filter')
    plt.plot(deriv, label='Derivative', marker='o', markersize=2, color='gold')
    plt.plot(square1, label='Square', marker='o', markersize=2, color='green')
    plt.plot(square2, label='Derivative and Square', marker='o', markersize=2, color='crimson')
    # Show the axis
    plt.axhline(0, color='black', linewidth=.5)
    plt.title('Derivative Filter')
    plt.legend(loc='lower left', framealpha=0.6)
    plt.show()
    
def main():
    print(path.basename(path.normpath(__file__)))
    #test_moving_average()
    test_derivative()

if __name__ == "__main__":
    main()
