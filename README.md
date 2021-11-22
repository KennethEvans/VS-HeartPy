# VS-HeartPy

VS-HeartPy is a Visual Studio project for exploring peak detection in real-time ECG's.The eventual goal is to implement this in an Android app using Java, but exploring seems easier using Python with Numpy and Matplotlib.

The device of interest is the Polar H10, which is capable of producing an ECG at 130 Hz. See <https://github.com/KennethEvans/KE.Net-ECG>.

These ECG's may be irregular. Many heart-rate detection schemes do not work on some of the ECG's considered. That is the reason for implementing another one. 

Most papers on peak detection seem to describe results but not the actual coded algorithms.

HeartPy (<https://github.com/paulvangentcom/heartrate_analysis_python>) is an open-source Python GitHub project and has been installed. It is well developed and has a lot of features. Some of the scripts in this project investigate its use. However, it appears that HeartPy does not give good peak detection for irregular heartbeats.

A Pan Thompkins (<https://en.wikipedia.org/wiki/Pan%E2%80%93Tompkins_algorithm> and <https://courses.cs.washington.edu/courses/cse474/18wi/labs/l8/QRSdetection.pdf>) algorithm has been considered. These authors have done a better job of explaining their implementation than most, but there is still a lot that needs to be figured out. One issue is that this algorithm was designed for 200 Hz, so the low and high-pass filter coefficients are not right for 130 Hz. It was also designed to be coded in assembly language using integer arithmetic. A newer Android application should not be restricted this way.

It should be noted that scipy.signal.find_peaks does seem to find the irregular peaks resonably well, but not perfectly for an irregular heartbeat. It is also designed for post-processing the entire file.