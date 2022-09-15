import os,sys
from scipy import stats
import numpy as np
import bottleneck as bn
from scipy import signal
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt 
#____________________________________5 Point Difference_____________________________________________________________________________
def pointDiff5(signal):
   result = np.zeros(len(signal))
   result[0] = signal[0]
   result[1] = signal[1]
   result[len(signal)-1] = signal[len(signal)-1]
   result[len(signal)-2] = signal[len(signal)-2]
   for n in range(2,len(signal)-2):
       result[n] = signal[n+2]+2*signal[n+1]-2*signal[n-1]-signal[n-2]
   r = np.multiply(result,(1/8.0)*256)
   return r
#__________________________________Set Threshold_____________________________________________________________________________________
def setThreshold(smoothed,N):
  # get the variance of the signal 
  # to check how often does the signal change to avoid considering low values 
  # set the variance as the window
  # save the peak of every window
  # get the average peak = Threshold
  peaks = []
  window = []
  j=0
  for i in range(len(smoothed)):
    window.append(smoothed[i])
    if( i % int(np.var(smoothed)) == 0):
      window = np.array(window)
      if(len(peaks)>1 and peaks[j-1]>np.amax(window)):
        peaks[j-1]= np.amax(window)
      else:
        peaks.append(np.amax(window))
        j+=1
      window = []
     
  T = sum(peaks)/float(len(peaks))
  return T
#________________________________Moving Average Window_____________________________________________________________________________________  
def movingAverageWindow(signal,N):
   result = []
   for i in range(len(signal)):
     result.append(np.sum(signal[i-N:i]))
   return np.multiply(result,(1/float(N)))
#___________________________________Filter Signal_____________________________________________________________________________________________
def filterSignal(array):
   #sampling rate
   fs = 256
   # Filter the signal using notch filter
   b, a = signal.iirnotch(50,30, fs)
   afternotch =  lfilter(b, a, array)

   # Filter the signal using BandPass Filter
   nyquist_freq = 0.5 * fs
   low = 0.5 / nyquist_freq
   high = 45 / nyquist_freq
   b, a = butter(1, [low, high], btype="band")
   afterbandPass = lfilter(b, a, afternotch)
  
   return afterbandPass
#_____________________________________QRS Detection_____________________________________________________________________________________________  
def QRS(signal,N):

    #filter signal using notch filter then bandpass filter
    filtered = filterSignal(signal)

    """# Plot signal before and after noise filtering
    fig = plt.figure()
    
    x1=fig.add_subplot(211)
    x1.plot(signal[0:1999])
    x1.title.set_text('ECG Signal Before Noise Filtering')

    x2=fig.add_subplot(212)
    x2.plot(filtered[0:1999])
    x2.title.set_text("ECG Signal After Noise Filtering")
    fig.tight_layout()

    # shift subplots down
    fig.subplots_adjust(top=0.85)
    plt.show()"""
    
    #calculate 5 point difference using the defined function above
    difference = pointDiff5(filtered)
 
    
    # Next, square the values
    squared = difference**2

    # Moving Average Window using defined function above
    smoothed = movingAverageWindow(squared,N)

    # Set Treshold using the defined function above
    T = setThreshold(smoothed,N)
    
    # Detect R waves 
    Timestamps = []
    window = []
    indices = []
    # Set a window of (250) 
    # 250 is chosen using trial and error
    # Take the max of every window, if it is greater than the threshold
    # Add its index in the array of smoothed signal in the timestamps array
    # setting a window helps to avoid multiple stars in a single peak
    for i in range(len(smoothed)):
      window.append(smoothed[i])
      indices.append(i)
      if(i%250==0):
        max = np.max(window)
        ind = np.argmax(window)
        if(max>=T):
          Timestamps.append(indices[ind])
        window=[]
        indices=[]


    
    RR = []
    # Get RR Intervals
    for i in range(1,len(Timestamps)):
      RR.append(Timestamps[i]-Timestamps[i-1])


   # Create indices for stars '*' in the plot
    markers = []
    for i in range(len(Timestamps)):
      if(Timestamps[i]<=1999):
        markers.append(Timestamps[i])

    # plotting the points  
    #plt.plot(smoothed[0:1999], marker="*",markevery=markers,markerfacecolor='red',markersize=10) 
    x = np.array(range(len(RR)))
    print(x)
    y = np.multiply(RR,1000/256)
    plt.plot(x,y) 
    # naming the x axis 
    plt.xlabel('Beat Number') 
    # naming the y axis 
    plt.ylabel('RR Interval in msec') 
    # giving a title to my graph 
    plt.title('RR Intervals') 
    #plt.xticks(np.arange(0,2250,250))
    #function to show the plot 
    plt.show()
    
    return Timestamps, RR
#_______________________________MAIN____________________________________________________________________________________________________________

with open('C:/Users/yasmi/Desktop/Semester 10/Biomedical/Assignment 1_31990/dataN.txt') as file_in:
      array = []
      for line in file_in:
        array.append(float(line))
QRS(array,1)

