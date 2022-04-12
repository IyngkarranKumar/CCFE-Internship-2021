''' The purpose of this script is to show Iyngkarran how to get 
	JET fast magnetics data.
	Alex Tinguely 2021-06-10
'''

### Import ###

#import numpy as np
import matplotlib.pyplot as plt
import sys


sys.path[:0]=['/jet/share/lib/python'];
from ppf import *
import getdat
sys.path[:0]=['/home/pgpuglia/AEAD/GUI/class'];
from INCAA_functions_all import INCAA_Iant_Vant


### Get data ###

# List of JET pulse numbers
pulses 	= [94700];
pulse 	= pulses[0]; 	# choose pulse

tStart 	= 47; 			# s, start time after which to get data, JET pulses start at 40 s
tStop	= 51; 			# s, end time before which to get data, JET pulses usually end at 60 s


# Here MD = Master Driver, frequency is in kHz
time_MD,frequency_MD,Iant_MD,Vant_MD,time,frequency,Iant,Vant = INCAA_Iant_Vant(pulse);

# Only consider antenna currents > threshold
Ithreshold 	= 2; # A
Iant2 		= np.abs(Iant[:,1]); # get absolute current in antenna 2 (because sometimes antenna 1 doesn't work)
boolI 		= Iant2 >= Ithreshold; #Indices greater than threshold

dt = 5e-7 

#frequency = frequency[boolI]
#time = time[boolI]
length = len(time)

start_bin = ((tStart-40)/20)*length
stop_bin = ((tStop-40)/20)*length    #Assume JET pulse of length 20s starting at 40s 
start_bin = round(start_bin)
stop_bin = round(stop_bin)
 
#frequency =frequency[start_bin:stop_bin]
#time = time[start_bin:stop_bin]

##Interpolating signal
interp_t = np.linspace(tStart,tStop,((tStop-tStart)/dt))
interp_f = np.interp(interp_t,time,frequency)

##Sending to grid
T_bin = 2.048e-3
N_bin = T_bin/dt






fig=plt.figure()

plt.plot(time,frequency)   #Plot amplitude 2 frequency vs time of second antenna
#plt.plot(interp_t,interp_f)
#plt.xlim((tStart,tStop))
plt.xlabel('time (s)');		# x-axis label
plt.ylabel('f (kHz)');				# y-axis label
plt.title(str(pulse));		# title

fig.show()

