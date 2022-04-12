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

tStart 	= 40; 			# s, start time after which to get data, JET pulses start at 40 s
tEnd	= 60; 			# s, end time before which to get data, JET pulses usually end at 60 s


# Here MD = Master Driver, frequency is in kHz
time_MD,frequency_MD,Iant_MD,Vant_MD,time,frequency,Iant,Vant = INCAA_Iant_Vant(pulse);

# Only consider antenna currents > threshold
Ithreshold 	= 2; # A
Iant2 		= np.abs(Iant[:,1]); # get absolute current in antenna 2 (because sometimes antenna 1 doesn't work)
boolI 		= Iant2 >= Ithreshold; #Indices greater than threshold

dt = 5e-7 

Aant = np.abs(np.multiply(Iant,Vant))**0.5
print(len(time))
print(len(Aant))

### Plot data ###
"""
fig 	= plt.figure();

plt.plot(time[boolI],frequency[boolI]);
plt.xlim((tStart,tEnd));	# x-axis limits
plt.xlabel('time (s)');		# x-axis label
plt.ylabel('f (kHz)');				# y-axis label
plt.title(str(pulse));		# title

fig.show(); 				# show figure
"""
