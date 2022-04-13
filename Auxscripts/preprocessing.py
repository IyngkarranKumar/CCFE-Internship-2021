"""
The purpose of this script is to provide a few simple functions to be used in the 
(pre?)processing of the fast magnetics data. Functions include spectogram, spectogram colour plot,
Coherence and average Coherence.
"""
import scipy.io as sio
from scipy import signal
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import itertools as it
import cmath
from copy import deepcopy

        
def spectogram(signal,N_bin=int(2**12),Fs=2e6,noverlap=2**10,
               start=40,stop=60,
               mode="complex",pad_to=None):
    
    
    
	if pad_to == None:
		pad_to = N_bin
    
	spect,freqs,time = mlab.specgram(signal,NFFT = N_bin,Fs = 2e6,
                                 noverlap=noverlap, mode=mode,pad_to=pad_to)                      
	
	time = np.linspace(start,stop,time.shape[0]) #Shift to start stop range
	time = time+(N_bin/(2*Fs)) #Shift to end of time bin
	
	return spect,freqs,time


def cross_spectrum(A,B,n,threshold=0):
	"""Generates CSD matrix for two spectogram matrices A and B.

	Parameters:
	----------
	A: Spectogram matrix of signal 1. Here we take spectogram matrix to be a matrix with A[ij]
	the fourier coefficient of jth frequency in ith time bin.

	B: Spectogram matrix of signal 2.

	n: Window size. #Should play around here

	threshold: Set elements less than threshold = 0. Default value set to 0
	as threshold applied when using avgCoherence function.

	Returns:
	-------

	Out: Cross spectral density (CSD)
	"""
			
	#data 	= A*np.conj(B)/np.abs(A)/np.abs(B);
	data 	= A*np.conj(B)/np.abs(np.multiply(A,np.conj(B))); 
	cross 	= deepcopy(data);
	denom 	= np.ones(shape=np.shape(data)); # denominator for averaging

	for i in range(n): 
		#jkth in cross is the sum of j(k-n)th to jkth in data
		cross[:,i+1:] += data[:,:-(i+1)]; # add data to cross by sliding matrix by one time bin
		denom[:,i+1:] += 1; # add 1 to denominator
		
	cross = cross/denom; # take average
	cross = np.abs(cross)  
			
			
            
	return cross

def Coherence(A,B,n=5,threshold = 0.4):
	#FUNCTION VOID
	
	
	COH = np.abs(cross_spectrum(A,B,n=n))
	COH[COH<threshold]=0
                    
	return COH
    
	#COH = np.divide(np.abs(cross_spectrum(A,B,n=n)**2),
	#np.multiply(cross_spectrum(A,A,n=n),cross_spectrum(B,B,n=n)))
    
    
    
    
def avgdCoherence(spectograms,threshold=0.6,n=5,tk=8):
	
	"""Average coherence matrix for probe files, with threshold applied.

	Parameters
	----------------------------
	probe_files: List of probe files
	threshold: Threshold values to apply to avg Coherence matrix
	n: Coherence matrix window size
	tk: Smaller tk means more Coherence matrices calculated


	Returns
	----------------------------
	avgCOH: Average coherence matrix, with threshold applied

	"""
	#freqs,time = spectogram(probe_files[0],start=start,stop=stop)[1:3]   #Freq and time array

	#First coherence matrix
	avgCOH = cross_spectrum(spectograms[0],spectograms[1],n=n)


	pairs = np.array(list(it.combinations(np.arange(0,len(spectograms),1),2)))  #All pairs
	take = np.arange(0,len(pairs),tk) #Take one in every 8
	pairs = pairs[take]
	for pair in pairs[1:]: #Skip first
		COH = cross_spectrum(spectograms[pair[0]],spectograms[pair[1]],n=n)

		avgCOH = avgCOH+COH
		#print(pair)
		
	avgCOH = np.abs(avgCOH)
	avgCOH = avgCOH/len(pairs)
	avgCOH[avgCOH<threshold]=0
	return avgCOH
	

	
	



