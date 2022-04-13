''' The  The purpose of this script is to get ICRH frequencies
	Alex Tinguely 2021-06-10
'''

if 1: ### Import ###

	import numpy as np
	import matplotlib.pyplot as plt
	import sys
	sys.path[:0]=['/jet/share/lib/python'];
	from ppf import *
	import getdat
	from copy import deepcopy


def get_ICRH(pulse,tStart,tStop,size):
	pulse=pulse

	dda 	= 'ICRH';
	antennas= ['A','B','C','D','E']; # note E = ITER-like antenna
	
	
	o_timebase= np.linspace(tStart,tStop,size); # interpolation timebase. size should be avgCOH.shape[1]
	timebase = deepcopy(o_timebase)
	

	F 	= []; # list for frequencies
	FU 	= []; # list for unique frequencies
	P 	= []; # list for powers
	empty_sig=0
	for antenna in antennas:
		
		# Get frequency
		f,x,tf,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
		ppfdata(pulse,'ICRH','FR'+antenna,seq=0,uid='jetppf',device='JET',\
		fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);

		# Get power 
		p,x,tp,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
		ppfdata(pulse,'ICRH','PRF'+antenna,seq=0,uid='jetppf',device='JET',\
		fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);
		
		
		try:
			
			F.append(np.interp(timebase,tf,f,left=np.nan,right=np.nan));
			FU.append(np.unique(f));
			P.append(np.interp(timebase,tp,p,left=0,right=0));
			
		except:
			empty_sig+=0
			F.append(np.zeros(timebase.size));
			FU.append(np.unique(f));
			P.append(np.zeros(timebase.size));
			
		
	if empty_sig==0:
		print(" \n Signals returned for all antennas")
	else:
		print(" \n No signal returned for "+str(empty_sig)+" antennas")
			


	# Determine beating frequencies
	Fthreshold 	= 1e3; # Hz, threshold for difference in frequency
	Pthreshold 	= 0.1e6; # W, threshold in power

	TB 			= []; # list for times
	FB 			= []; # list for beating frequencies

	for i in range(len(antennas)-1): # for all antennas except the last
		for j in range(i+1,len(antennas)): # for all other antennas
			
			if np.abs(FU[i]-FU[j]) < Fthreshold: continue; # if antennas have same frequency 
			
			boolP 	= np.logical_and(P[i] >= Pthreshold,P[j] >= Pthreshold); # get powers over threshold
			TB.append(timebase[boolP]); # get relevant timebase
			FB.append(np.abs(F[i][boolP]-F[j][boolP])); # get beating frequency
			
			
	# Determine aliasing frequencies
	Fsampling 	= 2e6; 	 # Hz or samples/second, sampling frequency
	FNyquist 	= Fsampling/2; # Hz, Nyquist frequency
	#Pthreshold 	= 0.1e6; # W, threshold in power

	TA			= []; # list for times
	FA 			= []; # list for aliasing frequencies

	for i in range(len(antennas)-1): # for all antennas except the last
		
		boolP 	= P[i] >= Pthreshold; # get powers over threshold
		FAL		= F[i] % Fsampling; # Hz, lower aliasing frequency
		FAU 	= Fsampling - FAL; # Hz, upper aliasing frequency
		
		TA.append(timebase[boolP]); # add relevant timebase
		FA.append(FAL[boolP]); # add lower aliasing frequency
		
		if sum(FAU-FAL) > 0:
			TA.append(timebase[boolP]); # add relevant timebase
			FA.append(FAU[boolP]); # add upper aliasing frequency
	
	for	k in range(len(FA)):	
		try:
			FA[k] =  np.interp(o_timebase,TA[k],FA[k])
		except:
			#If interpolation fails due to being out of range, set to 0
			FA[k] = np.zeros(o_timebase.size)
	for	k in range(len(FB)):
		try:
			FB[k] = np.interp(o_timebase,TB[k],FB[k])
		except:
			FB[k] = np.zeros(o_timebase.size)
		
			
	return FB,FA,TB,TA



