''' The purpose of this script is to get ICRH frequencies
	Alex Tinguely 2021-06-10
'''

if 1: ### Import ###

	import numpy as np
	import matplotlib.pyplot as plt
	import sys
	sys.path[:0]=['/jet/share/lib/python'];
	from ppf import *
	import getdat


if 1: ### Setup ###

	# List of JET pulse numbers
	pulses 	= [94700];
	pulse 	= pulses[0]; 	# choose pulse

	dda 	= 'ICRH';
	antennas= ['A','B','C','D','E']; # note E = ITER-like antenna

	tStart 	= 40; # s
	tEnd 	= 60; # s
	dt 		= 1e-3; # s, time resolution
	timebase= np.arange(tStart,tEnd,dt); # interpolation timebase

	F 	= []; # list for frequencies
	FU 	= []; # list for unique frequencies
	P 	= []; # list for powers


if 1: ### Get data ###

	for antenna in antennas:
		
		# Get frequency
		f,x,tf,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
		ppfdata(pulse,'ICRH','FR'+antenna,seq=0,uid='jetppf',device='JET',\
		fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);

		# Get power 
		p,x,tp,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
		ppfdata(pulse,'ICRH','PRF'+antenna,seq=0,uid='jetppf',device='JET',\
		fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);

		F.append(np.interp(timebase,tf,f,left=np.nan,right=np.nan));
		FU.append(np.unique(f));
		P.append(np.interp(timebase,tp,p,left=0,right=0));


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


### Plot data ###

fig,axs	= plt.subplots(3,1,sharex=True);

for k in range(len(antennas)):
	
	ax 	= axs[0]; # get first subplot axes
	ax.plot(timebase,F[k]*1e-6,label=antennas[k]);
	ax.set_ylabel('f (MHz)');							
	ax.set_title(str(pulse));							
	ax.legend();

	ax 	= axs[2]; # get third subplot axes
	ax.plot(timebase,P[k]*1e-6);
	ax.set_xlabel('t (s)');		
	ax.set_ylabel('P (MW)');	
	
for k in range(len(TB)):
	ax 	= axs[1]; # get second subplot axes
	ax.plot(TB[k],FB[k]*1e-3);
	ax.set_xlabel('t (s)');		
	ax.set_ylabel('df (kHz)');			

fig.show(); 				# show figure
	
