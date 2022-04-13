"""
The purpose of this script is to plot BAEs for Lidia Piron

### NOTE : need python 3.7
module unload python
module load python/3.7

"""

if 1: ##IMPORT
	import h5py as h5
	import numpy as np
	import pandas as pd
	import scipy.io as sio
	from copy import deepcopy
	import pandas as pd
	
	#from Auxscripts import readUAEDB as rDB
	from Auxscripts import readDB as rDB
	from Auxscripts import plotting as impl


if 1: ##SETUP
	DATABASE_I_path = "./DATABASE_BAE.h5";
	pulses 			= rDB.pulse_list(DATABASE_I_path);
	#pulses 			= np.array([95867]);
	
	NMin	= 0;
	NMax 	= 8;  
	MMin 	= NMin;
	MFactor = 2; # MMax = MFactor*NMax
	FMin	= 0; # kHz
	FMax 	= 30; # kHz
	
	savefig 	= 1;
	figfolder	= './figs_BAE/';

if 1: 
	
	for pulse in pulses:
		
		if savefig==1:
			figfile 	= figfolder+'BAE_'+str(pulse);
			savepath_coh= figfile+'_coh.png';
			savepath_amp= figfile+'_amp.png';
			savepath_tmn= figfile+'_tmn.png';
			savepath_pmn= figfile+'_pmn.png';
		else:
			temp 		= 'temp.png';
			savepath_coh=temp;
			savepath_amp=temp;
			savepath_tmn=temp;
			savepath_pmn=temp;
	
		avgCOH,TMN,PMN,avgAMPL,time,freqs = rDB.spect_reconstruct(pulse,DATABASE_I_path);
		
		if 1: # avg coh
			impl.image(avgCOH,time,freqs,min_f=FMin,max_f=FMax,
					   title=str(pulse)+": Average Coherence",savepath=savepath_coh);
		if 1: # avg amp
			impl.image(avgAMPL,time,freqs,min_f=FMin,max_f=FMax,cbarlog=True,
					   title=str(pulse)+": Average Amplitude",savepath=savepath_amp);
		if 0: # tmn	
			impl.image_n(TMN,time,freqs,NMax=NMax-1,NMin=NMin,min_f=FMin,max_f=FMax,
						 title=str(pulse)+": Toroidal Mode Number",savepath=savepath_tmn);
		if 0: # pmn
			impl.image_n(PMN,time,freqs,NMax=MFactor*NMax-1,NMin=MMin,min_f=FMin,max_f=FMax,
						 title=str(pulse)+": Poloidal Mode Number",savepath=savepath_pmn);
	
	
	

	
	
	
	
