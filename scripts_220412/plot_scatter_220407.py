"""
The purpose of this script is to plot histograms of data from the 
unstable AE databases.
Alex Tinguely 220405

"""

if 1: ##IMPORT
	import h5py as h5
	import numpy as np
	import pandas as pd
	import scipy.io as sio
	from copy import deepcopy
	import pandas as pd
	import matplotlib.pyplot as plt
	
	from Auxscripts import readUAEDB as rDB
	from Auxscripts import plotting as impl


if 1: ## SETUP
	folder 			= '/home/atinguel/work/fromIyngkarranK/';
	#DATABASE_II_path= folder+"/Identified_Unstable_modes.pkl"
	DATABASE_II_path= folder+"/Identified_Unstable_modes_AT_220406.pkl"	
	#DATABASE_II_fns = ['MEA_211230_0.pkl','MEA_211230_94743.pkl','MEA_211230_95948.pkl','MEA_211230_96866.pkl'];

	
if 1: ##DATABASE II
	
	db 		= pd.read_pickle(DATABASE_II_path);
	
	if 1: ### Print ###
		print('keys:')
		print(db.keys());
		print('averaging times:');
		print(db['TS'][0]); # averaging times for params
	
	"""
	Many columns in DATABASE II are themselves lists. Following function 
	returns array of data for chosen key and index
	e.g: key = "F", index=0, returns min frequency identified for all modes
	Nans and infs NOT removed however 
	
	keys=
	['JPN', 'n', 'm', 'T', 'F', 'A', 'NBLM', 'PTOT', 'TE0', 'angf_cxd6',
       'angf_cxg6', 'angf_cxh4', 'angf_cxhm', 'angf_cxkm', 'btnm', 'bvac',
       'elon', 'fD_kt5p', 'fH_ks3b', 'fH_kt5p', 'fHe3_kt5b', 'fHe4_kt5b',
       'fT_kt5p', 'gradn95', 'ipla', 'ka2-1', 'ka2-2', 'ka2-3', 'ka2-4',
       'ka2-5', 'ka3_ccds', 'n95', 'ne0', 'q0', 'q95', 's95', 'tril', 'triu',
       'TS']
	
	"""
	
	t_start	= rDB.get(db,key="T",index=0); # s, start time of mode
	t_mid 	= rDB.get(db,key="T",index=1); # s, median time
	t_end 	= rDB.get(db,key="T",index=2); # s, end time
	t_0 	= rDB.get(db,key="T",index=3); # s, time at start of pulse
	dt_Amax = rDB.get(db,key="T",index=4); # s, dt to reach max amplitude
	t_Amax 	= t_0+dt_Amax;
	
	f_min 	= rDB.get(db,key="F",index=0); # kHz, min frequency of mode
	f_med 	= rDB.get(db,key="F",index=1); # kHz, median freq
	f_max 	= rDB.get(db,key="F",index=2); # kHz, max freq
	f_Amax 	= rDB.get(db,key="F",index=3); # kHz, freq at max amplitude
	
	A_start = rDB.get(db,key="A",index=0); # T/s? amplitude at t_start
	A_mid 	= rDB.get(db,key="A",index=1); # amp at t_mid
	A_end 	= rDB.get(db,key="A",index=2); # amp at t_end
	Amax 	= rDB.get(db,key="A",index=3); # maximum amplitude
	
	g_sm 	= (A_mid-A_start)/(t_mid-t_start);
	g_me 	= (A_end-A_mid)/(t_end-t_mid);
	g_se 	= (A_end-A_start)/(t_end-t_start);
	g_sa 	= (Amax-A_start)/(t_Amax-t_start);
	
	n		= rDB.get(db,key='n'); # TMN
	m		= rDB.get(db,key='m'); # PMN
	
	fH 		= rDB.get(db,'fH_kt5p');
	fD 		= rDB.get(db,'fD_kt5p');
	fT		= rDB.get(db,'fT_kt5p');
	fHe3 	= rDB.get(db,'fHe3_kt5b');
	fHe4	= rDB.get(db,'fHe4_kt5b');

	#FTOT 	= 1*fH + 1*fD + 1*fT + 2*fHe3 + 2*fHe4;
	Aeff 	= 1*fH + 2*fD + 3*fT + 3*fHe3 + 4*fHe4;
	
	
if 1: ### Plot ###

	
	x 		= g_se; 	xlabel 	= 'growth rate (end - start)';
	x 		= g_sa; 	xlabel 	= 'growth rate (max - start)';
	
	'''	
	boolNaN = np.isnan(x);
	boolInf = np.isinf(x);
	boolTot = np.logical_and(~boolNaN,~boolInf);
	bins 	= np.linspace(np.min(x[boolTot]),np.max(x[boolTot]),26);
	bins 	= np.linspace(0,100,101);
	'''
	
	xlabel = 'beta_N';
	xlabel 	= 'JPN';
	xlabel 	= 'TE0';
	xlabel 	= 'ne0';
	xlabel 	= 'q0';
	xlabel 	= 'q95';
	x	 	= rDB.get(db,xlabel);
	x 		= Aeff; 	xlabel 	= 'Aeff';
	#x 		= Amax; 	xlabel 	= 'amp';
	#x 		= f_med; 	xlabel = 'f (kHz)';
	#x 		= n; 		xlabel = 'n';
	#x	 	= rDB.get(db,'PTOT')*1e-6; xlabel = 'ICRH (MW)';
	#x	 	= rDB.get(db,'NBLM')*1e-6; xlabel = 'NBI (MW)';
	
	#ylabel 	= 'ka3_ccds'; #y	 = rDB.get(db,ylabel);
	#ylabel 	= 'ka2-1';		
	#ylabel 	= 'ka2-2';
	#ylabel 	= 'ka2-3';
	#ylabel 	= 'ka2-4';
	#ylabel 	= 'ka2-5';		
	#y	 	= np.abs(rDB.get(db,ylabel));
	
	y 		= Amax; ylabel 	= 'max amp';
	
	# Frequency filter
	f_filter= [];
	f_filter= [0,100];
	f_filter= [100,300]; # kHz
	f_filter= [0,100,300,500]; # kHz
	
	if 1: ### Scatter ###
		fig 	= plt.figure();
		alpha 	= 0.2;
		
		# Filter
		if len(f_filter)>1:
			for i in range(len(f_filter)-1):
				f0 		= f_filter[i];
				f1 		= f_filter[i+1];
				boolF 	= np.logical_and(f_med>=f0,f_med<=f1);	
				#boolF 	= np.logical_and(f_min>=f0,f_max<=f1);	
				plt.plot(x[boolF],y[boolF],'o',alpha=alpha,label=str(f0)+'-'+str(f1)+' kHz ['+str(sum(boolF))+']');
		else:
			plt.plot(x,y);
			
		plt.xlabel(xlabel);
		plt.ylabel(ylabel);
		plt.legend();
		plt.grid();
		fig.tight_layout();
		fig.show();
		
	
