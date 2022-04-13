''' The purpose of this function is to make histrograms from the database
	Alex Tinguely 21-07-26
'''

if 1: ### Import ###
	import h5py as h5
	import numpy as np
	import matplotlib.pyplot as plt
	from random import randint
	
	
if 1: ### Setup ###
	filename= "/home/ht2059/Documents/first_run_failed/DATABASE_2_bad.h5";
	filename= "/home/ht2059/Documents/DATABASE_2.h5";
	
	fmin 	= 100e3; # Hz
	fmax 	= 300e3; # Hz
	
if 1: ### Read database ###
	
	db 		= h5.File(filename,'r');
	pulses 	= list(db.keys());
	params 	= list(db[pulses[0]].keys());
	
	
	
	#QUICK CASE
	
	pulses = pulses[:round(0.2*len(pulses))]
	
	print('...number of pulses: '+str(len(pulses)));
	
	
if 1: ### Histograms ###
	
	#for param in params:
	param = 'ipla';
	param = 'ne0';
	param = 'TMN spectogram';
	param = 'PTOT';
	param = 'NBLM';
	param = 'TE0';
	param = 'bvac';
	param = 'q0';
	param = 'q95';
	
	param = "TMN spectogram"
	
	if 1: 
		temp = np.array([]);
		print('... collecting data for '+param);
		
		for pulse in pulses:
			
			#if 'PMN' in param or 'TMN' in param:
			#	f 		= db[pulse]['Frequency array (Hz)'][:];
			#	boolf 	= np.logical_and(f>=fmin,f<=fmax);
			#	temp = np.concatenate((temp,db[pulse][param][boolf]));
			#else:
			temp = np.concatenate((temp,db[pulse][param]));
		
		temp = temp[~np.isnan(temp)];
		temp = temp[~np.isinf(temp)];
		
		if 'PMN' in param or 'TMN' in param:
			bins = np.arange(min(temp)-0.5,max(temp)+1.5,1);
			
		else:
			bins = 21;
			
			
		N,bins,patches = ax.hist(temp,bins=bins);
			
		
		colours = []

		for i in range(len(bins)):
			colours.append('#%06X' % randint(0, 0xFFFFFF))
			
		for i in range(len(bins)-1):
			patches[i].set_facecolor(colours[i])
			
		
		
		
		
		
		plt.xlabel(param);
		plt.ylabel('counts');
		
		plt.title('number of datapoints: '+str(len(temp)));
		plt.yscale('log');
		plt.xticks(ticks = np.arange(min(temp)-1,max(temp)+1,1))  
		
		ax = plt.subplot(111)
		plt.show(); #Show axes
			
