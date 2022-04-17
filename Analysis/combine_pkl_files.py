"""
The purpose of this script is to combine pkl files and store again 
unstable AE databases.


"""

if 1: ##IMPORT
	#import h5py as h5
	import numpy as np
	import pandas as pd
	#import scipy.io as sio
	#from copy import deepcopy
	#import pandas as pd
	#import matplotlib.pyplot as plt
	import pickle
	import sys
	sys.path[:0] = ['/home/atinguel/code'];
	from myFunctions import yymmdd
	
	#from Auxscripts import readUAEDB as rDB
	#from Auxscripts import plotting as impl


if 1: ### Setup ###
	folder 	= '/home/atinguel/work/fromIyngkarranK/';	
	fns 	= ['MEA_220405_0.pkl','MEA_220405_94743.pkl','MEA_220405_95948.pkl','MEA_220405_96866.pkl'];
	fn_new 	= 'Identified_Unstable_modes_AT_'+yymmdd()+'.pkl';
	
if 1: ### Store ###

	for fn in fns:
		temp 	= pd.read_pickle(folder+fn);
		#print(len(temp['n']));
		print(np.min(temp['JPN']));
		print(np.max(temp['JPN']));
		try: 	data 	= pd.concat([data,temp],ignore_index=True);
		except: data 	= temp;


	'''
	data 	= {};  # Create an empty dictionary
	for i in range(len(fns)):
		with open(folder+fns[i],'rb') as f:
			temp 	= pickle.load(f);
			#if i==0: 	data.update(temp);
			#else:
			for key in temp.keys():
				try:	data[key].append(temp[key]);
				except:	data[key] = temp[key]; #pass;
	'''
			
	#for fn in fns:
	#	with open(folder+fn,'rb') as f:
	#		temp.update(pickle.load(f));   # Update contents of f to the dictionary

	with open(fn_new,'wb') as handle:
		pickle.dump(data,handle,protocol=pickle.HIGHEST_PROTOCOL);
