"""
The purpose of this script is to set up a testing environment for accessing both UAE databases
Examples of function use are provided below


"""

if 1: ##IMPORT
	import h5py as h5
	import numpy as np
	import pandas as pd
	import scipy.io as sio
	from copy import deepcopy
	import pandas as pd
	
	from Auxscripts import readUAEDB as rDB
	from Auxscripts import plotting as impl


if 1: ##SETUP
	DATABASE_I_path = "/home/ht2059/Documents/Final Scripts and Databases/DATABASE_MAIN.h5"
	DATABASE_II_path = "/home/ht2059/Documents/Final Scripts and Databases/Identified_Unstable_modes.pkl"
	
	pules = rDB.pulse_list()
	params = rDB.param_list()

if 1: #SAMPLE PULSE - examples
	NMax = 8; NMin=-1
	pulse = 96851
	
	avgCOH,TMN,PMN,avgAMPL,time,freqs = rDB.spect_reconstruct(pulse)
	
	#impl.image(avgCOH,time,freqs,min_f=0,max_f=200,title="96851 Coherence")
	
	#impl.image(avgAMPL,time,freqs,min_f=0,max_f=200,cbarlog=True,title="96851 Amplitude")
	
	#impl.image_n(TMN,time,freqs,NMax=NMax-1,NMin=-1,min_f=0,max_f=200,title="96851 TMN")
	
	#impl.image_n(PMN,time,freqs,NMax=3*NMax-1,NMin=-1,min_f=0,max_f=200,title="96851 PMN")
	
	#Explore time arrays:
	
	DATABASE_I = h5.File(DATABASE_I_path,"r")  #pls pls pls do not open anything to write, only append. The DB will be wiped and you will have a slightly frustrated intern on your hands
	sub_db,index = rDB.find_pulse(pulse)
	pulse_entry = DATABASE_I[sub_db][str(pulse)]
	
	pulse_entry_B0 = np.array(pulse_entry["bvac"])
	
	
if 1: ##DATABASE II - examples
	
	DATABASE_II = pd.read_pickle(DATABASE_II_path)
	
	"""
	Many columns in DATABASE II are themselves lists. Following function returns array of data for chosen key and index
	
	e.g: key = "F", index=0, returns min frequency identified for all modes
	
	Nans and infs NOT removed however 
	"""
	mode_end_times = rDB.parray(key="T",index=2)
	
	mode_min_fs = rDB.parray(key="F",index=0)
	
	mode_B0_100ms = rDB.parray("bvac",index=2); mode_B0_250ms = rDB.parray("bvac",index=1)
	
	
	
	#Certain pulse
	
	
	
	
	
	
