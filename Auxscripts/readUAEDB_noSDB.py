"""
The purpose of this script is to provide a set of functions to interface with the UAE DATABASE (UAEDB).
The UAEDB is a hdf file, with the first set of keys being individual pulse numbers.
One can access the data stored in a certain pulse (94700 for example) in the following manner:

import h5py as h5

store = h5.File("DATABASE.h5","r")
data = store["94700"]

Data can then be extracted from the 94700 pulse in a similar method. data.keys() will show the available key names.

avgCOH - Averaged Coherence spectogram
TMN - Toroidal mode number spectogram
PMN - Poloidal mode number spectogram
avgAMPL - Averaged Fourier coefficients amplitude spectogram
"""

#-----IMPORT

import h5py as h5
import numpy as np
import pandas as pd
import scipy.io as sio
from copy import deepcopy
import pandas as pd



#For pulses stored in Database
pulse_data = sio.loadmat("/home/atinguel/work/forIyngkarranK/pulses_w_ICRH_93000-98006_Pthreshold500kW_Tthreshold500ms")
pulses = pulse_data["pulses"][0]
tStarts = pulse_data["tStart"][0]
tStops = pulse_data["tEnd"][0]

DATABASE_I_path = "/home/ht2059/Documents/Final Scripts and Databases/DATABASE_MAIN.h5"
DATABASE_II_path = "/home/ht2059/Documents/Final Scripts and Databases/Identified_Unstable_modes.pkl"





def pulse_list():
	"""
	Return list of pulse keys
	"""
	
	filename = DATABASE_I_path
	DATABASE = h5.File(filename,"r")
	
	pulses = list(DATABASE.keys())
	
	
		
	pulses = set(pulses)
	pulses = list(pulses)
	pulses = np.array(pulses)
	pulses = pulses.astype(int)
	pulses.sort()
	
	DATABASE.close()
		
	return pulses
	
def param_list():
	"""
	Return list of data keys
	"""
	
	filename = DATABASE_I_path
	DATABASE = h5.File(filename,"r")
	
	sub_db = DATABASE["sub DATABASE_1"]
	pulses = pulse_list()
	pulse_data = DATABASE[str(pulses[0])]
	param_list = list(pulse_data.keys())
	
	return param_list
	
	
def units():
	"""
	Return dict of units for some parameters
	"""
	d=dict()
	d["bvac"]="T";d["ipla"]="A";d["ne0"]="m-3";d["TE0"]="eV";d["NBLM"]="W";d["PTOT"]="W"
	d["n95"]="m-3";d["bvac"]="T";d["angf_cxd6"]="rads/s";d["angf_cxkm"]="rads/s"
	d["angf_cxg6"]="rads/s";d["angf_cxh4"]="rads/s";d["angf_cxhm"]="rads/s";
	
	return d
	


			
			
def spect_reconstruct(pulse):
	"""
	Reconstructs spectrograms stored in the database for a given pulse.
	Returns in order avgCOH,TMN,PMN,avgAMPL,time,freqs
	"""
	

	DATABASE_path = DATABASE_I_path
	DATABASE = h5.File(DATABASE_path,"r")
	

	data = DATABASE[str(pulse)]
	
	NZS = data["NZ indices"]
	NZS = np.array(NZS)
	NZS = (NZS[0],NZS[1])
	shape = np.array(data["Shape"])
	scaff = np.zeros(shape)
	avgCOHvals = np.array(data["Averaged Coherence"])
	TMNvals = np.array(data["TMN spectogram"])
	PMNvals = np.array(data["PMN spectogram"])
	avgAMPLvals = np.array(data["Averaged Amplitudes"])
	time = np.array(data["Time array (s)"])
	freqs = np.array(data["Frequency array (Hz)"])
		
		
		
	
	#avgCOH

	avgCOH = deepcopy(scaff)
	avgCOH[NZS] = np.array(avgCOHvals)
	
	#TMN

	TMN = deepcopy(scaff)
	TMN[:]=np.nan
	TMN[NZS] = np.array(TMNvals)
	
	#PMN

	PMN = deepcopy(scaff)
	PMN[:]=np.nan
	PMN[NZS] = np.array(PMNvals)
	
	#avgAMPL

	avgAMPL = deepcopy(scaff)
	avgAMPL[NZS] = np.array(avgAMPLvals)
	
	#time 

	time = np.array(time)
	
	#freqs

	freqs = np.array(freqs)
	
	DATABASE.close()
	

	
	return avgCOH,TMN,PMN,avgAMPL,time,freqs
	
	

def parray(key,index):
	"""
	Many columns in DATABASE II are themselves lists. This function returns array of data for chosen key and index
	
	e.g: key = "F", index=0, returns min frequency identified for all modes
	
	Nans and infs NOT removed however
	"""
	
	DATABASE_II = pd.read_pickle(DATABASE_II_path)
	
	key_data = DATABASE_II[key]
	
	param_data=[]
	for i in range(len(key_data)):
		plist = key_data.iloc[i]
		poi = plist[index]
		param_data.append(poi)
		
	param_data = np.array(param_data)
	return param_data
		 



#
			
	

		
	
	
	


		
	
	


