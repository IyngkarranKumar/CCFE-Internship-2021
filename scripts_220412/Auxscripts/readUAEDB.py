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


def help_UAEDB():
	print('''
	The following provides a description of the Unstable Alfven Eigenmode Database, covering the general database file structure as well as 
	some further general information
	
	file type: HDF5		file size = 9.6 GB	
	
	Pulses in database satisfy following criteria:
	
	ICRH power > 500 kW
	time of ICRH power > 100ms
	
	The Database is accessed in the following manner:
	
	Import the h5py module first to serve as an interface between the terminal and the HDF5 file. Then using the h5py.File() method, import the database.
	
	path = "/home/ht2059/Documents/DATABASE_main.h5"
	DATABASE = h5py.File(path,"r")
	
	IMPORTANT: DO NOT open the file to write ("w"). This will overwrite any data currently saved.
	
	The HDF5 acts as a nested dictionary. Calling   list(DATABASE.keys())   will return 5 keys:
	
	sub DATABASE_1			Contains pulse range 93130-93796 (80 pulses) and 93975-94577 (250 pulses)
	sub DATABASE_2			Contains pulse range 94578-95419 (400 pulses) 			
	sub DATABASE 3			Contains pulse range 93795-93977 (70 pulses) and 95421-996297 (400 pulses)
	sub DATABASE_4			Contains pulse range 96298-996848 (300 pulses) and 97499-97895 (200 pulses)
	sub DATABASE_5			Contains pulse range 96849-97498 (300 pulses) and 97896-98005 (60 pulses)
	
	Overall, the database was compiled from 2054 pulses
	
	
	One may access each individual database. An example of how to do so is below:
	
	: sub DATABASE_5 = DATABASE["sub DATABASE_5"]
	
	
	
	Upon entering sub Databases, type the following command (modified to the sub database in question) to see the list of pulses in the file:
	
	: pulses = list(sub_DATABASE_5.keys())
	: print(pulses)
	
	
	
	For each pulse, there are roughly 30 plasma and operational parameters time-series associated.
	

	
	
 
	
	
	
	
	
	
	'''
	)


def pulse_list():
	"""
	Return list of pulse keys
	"""
	
	filename = DATABASE_I_path
	DATABASE = h5.File(filename,"r")
	
	pulses = []
	
	for key in list(DATABASE.keys()):
		keys_subdb = list(DATABASE[key].keys())
		pulses=pulses+keys_subdb
		
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
	pulses = list(sub_db.keys())
	pulse = sub_db[str(pulses[0])]
	param_list = list(pulse.keys())
	
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
	

def find_pulse(pulse):
	
	"""
	Returns the sub database that a pulse is stored in, as well as its index location in the list of pulses in this sub database.
	"""
	pulse_data = sio.loadmat("/home/atinguel/work/forIyngkarranK/pulses_w_ICRH_93000-98006_Pthreshold500kW_Tthreshold500ms")
	pulses = pulse_data["pulses"][0]
	
	filename = DATABASE_I_path
	DATABASE = h5.File(filename,"r")
	sub_db_keydict = {}

	sub_db_keydict["sub_db1_keys"] = np.array(list(DATABASE["sub DATABASE_1"].keys()))
	sub_db_keydict["sub_db2_keys"] = np.array(list(DATABASE["sub DATABASE_2"].keys()))
	sub_db_keydict["sub_db3_keys"] = np.array(list(DATABASE["sub DATABASE_3"].keys()))
	sub_db_keydict["sub_db4_keys"] = np.array(list(DATABASE["sub DATABASE_4"].keys()))
	sub_db_keydict["sub_db5_keys"] = np.array(list(DATABASE["sub DATABASE_5"].keys()))
	DATABASE.close()
	
	sub_dbs = ["sub_db1_keys","sub_db2_keys","sub_db3_keys","sub_db4_keys","sub_db5_keys"]

	i=1
	for sub_db_keys in sub_dbs:
		
		
		arr = np.where(sub_db_keydict[sub_db_keys]==str(pulse))
		arr = arr[0]
		
		if len(arr)!=0:
			arr = list(arr)

			return "sub DATABASE_"+str(i),arr
		else:
			
			if i==5:
				print("ERROR: Pulse "+str(pulse)+" entered not found in any of the sub databases")
				return sub_db_keydict
				
			else:
				pass
				
		i+=1
			
			
def spect_reconstruct(pulse):
	"""
	Reconstructs spectrograms stored in the database for a given pulse.
	Returns in order avgCOH,TMN,PMN,avgAMPL,time,freqs
	"""
	

	DATABASE_path = DATABASE_I_path
	DATABASE = h5.File(DATABASE_path,"r")
	
	sub_db_path = find_pulse(pulse)[0]

	sub_db = DATABASE[sub_db_path]
	
	data = sub_db[str(pulse)]
	
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
		 

def get(data,key,index=0):
	
	key_data 	= data[key];
	param_data	= [];
	for i in range(len(key_data)):
		try: 	
			temp	= key_data.iloc[i];
			if isinstance(temp, (list, tuple, np.ndarray)):
				param_data.append(temp[index]);
			else: 
				param_data.append(temp);
		except: 
			param_data.append(np.nan);
	#print(param_data)	
	param_data = np.array(param_data);
	return param_data;

#
			
	

		
	
	
	


		
	
	


