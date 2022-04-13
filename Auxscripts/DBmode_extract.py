"""
The purpose of this script is to extract data about individual resonances seen across the 
pulses processed in the UAE DATABASE (UAEDB).
We aim to construct a Dataset in the UAEDB HDF5 File with the following structure:

	|	JPN	|	n	|	m 	| 	t	|	f	|	A	|	p1	|   p2	|				   |   pn	|
	|-------|-------|-------|-------|-------|-------|-------|-------|		...        |--------|
	|		|		|		|		|		|		|		|		|				   |		|
	|		|		|		|		|		|		|		|		|				   |		|
	



JPN - JET pulse number
n - Toroidal mode number 
m - Poloidal mode number
t - Tuple of times. Min,M,Max,t_a, 
	- Min,M,Max are lower quartile, median and upper quartile of the times that comprise mode in question
	- t_A is time at which maximum amplitude achieved
f - Tuple of frequencies. Min,M,Max,f_a, 
	- Min,M,Max are lower quartile, median and upper quartile of the frequencies that comprise mode in question
	- f_A is frequency at which maximum amplitude achieved
A - Max amplitude of mode
pi - Unsure on this still, but perhaps value of ith parameter at times in t?
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py as h5
from copy import deepcopy


from Auxscripts import readUAEDB as rDB
from Auxscripts import plotting as impl



DB_path = "/home/ht2059/Documents/DATABASE_MAIN.h5"

#KEY LISTS
pulses = rDB.pulse_list();pulses = np.array(pulses); pulses.astype(int)

data_keys = rDB.param_list()

NMax  = 8; NMin = -1; NRange = np.arange(NMin,NMax+1,1)
MMax = NMax*3;MMin=NMin; MRange = np.arange(MMin,MMax+1,1)

#FREE PARAMETERS
T_bin = 2.048e-3
df = 488
#0.075
#2100
ext_tgap = 0.0075//T_bin  #Time thresholds. How close must two indices be to be considered part of the same structure?
ext_fgap = 2100//df #Modes must be within 2100Hz of one another to be in same mode

#We define a significant structure as existing for 75 ms or more
SS_time_threshold = 75e-3  #Note that the sig_str parameter has not been optimised yet to pick out just AEs


size_threshold = 100 #To remove noise

check = 500 #How often clear non-SSs


#TEST:
test = True
k = 1
pulses = pulses[:k] 




mode_columns = ["JPN","n","m","t","f","A"]
PoIs = ["bvac","ne0"] #Parameters of interest
cols = mode_columns+PoIs
datfram = pd.DataFrame(columns = cols) #Empty dataframe


num_pulses = 1
pulses = [94700]
for pulse in pulses:
	
	if 1: #Getting data structures
		
		mode_store = pd.DataFrame(columns = cols)
		
		DATABASE = h5.File(DB_path,"r")
	
		sub_db_path = rDB.find_pulse(pulse)[0]

		sub_db = DATABASE[sub_db_path]
	
			
		

		
		
		
		
		
		pulse_data = sub_db[str(pulse)] 
		NZS = pulse_data["NZ indices"];NZS = deepcopy(np.array(NZS));NZS = np.transpose(NZS)
		TMN_vals = pulse_data["TMN spectogram"];TMN_vals = deepcopy(np.array(TMN_vals))
		Shape = pulse_data["Shape"];Shape = np.array(Shape)
		time = pulse_data["Time array (s)"];time = deepcopy(np.array(time))
		param_dict = {}
		for param in PoIs:
			param_dict[param]=np.array(pulse_data[param])
			
		DATABASE.close()
		
		SS_bin_threshold = SS_time_threshold//Shape[1] 
		
		
		
		PMN,avgAMPL = rDB.spect_reconstruct(pulse)[2:4] #PMN and avgAMPL spects
		#PMN_vals = deepcopy(np.array(pulse_data["PMN spectogram"]))
		#avgAMPL = deepcopy(np.array(pulse_data["Averaged Amplitudes"]))
		
		tStart = time[0]; tStop = time[-1]   #Start and end times of PULSE
		

	
	val_dict = dict()  #Empty dictionary to store 
	
	SIG_STRS_ALL=[]
	
	for n in NRange: #Filtering by TMN
		
		
		if 1: #Prepare indice array
			
			print("\n\n\nFiltering for mode number "+str(n))
			booln = TMN_vals==n 
			NZS_n = NZS[booln] #NONZERO indices for some pulse with toroidal mode number n. First columns are ordered frequency bins
			
			#Mode separation
			NZS_n = np.flip(NZS_n,1)
			NZS_argsort = np.argsort(NZS_n[:,0]) #Sorts NZS along
			NZS_n = NZS_n[NZS_argsort] #Now NONZERO indices for some pulse, but first column are time ordered bins
		
		sig_strs = []
		
		if 1: #Time sep first, followed by freq
			"""
			
			for i in range(len(NZS_n)): #Time separation of all indices of this mode number. sig_str is list of lists.
			
				if len(NZS_n) == 0: #No indices for n
					continue
					
				
				if i==0:
					sig_strs.append([])
					sig_strs[-1].append(NZS_n[i])
				
				else:
					if (NZS_n[i][0]-NZS_n[i-1][0])>=ext_tgap:    #If difference greater than threshold 
									sig_strs.append([])    #Add new list
									sig_strs[-1].append(NZS_n[i])    #Append i to that new list
					else: 
						sig_strs[-1].append(NZS_n[i])    #Add indice to existing list
						
			print(len(sig_strs))
			ii=0
			sig_strs = [np.asarray(sig_str) for sig_str in sig_strs]
			
			
			if 1: #Ensuring our list only contains structures existing in time for over 7.5 ms
				pre = len(sig_strs)
				print(type(sig_strs[0]))
				sig_strs = [sig_str for sig_str in sig_strs if np.subtract(np.amax(sig_str[:,0]),np.amin(sig_str[:,0]))>SS_bin_threshold]
				post = len(sig_strs)
				print("From a list of "+str(pre)+" structures, upon applying threshold filter, "+str(post)+" remain")  
				
			
			#sig_strs now contains all significant TIME structures. sig_strs[i] are the indices of ith significant time structure
			#Now separate in frequency space
			
			if 1: #Frequency separation
			
				n_ts = len(sig_strs) #Number of time sig strs
				for t_struct in sig_strs[:n_ts]:
					t_struct = np.array(t_struct)
					f_sort = np.argsort(t_struct[:,1]) #Order by frequencies
					t_struct = t_struct[f_sort] #Time structure now ordered by frequencies. Still arranged as (t,f)
					for i in range(len(t_struct)):
								if i==0:
									sig_strs.append([])
									sig_strs[-1].append(t_struct[i])  #Create bin for first case
								else: 
									if (t_struct[i][1]-t_struct[i-1][1])>=ext_fgap:
										sig_strs.append([])
										sig_strs[-1].append(t_struct[i])
									else:
										sig_strs[-1].append(t_struct[i])
								
					
			sig_strs = sig_strs[n_ts+1:] #Once time structure has been separated in frequency space, can remove original unseparated time structures (first n_ts)
			
			
			
		"""
		
		
		
		
		
		pre = len(sig_strs)
		sig_strs = [np.asarray(sig_str) for sig_str in sig_strs]
		sig_strs = [sig_str for sig_str in sig_strs if np.subtract(np.amax(sig_str[:,0]),np.amin(sig_str[:,0]))>SS_bin_threshold]
		sig_strs = [sig_str for sig_str in sig_strs if len(sig_str)>size_threshold]
		post = len(sig_strs)
		print("From a list of "+str(pre)+" structures, upon applying threshold filter, "+str(post)+" remain")
		#Now sig_strs[i] is ith FREQUENCY AND TIME separated structure
		#All that is left to do is computation of column headings listed in documentation at top, for each sig structure (SS)
		
		SIG_STRS_ALL.append(sig_strs) #Keep list of sig_strs for each n
		
		if 1: #Finding defining properties (t,f,m,A,parameters) of each SS in sig _strs, then storing in datfram
			print("For TMN "+str(n)+" there are "+str(len(sig_strs))+" structures identified to persist for "+str(SS_time_threshold)+" ms or longer")
			for SS in sig_strs:
				time_in = np.array([inds[0] for inds in SS])  # time array of all indices in SS
				frequencies_in = np.array([inds[1] for inds in SS]) #S freq array of all indices in SS
				
						
				tb_min = np.min(time_in);tb_Med = np.quantile(time_in,0.5);tb_max=np.max(time_in)
				fb_min = np.min(frequencies_in);fb_Med = np.quantile(frequencies_in,0.5);fb_max = np.max(frequencies_in)
				t_find = lambda tb: tStart+tb*T_bin
				f_find = lambda fb: (fb*df)/1000     #kHz
				t_min = round(t_find(tb_min),4);t_Med = round(t_find(tb_Med),4);t_max=round(t_find(tb_max),4)
				f_min = round(f_find(fb_min),2);f_Med = round(f_find(fb_Med),2);f_max=round(f_find(fb_max),2)
				#Finds min,med,max of time and freq values
				
				inds_n_SS = (frequencies_in,time_in) #tuple of freq time indices for this n,SS combination
				
				avgAMPL_n_SS = avgAMPL[inds_n_SS]
				A = round(np.max(avgAMPL_n_SS),2) #Max amplitude in this structure
				ind_max = np.argmax(avgAMPL_n_SS)
				t_A = t_find(inds_n_SS[1][ind_max])-tStart;f_A = f_find(inds_n_SS[0][ind_max])  #t_A is time since mode appeared

				t=(t_min,t_Med,t_max,t_A)
				f=(f_min,f_Med,f_max,f_A)
				

				PMN_n_SS = PMN[inds_n_SS] #PMNS for this structure
				ms=[]
				for m in MRange:
					boolm = PMN_n_SS==m
					PMN_n_m_SS = PMN_n_SS[boolm]  
					mlen = len(PMN_n_m_SS)
					if (mlen/len(PMN_n_SS)) > 0.3:
						#If PMN m comprises greater than 20 pc of pmns in structure
						ratio = mlen/len(PMN_n_SS)
						ms.append((m,ratio))
						
				m = tuple(ms) #Tuple of tuples
			
			
			#GETTING PARAMETERS:
			#NEED ONLY GET PARAMETERS OF INTEREST, PoIs defined above

			
				diff = np.abs(time - t_Med)
				p_index = np.argmin(diff) #Index of minimum value of diff
				parameters_ofInt = []
				
				for param in PoIs:
					param_data = param_dict[param]
					param_tMed = param_data[p_index]
					parameters_ofInt.append(param_tMed)
					
				
				#Populate temp_dataframe for a SS
				
				vals = [str(pulse),n,m,t,f,A]+parameters_ofInt #Values for a SS of mode number n
				temp_dict = dict()
				i=0
				for col in cols:
					
					temp_dict[col]=vals[i]
					i+=1
					
				datfram=datfram.append(temp_dict,ignore_index=True)
				
				#print(temp_dict)
				temp_dict.clear()
		
	print("Finished "+str(num_pulses)+" pulse(s)")
	
	if test==True:
		if num_pulses==k:
			break
	num_pulses+=1

					
				
				
def rough_plot():
	for i in range(len(datfram["t"])):
		plt.plot([datfram["t"][i][0],datfram["t"][i][2]],[datfram["f"][i][0],datfram["f"][i][0]])
		plt.show()
			
			
		
		
		
		
			
			
			
			
			
		
		
		
	
	
	
	
	
	
	
	
	
	
	


	
	

