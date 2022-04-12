"""
The purpose of this script is to extract data about individual modes seen across the 
pulses processed in the UAE DATABASE (UAEDB).
We aim to construct a DataFRame in the UAEDB HDF5 File with the following structure:

	|	JPN	|	n	|	m 	| 	t	|	f	|	A	|	p1	|   p2	|				   |   pn	|
	|-------|-------|-------|-------|-------|-------|-------|-------|		...        |--------|
	|		|		|		|		|		|		|		|		|				   |		|
	|		|		|		|		|		|		|		|		|				   |		|
	
JPN - JET pulse number
n - Toroidal mode number (array). If toroidal mode number n comprises over 60% of indices in mode, included (as well as percentage of indices that have this n).
m - Poloidal mode number (array). Same as above, except threshold to be included in 25%
t - Tuple of times. t_min,t_med,t_max,t_A and tStart, 
	- Min,M,Max are lower quartile, median and upper quartile of the times that comprise mode in question
	- t_A is time for pulse to first reach max amplitude
	- tStart is starting time of pulse. Useful for calculations of mode time interval length
f - Tuple of frequencies. f_min,f_med,f_max,f_A
	-Defined analogously to those above
A - Average amplitude of mode
pi - ith parameter in PoI list averaged over mode and 0.5/0.25.0.1 s prior to mode first being seen


I use the term `structures' often in this script. These are groups of indices in close proximity in FT space that COULD be an unstable mode. 
Many structures are removed during the filtering process
"""


if 1:
	import numpy as np
	import matplotlib.pyplot as plt
	import pandas as pd
	import h5py as h5
	from copy import deepcopy

	from Auxscripts import merge as merge
	from Auxscripts import readUAEDB as rDB
	from Auxscripts import plotting as impl
	
	import sys
	sys.path[:0] = ['/home/atinguel/code'];
	from myFunctions import yymmdd


if 1: 
	#KEY LISTS
	pulses 		= rDB.pulse_list();
	pulses 		= np.array(pulses); 
	pulses.astype(int);

	data_keys 	= rDB.param_list()

	NMax  		= 8; 
	NMin 		= -1; 
	NRange 		= np.arange(NMin,NMax+1,1)
	MMax 		= NMax*3;
	MMin 		= NMin; 
	MRange 		= np.arange(MMin,MMax+1,1)
	
	mode_columns = ["JPN","n","m","T","F","A"] 


if 1:  #FILTERS AND FREE PARAMETERS
	# Regarding setting of free parameters(thresholds), I have found from a weeks worth of playing around that there
	# does not seem to be a `perfect set' of parameters that works best for all pulses. 
	# Maybe parameters CAN be optimised for individual frequency bands? (only had time to experiment with this to a basic degree)
	
	#pulse_start 		= 97898;#97866;#96780;#95844;#94703;#0; # 0 if start from first pulse
	#pulse_start 		= 96866;#95948;#94743;#0; # 0 if start from first pulse
	pulse_start 		= 0; # 0 if start from first pulse
	
	#Frequency range under consideration
	MAX_F 				= 500e3; # Hz
	MIN_F 				= 0; # Hz
	f_split 			= 50e2; # Hz, to apply different thresholds above/below this frequency

	T_bin 				= 2.048e-3; # s
	df 					= 488; # Hz?
	
	ext_tgap 			= 0.04//T_bin;  #Time thresholds. How close must two indices be to be considered part of the same time structure.
	ext_fgap 			= 1500//df; 	#Frequency thresholds. How close must two indices be to be considered part of the same frequency structure. 

	time_filter			= True;
	SS_time_threshold 	= 75e-3; #s  #Consider modes only of length greater than 75ms
	
	size_filter 		= True;
	size_threshold 		= 100; #To remove noise. With the addition of the area filter, essentially not in use

	percentile_filter 	= True; 
	min_tperc 			= 0.05;  #Consider indices for each structure bounded by this percentile range
	max_tperc 			= 0.95;
	min_fperc 			= 0.05;
	max_fperc 			= 0.95;

	density_filter 		= True;
	density_threshold 	= 0.1; #To remove noise. This filter seems to have greatest effect on mode identification

	area_filter 		= True;
	area_threshold 		= 500; # units?

	ICRH_filter 		= True
	ICRH_f 				= 6; # units?
	
	merging 			= True; #Turn on merging function. If True, merges many identified modes that may be one fishbone mode/cut up by sawtooth crashes.
	
	#Would experiment with this on and off
	
	t_o 				= 0.1; # 
	f_o 				= 0.05; #Stitches together modes if overlap by 10% in time and 5% in f space. 

	TMN_threshold		= 0.6; # If toroidal mode number n comprises over X% of indices in mode, included (as well as percentage of indices that have this n).
	PMN_threshold 		= 0.25; # Same as above, except threshold to be included is Y%

	TS 					= [0,50e-3];# [0.5,0.25,0.1]; #s, times over which to average parameters of interest


if 1:#Auxiliary functions. Index density, structure area and plotting
	def index_density(sig_str):
		
		sig_str = np.array(sig_str)
		SS_time = sig_str[:,0];SS_freq = sig_str[:,1]
		t_0 	= np.amin(SS_time);t_100 = np.amax(SS_time)
		f_0 	= np.amin(SS_freq);f_100 = np.amax(SS_freq)
		t_10 	= np.quantile(SS_time,0.1);t_90 = np.quantile(SS_time,0.9)
		f_10 	= np.quantile(SS_freq,0.1);f_90 = np.quantile(SS_freq,0.9)
		f_t_area= (t_90-t_10)*(f_90-f_10) #Area in ft space
		ind_num = 0.81*len(sig_str)
		ind_dens= ind_num/f_t_area
		dens_hold.append((n,round(t_0,2),round(t_100,2),round(f_0,2),round(f_100,2),ind_dens))
		
		return ind_dens
		
	def structure_area(sig_str):
		sig_str = np.array(sig_str)
		SS_time = sig_str[:,0];SS_freq = sig_str[:,1]
		t_0 	= np.amin(SS_time);t_100 = np.amax(SS_time)
		f_0 	= np.amin(SS_freq);f_100 = np.amax(SS_freq)
		t_5 	= np.quantile(SS_time,0.05);t_95 = np.quantile(SS_time,0.95)
		f_5 	= np.quantile(SS_freq,0.05);f_95 = np.quantile(SS_freq,0.95)
		t_50 	= np.quantile(SS_time,0.50);f_50 = np.quantile(SS_freq,0.50)
		t_len 	= t_95-t_5; f_len = f_95-f_5
		structure_area = f_len*t_len
		
		return structure_area
		
	
	#PLOTTING FUNCTIONS:
	#call_line is line plot to see identified modes

	def call_line(mode=None):
		if mode==None:
			for nstruct in SIG_STRS_ALL:
				colour = np.array((np.random.choice(range(256), size=3)));colour = colour/256
				impl.line_plot(Shape,nstruct,min_f=MIN_F,max_f=MAX_F,colour=colour,tStart=tStart,tStop=tStop,pulse=pulse)
		else:
			nstruct=SIG_STRS_ALL[mode+1]
			colour = np.array((np.random.choice(range(256), size=3)));colour = colour/256
			impl.line_plot(Shape,nstruct,min_f=MIN_F,max_f=MAX_F,colour=colour,pulse=pulse)



if 1: #MAIN SETUP
	
	#The FILTERS AND FREE PARAMETERS section will also need to be checked
	
	#TEST: If true will break after k iteration, and return DataFrame of modes from these pulses
	test 	= True;
	k 		= 500;
	pulses 	= rDB.pulse_list()
	 

	#THE PATH TO DATBASE 1:
	folder 		= '/home/atinguel/work/fromIyngkarranK/';
	DB_path 	= folder+"DATABASE_MAIN.h5"

	#storing_file= folder+"test.pkl"  #Pkl file worked well for me when storing pd DataFrame
	#storing_file= folder+'MEA_211229.pkl'  #Pkl file worked well for me when storing pd DataFrame
	#storing_file= folder+'MEA_211230_'+str(int(pulse_start))+'.pkl'  #Pkl file worked well for me when storing pd DataFrame
	storing_file= folder+'MEA_'+yymmdd()+'_'+str(int(pulse_start))+'.pkl'  #Pkl file worked well for me when storing pd DataFrame
	logs_file 	= "log_file.txt"
	
	image_save 	= True; #if true, saves images of identified mode line plot to save folder
	image_save_freq = 2;
	save_folder = "/home/atinguel/work/fromIyngkarranK/figs/"

	#Parameters of interest. Can add in any that are in Database 1. Averaged over mode time (and ~10^2 ms prior)	
	PoIs 		= [  'NBLM',
					 'PTOT',
					 'TE0',
					 'angf_cxd6',
					 'angf_cxg6',
					 'angf_cxh4',
					 'angf_cxhm',
					 'angf_cxkm',
					 'btnm',
					 'bvac',
					 'elon',
					 'fD_kt5p',
					 'fH_ks3b',
					 'fH_kt5p',
					 'fHe3_kt5b',
					 'fHe4_kt5b',
					 'fT_kt5p',
					 'gradn95',
					 'ipla',
					 'ka2-1',
					 'ka2-2',
					 'ka2-3',
					 'ka2-4',
					 'ka2-5',
					 'ka3_ccds',
					 'n95',
					 'ne0',
					 'q0',
					 'q95',
					 's95',
					 'tril',
					 'triu'];
					
	cols 		= mode_columns+PoIs

	IDENTIFIED_MODES = pd.DataFrame(columns = cols) #Empty dataframe for mode characteristics
	
	
	if pulse_start==0:
		#Starting index of pulses
		pulse_counter=0
	else:
		#OR (comment out one not in use)
		try:
			pulse_counter = int(np.where(pulses==pulse_start)[0][0])
		except:
			raise Exception("Entered pulse not found in Database I")

num_pulses 		= 0; 

for pulse in pulses[pulse_counter:]: #MAIN LOOP
	
	times 		= rDB.spect_reconstruct(pulse)[5]
	tStart 		= times[0];tStop=times[-1]
	STORE 		= np.zeros([1000,60]) #Temporary storage of mode properties. There should not be 1000 identified modes (or more) per pulse.
	
	print("\n\n\n\n Extracting modes for JPN "+str(pulse)+", "+str(num_pulses))
	
	number_modes=0;
	
	if 1: #Getting data structures
		
		DATABASE 	= h5.File(DB_path,"r")
		sub_db_path = rDB.find_pulse(pulse)[0]
		sub_db 		= DATABASE[sub_db_path]
		
		pulse_data 	= sub_db[str(pulse)] 
		NZS 		= pulse_data["NZ indices"];NZS = deepcopy(np.array(NZS));NZS = np.transpose(NZS)
		TMN_vals 	= pulse_data["TMN spectogram"];TMN_vals = deepcopy(np.array(TMN_vals))
		Shape 		= pulse_data["Shape"];Shape = np.array(Shape)
		time 		= pulse_data["Time array (s)"];time = deepcopy(np.array(time))
		param_dict 	= {};
		for param in PoIs:
			try: 	param_dict[param]	= np.array(pulse_data[param])
			except: param_dict[param] 	= time*np.nan;
			
		DATABASE.close()
		
		SS_timebin_threshold	= SS_time_threshold//Shape[1]
		
		TMN,PMN,avgAMPL 		= rDB.spect_reconstruct(pulse)[1:4] #PMN and avgAMPL spects
		#PMN_vals = deepcopy(np.array(pulse_data["PMN spectogram"]))
		#avgAMPL = deepcopy(np.array(pulse_data["Averaged Amplitudes"]))
		
		tStart 		= time[0]; 
		tStop 		= time[-1]   #Start and end times of PULSE
		
	mode_dict 		= dict()  #Empty dictionary to store 
	SIG_STRS_ALL 	= []  #To store significant structures
	dens_hold 		= [] #Stores structure densities. For testing mainly.
	PMN_hold		= [] #Stores structure poloidal mode number arrays
	
	if 1: #Removal of indices with f > 500kHz. We have lower confidence in results outside this range (mostly ICRH/noise)
		#print(len(NZS))
		
		f_high_bin 	= MAX_F//df; 
		f_low_bin 	= MIN_F//df 
		freqs_all 	= NZS[:,0] #Currently indices as ft
		boolH 		= freqs_all >= f_high_bin; 
		boolL 		= freqs_all <= f_low_bin
		boolfreqs 	= np.logical_or(boolH,boolL)
		NZS 		= NZS[~boolfreqs]
		TMN_vals 	= TMN_vals[~boolfreqs]

	
	for n in NRange: #Filtering by TMN
		
		if 1: #Prepare indice array
			
			#print("\n\n\nFiltering for mode number "+str(n))
			booln 	= TMN_vals==n 
			NZS_n 	= NZS[booln] #NONZERO indices for some pulse with toroidal mode number n. First columns are ordered frequency bins
			
			#Mode separation
			NZS_n 	= np.flip(NZS_n,1) #tf, frequency ordered

		#NOTE. I HAVE FOUND F SEP FIRST, THEN T SEP TO BE MORE EFFECTIVE THAN T THEN F. 

		sig_strs 	= []
		for i in range(len(NZS_n)): #Freq separation of all indices of this mode number. sig_str is list of lists.
			
			if len(NZS_n) == 0: #No indices for n
				continue
			
			if i==0:
				sig_strs.append([])
				sig_strs[-1].append(NZS_n[i])
			else:
				if NZS_n[i][1]<f_split/df:
					if (NZS_n[i][1]-NZS_n[i-1][1])>ext_fgap:    #If difference greater than threshold 
							sig_strs.append([])    				#Add new list
							sig_strs[-1].append(NZS_n[i])    	#Append i to that new list
					else: 
						sig_strs[-1].append(NZS_n[i])
					
				
				if (NZS_n[i][1]-NZS_n[i-1][1])>ext_fgap:    	#If difference greater than threshold 
							sig_strs.append([])    				#Add new list
							sig_strs[-1].append(NZS_n[i])    	#Append i to that new list
				else: 
					sig_strs[-1].append(NZS_n[i])    #Add indice to existing list
					
		#print(len(sig_strs))
		ii=0
		sig_strs 	= [np.asarray(sig_str) for sig_str in sig_strs] #numpy arrays


		if 1: #Time separation
		
			n_ts 	= len(sig_strs) #Number of f sig strs
			for f_struct in sig_strs[:n_ts]:
				f_struct 	= np.array(f_struct)
				t_sort 		= np.argsort(f_struct[:,0]) #Order by times
				f_struct 	= f_struct[t_sort] #F structure now ordered by time. Still arranged as (t,f)
				for i in range(len(f_struct)):
					if i==0:
						sig_strs.append([])
						sig_strs[-1].append(f_struct[i])  #Create bin for first case		
					else: 
						if (f_struct[i][0]-f_struct[i-1][0])>ext_tgap:
							sig_strs.append([])
							sig_strs[-1].append(f_struct[i])
						else:
							sig_strs[-1].append(f_struct[i])
								
		sig_strs = sig_strs[n_ts+1:] #Once all frequency structures have been separated in time space, can remove original unseparated frequency structures (first n_ts)


		if 1: #FURTHER FILTERING. TIME, SIZE, PERCENTILE, DENSITY, AREA AND ICRH FILTERS
			pre 			= len(sig_strs)
			sig_strs 		= [np.asarray(sig_str) for sig_str in sig_strs]
			
			if time_filter==True:
				sig_strs 	= [sig_str for sig_str in sig_strs if np.subtract(np.amax(sig_str[:,0]),np.amin(sig_str[:,0]))>SS_timebin_threshold]
				
			if size_filter==True:
				sig_strs 	= [sig_str for sig_str in sig_strs if len(sig_str)>size_threshold]
			
			post 			= len(sig_strs)
			#print("From a list of "+str(pre)+" structures, upon applying time and size threshold filter, "+str(post)+" remain")
			#Now sig_strs[i] is ith FREQUENCY AND TIME separated structure
		

			if percentile_filter==True: #Calculation of min_perc and max_perc percentiles and removal of indices outside of this range

				i=0
				for sig_str in sig_strs: #So MaNy fOr LoOpS
					
					sig_str = np.array(sig_str)
					o_len 	= len(sig_str)
					SS_time = sig_str[:,0];SS_freq = sig_str[:,1]
					t_5 	= np.quantile(SS_time,min_tperc);t_95 = np.quantile(SS_time,max_tperc)
					f_5 	= np.quantile(SS_freq,min_fperc);f_95 = np.quantile(SS_freq,max_fperc)
					sig_str = [sig_str for sig_str in sig_str if sig_str[0]>t_5 and sig_str[0]<t_95]
					sig_str = [sig_str for sig_str in sig_str if sig_str[1]>f_5 and sig_str[1]<f_95]
					if len(sig_str)!=0:
						sig_strs[i]=sig_str
					else:
						sig_strs[i]="remove"
					pc_rem 	= ((o_len-len(sig_str))/o_len)*100
					#print("\n"+str(round(pc_rem))+" percent of indices removed from sig str by percentile crop")
					
					i+=1
						
				sig_strs 	= [sig_str for sig_str in sig_strs if sig_str!="remove"]		
				

			if density_filter==True:
				sig_strs_dens=[]
				i=0
				for sig_str in sig_strs:
					sig_str 	= np.array(sig_str)
					f_50_bin 	= f_split/df;#50000/df
					if np.mean(sig_str[:,1])<f_50_bin:
						#if index_density(sig_str)<0.05:	pass
						if index_density(sig_str) < density_threshold:	pass
						else:	sig_strs_dens.append(sig_str)
					else:
						#if index_density(sig_str)<0.1:	pass
						if index_density(sig_str) < density_threshold:	pass
						else:	sig_strs_dens.append(sig_str)
							
				sig_strs 	= sig_strs_dens
				

			if area_filter==True:
				sig_strs 	= [sig_str for sig_str in sig_strs if structure_area(sig_str)>area_threshold]
			
		
		#ICRH remove
			if ICRH_filter == True:
		
				sig_strs = [sig_str for sig_str in sig_strs if np.amax(np.array(sig_str)[:,1])-np.amin(np.array(sig_str)[:,1])>ICRH_f]
			
			
				sig_strs = [sig_str for sig_str in sig_strs if len(sig_str)!=0]
				
				
				
		SIG_STRS_ALL.append(sig_strs) #Add sig_strs of mode number n to main group of SIG_STRS_ALL
		number_modes+=len(sig_strs)

	
	
	#MERGING low freq structures
	#It was observed that low frequency structures were being cut up into ~10-20 modes by sawtooth crashes/fishbone modes?
	#This was an attempt to resolve that (which works decently)
	
	if merging==True:
		print("Pre merge "+str(number_modes))
		
		if 1: #Change SIG_STRS_ALL structure
			SIG_STRS_ONE 	= []
			for nstruct in SIG_STRS_ALL:
				for struct in nstruct:
					SIG_STR_ONE = SIG_STRS_ONE.append(np.array(struct))
			SIG_STRS_ALL 	= deepcopy(SIG_STRS_ONE)
			SIG_STRS_ONE 	= None
		
		SIG_STRS_M			= [];
		SIG_STRS_NM 		= []
		for STR in SIG_STRS_ALL:
			STR 			= np.array(STR)
			mean_f 			= np.mean(STR[:,1]);#print(mean_f)
			if mean_f>f_split/df:#50000/df:
				SIG_STRS_NM.append(STR)
			else:
				SIG_STRS_M.append(STR)
		#print(len(SIG_STRS_ALL),len(SIG_STRS_M),len(SIG_STRS_NM))
		
		if len(SIG_STRS_M)>1:
			SIG_STRS_M 	= merge.merge(SIG_STRS_M,t_o=t_o,f_o=f_o)
		else:			pass
		
		#SIG_STRS_ALL 	= SIG_STRS_M+SIG_STRS_N
		SIG_STRS_ALL 	= SIG_STRS_M+SIG_STRS_NM
		print("\n Post merge "+str(len(SIG_STRS_ALL)))
			

	#Save images
	if image_save==True:
		if num_pulses%image_save_freq==0:
			plt.clf()
			call_line()
			plt.savefig(save_folder+"/"+str(pulse)+"_identified modes")
			print("\n "+"SAVING IDENTIFIED MODES PLOT")
	
	if 1: #Finding defining properties (t,f,m,A,parameters) of each SS in sig _strs, then storing in datfram
		#print("For TMN "+str(n)+" there are "+str(len(sig_strs))+" structures identified to persist for "+str(SS_time_threshold)+" ms or longer")
		i=0
		for SS in SIG_STRS_ALL:
			
			if 1: #time, frequencies and amplitude
				time_in 		= np.array([inds[0] for inds in SS])  # time array of all indices in SS
				frequencies_in 	= np.array([inds[1] for inds in SS]) #S freq array of all indices in SS
						
				tb_min 			= np.min(time_in);
				tb_Med 			= np.quantile(time_in,0.5);
				tb_max 			= np.max(time_in)
				fb_min 			= np.min(frequencies_in);
				fb_Med 			= np.quantile(frequencies_in,0.5);
				fb_max 			= np.max(frequencies_in)
				t_find 			= lambda tb: tStart+tb*T_bin
				f_find 			= lambda fb: (fb*df)/1000     #kHz
				t_min 			= round(t_find(tb_min),4);
				t_Med 			= round(t_find(tb_Med),4);
				t_max			= round(t_find(tb_max),4);
				i_min			= np.argmin(np.abs(time_in-t_min)); # get index for time
				i_Med			= np.argmin(np.abs(time_in-t_Med));
				i_max			= np.argmin(np.abs(time_in-t_max));
				f_min 			= round(f_find(fb_min),2);
				f_Med 			= round(f_find(fb_Med),2);
				f_max 			= round(f_find(fb_max),2);
				#Finds min,med,max of time and freq values
							
				inds_n_SS 		= (frequencies_in,time_in) #tuple of freq time indices for this n,SS combination
				avgAMPL_n_SS 	= avgAMPL[inds_n_SS]
				A			 	= np.round(np.mean(avgAMPL_n_SS),8) #AVG AMPLITUDE of structure
				#A_min 			= np.min(avgAMPL_n_SS);
				#A_Med 			= np.quantile(avgAMPL_n_SS,0.5);
				#A_max 			= np.max(avgAMPL_n_SS);
				A_min 			= round(avgAMPL_n_SS[i_min],8); # get at appropriate time
				A_Med 			= round(avgAMPL_n_SS[i_Med],8);
				A_max 			= round(avgAMPL_n_SS[i_max],8);
				ind_max 		= np.argmax(avgAMPL_n_SS)
				t_A 			= t_find(inds_n_SS[1][ind_max])-tStart;
				f_A 			= f_find(inds_n_SS[0][ind_max])
				t				= [t_min,t_Med,t_max,round(tStart,2),t_A]
				f				= [f_min,f_Med,f_max,f_A]
				A 				= [A_min,A_Med,A_max,round(np.max(avgAMPL_n_SS),8)]
				
			if 1: #TMNS
				TMN_SS 			= TMN[inds_n_SS] #PMNS for this structure
				ns				= []

				for n in NRange:
					booln 		= TMN_SS==n
					TMN_n_SS 	= TMN_SS[booln] 
					nlen 		= len(TMN_n_SS)
					if (nlen/len(TMN_SS)) > TMN_threshold:#0.60:
						ratio 	= round(nlen/len(TMN_SS),2)
						ns.append((n,ratio))
						
				n = np.array(ns); n = np.ravel(n); n = list(n) #Mafs
				

			if 1: #PMNS
				PMN_n_SS 		= PMN[inds_n_SS] #PMNS for this structure
				ms 				= []

				for m in MRange:
					boolm 		= PMN_n_SS==m
					PMN_n_m_SS 	= PMN_n_SS[boolm]  
					mlen 		= len(PMN_n_m_SS)
					if (mlen/len(PMN_n_SS)) > PMN_threshold:#0.25:
						#If PMN m comprises greater than 25 pc of pmns in structure
						ratio 	= round(mlen/len(PMN_n_SS),2)
						ms.append((m,ratio))
						
				#m = tuple(ms) #Tuple of tuples
				#PMN_hold.append((n,tb_min,tb_max,fb_min,fb_max,m))
				m = np.array(ms); 
				m = np.ravel(m); 
				m = list(m) #Mafs
				
			if 1: #Parameters averaged over time interval
		
				#GETTING PARAMETERS:
				#NEED ONLY GET PARAMETERS OF INTEREST, PoIs given above
			
				diff 			= np.abs(time - t_Med)
				p_index 		= np.argmin(diff) #Index of minimum value of diff
				parameters_ofInt= []
				
				for param in PoIs:
					#TS = [0.5,0.25,0.1]
					param_avgs = []
					for i in TS:
						param_data 	= param_dict[param]
						boolL 		= time < (t_min-i); 
						boolR 		= time > t_max
						boolAll 	= np.logical_or(boolL,boolR)
						param_data 	= param_data[~boolAll] #Params in mode t interval
						param_data 	= param_data[~np.isnan(param_data)];
						param_data 	= param_data[~np.isinf(param_data)];
						param_avg 	= np.mean(param_data)
						param_avgs.append(param_avg)

					parameters_ofInt.append(param_avgs)
				
			
			if 0: #POPULATE STORE MATRIX
				
				vals 		= np.zeros(STORE.shape[1]) #Initialise storing array
				mode_props 	= [n,np.nan]+m+[np.nan]+t+[np.nan]+f+[np.nan]+[A,np.nan]+parameters_ofInt
				mode_props 	= np.asarray(mode_props)
				vals[:len(mode_props)] = mode_props
				STORE[number_modes] = vals #Store values in matrix
				number_modes+=1
			
			if 1:
				temp_dict		={}
				temp_dict["JPN"]=pulse
				temp_dict["n"]	=n
				temp_dict["m"]	=m
				temp_dict["T"]	=t
				temp_dict["F"]	=f
				temp_dict["A"]	=A
				temp_dict["TS"]	=TS
				i=0
				for param in PoIs:
					temp_dict[param]=parameters_ofInt[i]
					i+=1
				
				#print("F range: "+str(f))	
				IDENTIFIED_MODES = IDENTIFIED_MODES.append(temp_dict,ignore_index=True)
				number_modes+=1

				temp_dict.clear()
					
					
	STORE= np.delete(STORE,np.where(~STORE.any(axis=1))[0], axis=0)	
	print(str(len(SIG_STRS_ALL))+" modes identified for pulse "+str(pulse))
		
		
	if 1: #SAVE AND UPDATE LOGS
		with open(logs_file,"a") as logs:
			logs.write("\n\n\n\n Saved individual mode data from "+str(pulse))
			logs.write("\n "+str(number_modes)+ " structures identified")
		

	
	
	
	
	if test==True: #TEST CASE
		if num_pulses==k:
			call_line();
			#plt.show()
			break
	
	num_pulses+=1
	pulse_counter+=1
	

	
print("FINISHED")
print("\n SAVING TO PICKLE FILE")
IDENTIFIED_MODES.to_pickle(storing_file)
		
		
