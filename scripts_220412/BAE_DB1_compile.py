"""

##Please see help in UAE_DB interfacing script for more info.

The following scripts, and the accompanying ones provided in the folder Auxscripts have been created for the
compilation of a DATABASE I type file.

Copied from UAE_DB1_compile.py

AT 210908: Try to compile a database of BAEs for L Piron

"""


if 1: ##IMPORT
	import sys
	sys.path[:0]=['/jet/share/lib/python'];
	from ppf import *
	import getdat
	import numpy as np
	import scipy
	from scipy import io as sio
	from scipy import interpolate as interp
	from scipy.sparse import csr_matrix
	sys.path[:0]=['/home/pgpuglia/AEAD/GUI/class'];
	from INCAA_functions_all import INCAA_Iant_Vant
	import matplotlib.pyplot as plt
	from copy import deepcopy
	import pandas as pd
	from pandas import HDFStore
	import h5py as h5
	import glob
	#import pdb #debugging
	#pdb.set_trace() #debugging

	from Auxscripts import preprocessing as prep
	from Auxscripts import plotting as impl
	from Auxscripts import new_mode_number_calculation as nmnc
	from Auxscripts import get_ICRH_signal as get_ICRH
	from Auxscripts import get_parameters as get_params


#-----------------------NODE AND PROBE SETTINGS-----------------------------------
if 1:
	node 					= 'DA/C1M-';    #Node from which to collect data from

	#test_pr 				= False; 
	#k						= 8;

	#Toroidal probes and angles (radians)
	toroidal_probes_main 	= ['T001', 'T002', 'T006', 'T007', 'T008', 'T009',
							   'H301', 'H302', 'H303', 'H304', 'H305', 'H306', 'H307'];#, 'H308'];
	# NOTE: only 11 probes used for TMN calculation, not H306,7,8 at same location as H305
	phi_t 					= [3.0,42.2,182.9,222.2,257.1,290.4,77.0,92.9,103.1,108.4,110.4]; 
	##toroidal_probes_main 	= ['H301','H302','H303','H304','H305',
	##						   'T001','T002','T006','T007','T008','T009'];
	##phi_t 					= [77.0,92.9,103.1,108.4,110.4,
	##						   3.0,42.2,182.9,222.2,257.1,290.4]; # degrees
	phi_t_main 				= np.radians(phi_t) #radians
	n_torprobes_main 		= len(toroidal_probes_main)

	#Poloidal probes and location (R,Z).
	poloidal_probes_main 	= ["PP801","PP802","PP803","PP804","PP805","PP806","PP807"]
	R_p_main 				= np.array([3.766,3.843,3.952,3.983,3.863,3.673,3.450]) #   Radial location of probes. 
	Z_p_main 				= np.array([1.249,1.075,0.694,0.184,-0.383,-0.770,-1.057]) #Vertical axis location of probes
	m_var_power 			= 2; # power of the distance between probe and plasma centroid to use in variance
	
	# JET parameters
	R0 						= 3; # m, normal major radius
	Z0 						= 0.2; # m, normal minor radius
	
	#Angle measurement uncertainty
	#thetaUNC 				= np.radians(30)
	n_polprobes_main 		= len(poloidal_probes_main)
	probes_main 			= toroidal_probes_main+poloidal_probes_main

##------------------------------------------------------------------------


#----------------PLASMA AND OPERATIONAL PARAMETERS----------------------
if 1:
	#To add param, follow structure below
	params = [	['magn','bvac','jetppf'  ], 								# Vacuum toroidal field at R=2.96 m (T)
				['magn','ipla','jetppf'  ], 								# Plasma current (A)
				['hrtx','ne0' ,'jetppf'  ], 								# Thomson scattering, electron density on axis, m-3 - chosen by Alex
				['HRTX','TE0' ,'jetppf'  ], 								# Thomson scattering, electron temperature on axis, eV?
				['efit','q95' ,'jetppf'  ], 								# q at 95% flux surface - chosen by Alex
				['efit','s95' ,'jetppf'  ], 								# shear (r/q)*(dq/dr) at 95% flux surface
				['efit','elon','jetppf'  ], 								# elongation - chosen by Alex
				['NBI' ,'NBLM','jetppf'  ], 								# Neutral beam power, W
				['ICRH','PTOT','jetppf'  ], 								# The total ICRH power coupled to the plasma, W
				['efit','triu','jetppf'  ], 								# upper triangularity
				['efit','tril','jetppf'  ], 								# lower triangularity
				['efit','q'   ,'jetppf'  ,'q0'],							# q profile, but just get q0
				['efit','btnm','jetppf'  ], 								# beta normalized (1)
				['ks3b','hthd','jetppf'  ,'fH_ks3b'], 						# this is H/(H+D+T) so ~ n_H/n_e
				['kt5p','hthd','jetppf'  ,'fH_kt5p'], 						# this is H/(H+D+T) so ~ n_H/n_e
				['kt5p','dthd','jetppf'  ,'fD_kt5p'], 						# this is D/(H+D+T) so ~ n_D/n_e
				['kt5p','tttd','jetppf'  ,'fT_kt5p'], 						# this is T/(H+D+T) so ~ n_T/n_e
				['kt5b','che' ,'edelabie','fHe3_kt5b'], 					# this is n_He/n_e, so max is 0.5
				['kt5b','che' ,'edelabie','fHe4_kt5b'], 					# this is n_He/n_e, so max is 0.5
				['hrts','n95' ,'jetppf'  ], 								# Thomson scattering, electron density at q95, m-3 - chosen by Alex
				['hrts','gradn95','jetppf'], 								# Thomson scattering, gradient of n95
				['cxd6','angf','jetppf','angf_cxd6'], 						# Rotation angular frequency (rad/s)
				['cxg6','angf','jetppf','angf_cxg6'], 						# Rotation angular frequency (rad/s)
				['cxh4','angf','jetppf','angf_cxh4'], 						# Rotation angular frequency (rad/s)
				['cxhm','angf','jetppf','angf_cxhm'], 						# Rotation angular frequency (rad/s)
				['cxkm','angf','jetppf','angf_cxkm'], 						# Rotation angular frequency (rad/s)
				];
	nParams 		= len(params);
	param_cols 		= [p[1] for p in params]

	parameters_df 	= pd.DataFrame(columns=param_cols)


#----------------SCRIPT FUNCTIONALITIES & FREE PARAMETERS------------------
if 1:

	#tStart = float(input("Enter start time: "))
	#tStop = float(input("Enter end time: "))
	#pulses = list(input("Enter Pulse list: "))


	#---------------Generate list of spectograms and averaged Coherence spectogram. Required true for all other functionalities
	spctAC = True;

	threshold 		= 0.6     	#avgCOH coherence threshold
	take_freq 		= 4     	#In list of NC2 pairings of N probe spectrograms, how often to choose pairing to include in avgCOH calc. Low take_freq --> stronger Coherence
	coh_n 			= 5   		#Time bin overlap for coherence matrix creation
	T_bin 			= 2.048e-3 	#s, time bin for spectrogram
	df 				= 1//T_bin 	#Hz
	min_fres 		= 1/T_bin 	#Minimum resolvable frequency from Fourier Transform, numerator is # cycles
	#-----------------------------

	#------------Antenna signal remove 
	antenna_remove	= True
	Athreshold 		= 0.25; 	# If amplitudes averaged over bin is greater than this threshold, consider AE and do not remove
	fb_width0 		= 1500  #Frequency width of antenna signal (Hz)
	tb_width0 		= 0.01  #Time bin width of antenna signal (s)
	Ithreshold 		= 2;   	#Antenna current threshold (we only consider input currents of 2A or greater) (A)
	fpres 			= 1     #Used to retain signal that intersects with TAEs
	tpres 			= 25	#Used to retain signal that intersects with TAEs
	af_threshold 	= threshold 	#Removes bins that have not been completely cut out by removing reconstructed antenna signal
	#------------------

	#-----------ICRH signal remove
	ICRH_remove		= False;#True
	#---------------------

	#------------Toroidal mode number spectogram
	tmn				= True
	NMax 			= 8; 
	NMin 			= 0;#-1;
	#-----------------------

	#--------------Poloidal mode number
	pmn				= True;
	fixed_centroid 	= False;
	MMaxFactor 		= 4; # Let MMax be this factor times the NMax
	MMax 			= MMaxFactor*NMax;
	MMin			= NMin
	#-------------------

	#----------------Averaged Amplitude spectgram
	amplitude_spect	= True
	#------------------------
	
	#----------------Saving spectrograms
	image_save		= True
	#------------------------

	#----------------Acquire plasma parameters
	get_parameters 	= True
	#------------------------

	#-------------------Mode Properties extraced to Database
	#extract is general extraction (can extract individual modes and/or spectograms),
	#storing is storing data in hdf5 file, periodic save - Save every save_freq pulses.
	#Following recent edits, would personally NOT CHANGE three below
	extract			= True; 
	storing			= True;#False; 
	periodic_save	= True;

	#Testing
	test_pulse		= False  	#CHOOSE ONE PULSE
	test_pr 		= False  	#CHOOSE RANGE OF PROBES
	k 				= 2 		#No of test probes

	#-------------------------
	#If feeling lazy
	ALL 			= False;
	#if ALL==True:	spect_AC=True;antenna_remove=True;ICRH_remove=True;tmn=True;pmn=True;amplitude_spect=True;get_parameters=True;image_save=True;extract=True;storing=True;periodic_save=True
	if ALL==True:	spect_AC=True;antenna_remove=True;ICRH_remove=False;tmn=True;pmn=True;amplitude_spect=True;get_parameters=True;image_save=True;extract=True;storing=True;periodic_save=True

	#To skip code. This is used when most probes do not return signal. KEEP AS FALSE
	SKIP 			= False


##--------------------------MAIN SETUP----------------------------------------------

if 1: #Pulse setup prior to main loop. Option for one test pulse here
	
	"""
	IMPORTANT:
		#To alter the fast magnetics probe config, see -NODE AND PROBE SETTINGS
		#To alter the parameter time-series incorporated into the DB, see -PLASMA AND OPERATIONAL PARAMETERS
		#To alter the main script functionalities and associated free parameters, see -SCRIPT FUNCTIONALITIES & FREE PARAMETERS
	"""
	
	#Saving settings
	#storing_file 	= "DATABASE_testingBAE.h5" 							#HD5 file path
	storing_file 	= "DATABASE_BAE.h5" 								#HD5 file path
	save_freq 		= 2; 												#Save to hdf every `save_freq` pulses
	#logs_file 		= "./pulses_completed_testing_BAE.txt" 				#Logs path
	logs_file 		= "./pulses_completed_BAE.txt" 						#Logs path

	#Image folder #Save every image_save_freq pulse. Choose frequency range for PLOTS ONLY
	image_folder 	= "./figs/"
	image_save_freq = 1;
	im_minf			= 0; # kHz
	im_maxf			= 500; # kHz 
	
	
	#Pulses to analyse here. 
	#pulse_data 		= sio.loadmat("/home/atinguel/work/forIyngkarranK/pulses_w_ICRH_93000-98006_Pthreshold500kW_Tthreshold500ms")
	#pulses 			= pulse_data["pulses"][0]
	#tStarts 		= pulse_data["tStart"][0]
	#tStops 			= pulse_data["tEnd"][0]
	
	# from L Piron email, September 2021
	pulses0 		= np.array([95867,95871]);
	pulses1			= np.array([95601,95602,95603,95604,95865,95867,95869,95870,95871,95872,95873,95874,95875,
								95877,95878,95879]);
	pulses2			= np.array([98066,98067,98069,98071,98088,98089,98090,98091,98094,98095,92018,92019,90594,
								90595,90596,90598,90599,90600,90601,90602,90603,90611,90612,90613,90614,90615,
								90616,90617,90618,90619,90620,90621,90624,90625,90675,90676,90679,90680,90681,
								90684,90685,90686,90710,85974,85975,85976,85977,91025,91026,91028,91029,91030,
								91031,91032,91033,91036,91037,91038,91039,91040]);
	
	pulses 			= pulses1; 
	#pulses 			= np.array([95867]);
	#pulses 			= np.concatenate((pulses0,pulses1,pulses2));
	#pulses 			= np.concatenate((pulses1,pulses2));
	pulses 			= np.unique(pulses);
	
	pulse0 			= 95879;#95872;#95869;#95865;#95603; # start from this pulse
	pindx 			= np.argmin(abs(pulses-pulse0));
	pulses 			= pulses[pindx:];
	
	tStarts 		= np.zeros(len(pulses))+40; # s
	tStops 			= np.zeros(len(pulses))+60; # s
	
	#tStarts 		= np.zeros(len(pulses))+45; # s
	#tStops 			= np.zeros(len(pulses))+45.1; # s

	#CHOOSE
	pulse_counter 	= 0  #Set to k-1 to start at kth pulse
	#OR
	#pulse_counter 	= np.where(pulses==94700)[0][0] 
	
	test_pulse 		= False
	if test_pulse	==True:
		pulse 		= 96851
		pulse_counter= np.where(pulses==pulse)[0][0] #1 specific pulse
	test_timings 	= False #True if want to compile for test_time_length interval (for debugging mainly)
	test_time_length= 0.1
	
	
if storing: ### Check if files exist ###
	if len(glob.glob(storing_file))==0: # if no storing file
		hFile 		= h5.File(storing_file,'w');
		hFile.close();
	if len(glob.glob(logs_file))==0: # if no logs file
		lFile 		= open(logs_file,'w');
		lFile.close();
	

#MAIN LOOP##############################################################
	
num_pulses				= 1 	#Number pulses iterated through
pulses_complete 		= [] 	#List of completed pulses

for pulse in pulses[pulse_counter:]:
#for pulse in pulses[pulse_counter:pulse_counter+2]:			
	#Probes are redefined from last run (some were removed if no signal returned)
	if 1:
			
		toroidal_probes = deepcopy(toroidal_probes_main)
		phi_t 			= deepcopy(phi_t_main)
		ntorprobes 		= len(toroidal_probes)

		poloidal_probes = deepcopy(poloidal_probes_main)
		R_p 			= deepcopy(R_p_main)
		Z_p 			= deepcopy(Z_p_main)
		npolprobes 		= len(poloidal_probes)
		
		probes 			= toroidal_probes+poloidal_probes
				
	print("\n\n\n\n\n\n\n\n\nPulse "+str(pulse_counter)+" of "+str(len(pulses))+" pulses.")
	print("\n\nCURRENT JPN: "+str(pulse))
	
	#Defining pulse dependent parameters (tStart, tSTop etc)
	if 1:
		tStart 			= tStarts[pulse_counter]
		tStop 			= tStops[pulse_counter]
		
		if test_timings==True:
			#TEST TIMES
			tStop_orig 	= tStop
			tStop 		= tStart+test_time_length
	

#------------------Spectrogram and avgCOH creation---------------------------

	if spctAC==True:
		
		i=0
		
		if test_pr==True:	probes 		= probes[0:k] #For quick testing
		else:				pass
			
		#These are to find common range in which ALL probes record data
		t_hold			= [0,0]
		probe_tstarts	= []
		probe_tstops	= []
		
		#Spectogram list
		spectograms		= []
		
		#i=0	
		#Processing data from one probe at a time		
		print("\n\nGETTING SIGNALS AND CREATING SPECTOGRAMS FOR "+str(len(probes))+" PROBES FOR TIME INTERVAL: "+str(round(tStart,2))+" - "+str(round(tStop,2)))
		
		sampling_freqs 		= []
		for probe in probes:  #Loops through ALL probes for spectogram list
		
			signalName 		= node+"-"+probe
			times 			= "TSTART ="+ str(round(tStart,2))+ " TEND = "+str(round(tStop,2))
			signalName_time = node+probe+'?TSTART='+str(round(tStart,2))+'&TEND='+str(round(tStop,2))
		
			signal,time_r,_,__,___,ier 	= getdat.getdat8(signalName_time,pulse)    #Signal from probe,JPN during given time interval.
			
			
			if len(signal)==0:
				print("\n NOTE: N0 signal returned for "+str(signalName)+" during time interval "+times)
				if i<len(toroidal_probes):
					ntorprobes-=1
					phi_t.pop(i); # edit from I Kumar 210913
				if i >=(len(probes)-len(poloidal_probes)):
					n_polprobes-=1
					R_p1.pop(i-len(toroidal_probes))
					Z_p1.pop(i-len(toroidal_probes)) #Remove coordinates from list
					
				continue
				
			dt 			= np.mean(np.diff(np.array(time_r)));#[0:100]))
			fs 			= np.round(1/dt,-6)
			sampling_freqs.append(int(fs))
			#print("\n Probe Sampling frequency is "+str(fs))
				
			probe_tstarts.append(time_r[0]);
			probe_tstops.append(time_r[-1])
			
			if fs!=2e6:
				#Fit to 2MHZ sampling freq. We find all spectrograms assuming fs = 2MHz. Then, pulses with probes at 1MHz are fitted to 0-500kHz range later in script (perhaps better way to do this)
				dt 		= 5e-7
				time_r_n= np.arange(time_r[0],time_r[-1],dt)
				signal 	= np.interp(time_r_n,time_r,signal)
								
			spect,freqs = prep.spectogram(signal,start=tStart,stop=tStop,Fs=2e6)[0:2]
			spectograms.append(spect)
						
			#Update t_hold if need be
			if i==0:
				t_hold[0]	=time_r[0]
				t_hold[-1]	=time_r[-1]
			else:
				fst = time_r[0];
				lst	= time_r[-1]
				probe_tstarts.append(fst);
				probe_tstops.append(lst)
				
				if fst>t_hold[0]:	t_hold[0]=fst
				if lst<t_hold[-1]:	t_hold[-1]=lst
					
			i+=1		
								
		
		if len(spectograms)<=15:	
			print("\n\n Insufficient spectogram number")						
			SKIP 	= True
			if SKIP == True:	#SKIP THE REST OF THE ALGORITHM				
				#`Turn off' the rest of the algorithm
				spect_AC=False;antenna_remove=False;ICRH_remove=False;tmn=False;pmn=False;amplitude_spect=False;get_parameters=False;image_save=False;extract=False;storing=False;periodic_save=False
				print("\n\n SKIPPING REST OF ALGORITHM")
			
		
		if SKIP==False:
			t_Start = t_hold[0]
			t_Stop 	= t_hold[-1]
		#These are times in which ALL probes record signal			


			#NOW WE MUST ENSURE THAT ALL SPECTOGRAMS ARE SAME SHAPE
			if 1:	
				ii=0
				tStart 			= t_hold[0];
				tStop 			= t_hold[-1] #Overlapping time interval
				print("\nTime interval in which all probes record signal is "+str(round(t_hold[0],2))+" - "+str(round(t_hold[-1],2)))
				print("\nFitting spectograms to this range")
				for spect in spectograms:
					spect_tl 	= (spect.shape[1])
					spect_tint 	= probe_tstops[ii]-probe_tstarts[ii]
					spect_tres 	= spect_tint//spect_tl #Time interval covered per bin in spectogram
					stt_bin 	= ((tStart-probe_tstarts[ii])/spect_tint)*spect_tl
					stt_bin 	= math.floor(stt_bin)
					stp_bin 	= ((tStop-probe_tstarts[ii])/spect_tint)*spect_tl
					stp_bin 	= math.floor(stp_bin)
					spect 		= spect[:,stt_bin:stp_bin+1]  #Cut spectogram to range
					spectograms[ii]=spect
					ii+=1
										
				#For small discrepancies
				#jj=0
				tbases 	= [spect.shape[1] for spect in spectograms]
				t_l 	= min(tbases)
				jj=0
				for spect in spectograms:
					if spect.shape[1]==t_l:
						jj+=1
						continue
					else:
						aa = t_l-spect.shape[1]
						aa = np.abs(aa)
						if aa<=2:
							print("Small spectogram shape error. Resolving...")
							spectograms[jj]=spect[:,:-aa] #If small error, truncate					
						else:
							raise ValueError("Spectogram shape error")
						jj+=1
				
					 
			##We take coherence of all probes with first, then divide by number of probes for average coherence matrix

			spectograms = np.array(spectograms)
			print("\n\n\nCreating averaged Coherence Matrix from "+str(len(spectograms))+" probes.")

			#avgCOH 	= prep.avgdCoherence(spectograms,threshold=0.6,n=coh_n,tk=int(take_freq))
			avgCOH 		= prep.avgdCoherence(spectograms,threshold=threshold,n=coh_n,tk=int(take_freq))

			#Finally truncate time to t_hold
			time 		= np.linspace(t_hold[0],t_hold[-1],avgCOH.shape[1])
				
			#We now have Coherence matrix for all probes in pulse, as well as list of spectograms.
			#avgCOH_nt = deepcopy(avgCOH)
			#avgCOH[avgCOH<threshold] = 0
		
			print("\n Average coherence matrix with threshold found")
			#print("\n"+str(sampling_freqs))
			
			t_bin 		= ((time[-1]-time[0])/(avgCOH.shape[1])) #Effective T_bin. 


	#NOTE THAT SPECTOGRAMS IN FORM (probe,frequencies,time)
	##-----------------------------------------------------------------------


	#---------------------Antenna frequency and ICRH remove---------------------------
	if antenna_remove==True:
		
		fb_width 		= int(fb_width0//df)   	#Frequency width of reconstructed antenna signal
		tb_width 		= int(tb_width0//T_bin)  #Time bin width of reconstructed antenna signal
		
		if antenna_remove == spctAC == True:	pass
		else:	raise Exception("Antenna removal requires spectogram calculation first.")
		
		print("\n\n\nREMOVING ANTENNA SIGNAL"); ant_rem_proceed=True
		
		sys.path[:0]=['/home/pgpuglia/AEAD/GUI/class'];
		from INCAA_functions_all import INCAA_Iant_Vant
		try: 
			time_MD,frequency_MD,Iant_MD,Vant_MD,t,f,Iant,Vant = INCAA_Iant_Vant(pulse);
			ant_rem_proceed=True		
		except:
			print("\n No antenna signal returned for JPN: "+str(pulse))
			ant_rem_proceed=False

		if ant_rem_proceed==True:
			
			Iant2 		= np.sum(np.abs(Iant),axis=1); 
			boolI 		= Iant2 >= Ithreshold; #Indices greater than threshold
			f 			= f[boolI]
			t 			= t[boolI];#print(t[0]);print(t[-1])
			ant_rem_proceed=True
			
			#Even if antenna signal present, may not be in same interval as ICRH range:
			
			boolStart 	= t < tStart
			boolStop 	= t > tStop
			boolAll 	= np.logical_or(boolStart,boolStop) #Boolean, with Trues being t elements that are outside range tStart-tStop
			t 			= t[~boolAll] #t relevant to this time window
			f 			= f[~boolAll]
			
			if len(t)==0:
				print("\n No antenna signal for "+str(tStart)+" - "+str(tStop)+" range.")
				ant_rem_proceed = False
			
			
		if ant_rem_proceed==True:
			t0 			= t[0]
			tn1 		= t[-1]
			
			#Interpolate t to avgCOH timebase
			tlen 		= tn1-t0; nt_bins = int(tlen//t_bin)  #Number of bins in avgCOH resolution
			t_interp 	= np.linspace(t0,tn1,nt_bins)
			f 			= np.interp(t_interp,t,f)
			

			#Interpolation to fit number of time bins in interval

			"""
			##Restricting to times of interest
			start_bin = int(round(((tStart-t0)/(tn1-t0))*(t.shape[0])))
			stop_bin = int(round(((tStop-t0)/(tn1-t0))*(t.shape[0])))
			if (start_bin<0) or (stop_bin<0):
				raise Exception("Invalid time interval entered for antenna frequency. Strong antenna signal ranges from "+str(t0)+" - "+str(tn1))
			

			t=t[start_bin:stop_bin] #Restrict f and t to bins off interest
			f=f[start_bin:stop_bin]

			interp_t = np.linspace(tStart,tStop,avgCOH.shape[1])    #Interpolate to give same shape as avgCOH
			interp_f = np.interp(interp_t,t,f)
			"""

			#R = np.zeros([freqs.shape[0],interp_t.shape[0]])  #Array to hold antenna spectogram for all times
			#Create antenna signal spectogram		
			R 			= np.zeros(avgCOH.shape)
			start_bin 	= int((t0-tStart)//t_bin);
			stop_bin	= int((tn1-tStart)//t_bin)
			
			j=1
			r_len 		= R.shape[1]
			for j in range(1,4): #j is harmonic multiplier, going up to 3
				for i in range(start_bin,stop_bin):

					#Assign antenna freq for each time bin
					try:	frequency = (f[i-start_bin])*1e3  #As f returned in kHz. 
					except: frequency = (f[i-1-start_bin])*1e3

					n_f 	= round(frequency/df)*df  #Round to closest frequency in freqs array
					f_bin 	= round((n_f/1e6)*avgCOH.shape[0]) #Find corresponding bin number
					f_bin 	= j*f_bin   #Harmonic multiplier
								
					#Check the surrounding time bins in avgCOH. As TAEs present for a significantly longer time than ant_f, use this to differentiate. 
					#Look at surrounding points in avgCOH; if all ~0.9, indicates a TAE, so do not add to remove spectogram
					
					A = avgCOH[int(f_bin-fpres):int(f_bin+fpres),int(i-tpres):int(i+tpres)]
					#if np.sum(A) >= (0.5*(A.size))*0.5:    continue;	#If given fraction of surrounding structure is signal, then most likely AE. Leave this region of R as zero (remove nothing)
					if np.sum(A) >= A.size*Athreshold:	continue;	#If given fraction of surrounding structure is signal, then most likely AE. Leave this region of R as zero (remove nothing)
						
					try:	R[int(f_bin-fb_width):int(f_bin+fb_width),int(i-tb_width):int(i+tb_width)]=1   #Sets bin and surrounding points to scale
					except:	R[int(f_bin),int(i)]=1   #Sets bin and surrounding ones to 1
					
			R_F 	= np.multiply(avgCOH,R) #Scale antenna frequency
			avgCOH 	= avgCOH-R_F
			avgCOH[avgCOH<af_threshold]=0  #Negative values set to 0 too, which we should obtain when no harmonics
			print("\n Antenna signal removed")
			
		if ant_rem_proceed==False:
			np.zeros(avgCOH.shape)


	if ICRH_remove 	== True:
		if spctAC 	== False:	raise Exception("ERROR. Averaged Coherence calculation required for ICRH removal")
		print("\n\n\nREMOVING ICRH SIGNALS")
		
		#Pull ICRH signals
		FB,FA,TB,TA 			= get_ICRH.get_ICRH(pulse,tStart,tStop,avgCOH.shape[1])  #Frequency signals for avgCOH, interpolated for timebase already
		FB 						= np.array(FB);FA = np.array(FA)
		ICRH_signals			= np.vstack((FB,FA))
		ICRH_signals_bins		= np.rint(ICRH_signals/df)
		R_i 					= np.zeros([2*avgCOH.shape[0],avgCOH.shape[1]]) #Initialize ICRH matrix
		
		for iC in range(avgCOH.shape[1]):
			bin_col 			= ICRH_signals_bins[:,iC]    
			bin_col[bin_col>avgCOH.shape[0]] = 0  #Ignore if out of frequency range
			bin_col 			= bin_col.astype(int)
			R_i[:,iC][bin_col]	= 1; 
			R_i[:,iC][bin_col-1]= 1;  
	
		fbin_nyq 	= int(2*avgCOH.shape[0]//2)
		R_i 		= R_i[:fbin_nyq,:]    #Only want frequency bins below F_nyquist
		R_I 		= np.multiply(R_i,avgCOH) #scaling
		avgCOH 		= avgCOH-R_I #removal
		avgCOH[avgCOH<threshold] = 0
		print(" ICRH signals removed")
		
		avgCOH_sp = csr_matrix(avgCOH)
			
	#-----------------------------------------------------------------------------

			
	#---------------------------TOROIDAL MODE NUMBER SPECTROGRAM--------------------	
	if tmn==True:
		
		rav 			= np.ravel(avgCOH); 
		nzs 			= np.nonzero(rav)[0];
		
		if tmn == spctAC == True:	pass
		else:	raise Exception("Toroidal mode number calculation requires individual spectogram calculation.")
		
		print("\n\n\nTOROIDAL MODE NUMBER CALCULATION")
		
		shape 			= avgCOH.shape
		avgCOH 			= np.ravel(avgCOH)
		indices 		= np.nonzero(avgCOH)[0] 
		avgCOH 			= avgCOH.reshape(shape)
		indices 		= indices.astype(int)    #Non zero elements of avgCOH in 1D array
		
		#Consider only nonzero indices
		#Nonzero indices list
		i=0
		#print(len(indices),len(phi))
		PHASE 			= np.zeros([len(indices),len(phi_t)])   #Initialize phase matrix. 
		for spect in spectograms[:len(phi_t)]:
			ravel 		= np.ravel(spect)
			nonzero 	= ravel[indices]   #Fourier coefficients of interest
			angle 		= np.angle(nonzero)
			PHASE[:,i]	= angle
			i+=1
		
		#NEST= mnc.calculate_nmode(avgCOH,"toroidal",phi_t,PHASE,NMax,Ip=1)[2]
		PHASE 			= PHASE[:,np.newaxis,:]
		NEst 			= nmnc.calculate_mode_num(phi_t,PHASE,nzs,NMax,NMinAllowed=NMin,Ip=1)[2] #TMN calculation lx1
		
		#Broadcasting back to avgCOH shape
		NEST 			= np.zeros([avgCOH.size,1]);
		NEST[:]			= np.nan
		NEST[nzs]		= NEst
		NEST 			= np.reshape(NEST,avgCOH.shape)
		print("\n Toroidal mode number calculation complete")
		
		NEST_sp 		= csr_matrix(NEST)

	else:
		NEST 			= None
		NEST_sp 		= None	
			
	##-----------------------------------------------------------------------------
	
	
	#---------------------------POLOIDAL MODE NUMBERS-new----------------------------
	if pmn==True:
		if test_pr==True:	raise Exception("ERROR. Poloidal number calculation will not work in test probe case")
		
		#For fixed time bin
		
		if pmn == spctAC == True:	pass
		else:	raise Exception("Poloidal mode number calculation requires individual spectogram calculation. I'm not angry... just dissapointed.")
			
		if fixed_centroid==True:
			#THIS PART OF CODE NOT UP TO DATE
			#Same as TMN case
			Rc 				= R0;#3; # m
			Zc 				= Z0;#0.2; # m
		
			print("\n\n POLOIDAL MODE NUMBER CALCULATION")
			shape 			= avgCOH.shape
			avgCOH 			= np.ravel(avgCOH)
			indices 		= np.nonzero(avgCOH)[0]
			avgCOH 			= avgCOH.reshape(shape)
			indices 		= indices.astype(int)
			
			phi_p 			= np.arctan2(Z_p-Zc,R_p-Rc)   #Assuming plasma centroid located at centre of tokamak and unchanging in time
			PHASE 			= np.zeros([len(indices),len(R_p)])
			i=0
			for spect in spectograms[-(len(poloidal_probes)):]:
				ravel 		= np.ravel(spect)
				nonzero 	= ravel[indices]   #Fourier coefficients of interest
				angle 		= np.angle(nonzero)  #Poloidal phase angles
				PHASE[:,i]	= angle
				i+=1
				
			dZ2 			= (Z_p-Zc)**2; 
			dR2 			= (R_p-Rc)**2; 
			var 			= np.power(dZ2+dR2,m_var_power)
			MEst			= nmnc.calculate_nmode(phi_p,PHASE,MMax,NMinAllowed=MMin,Ip=1,var=var)[2]    #axis="toroidal" as probes at fixed location in space over time, similar to toroidal probes
			
			
		else:    #A time varying centroid is assumed
			
			print("\n\n\nPOLOIDAL MODE NUMBER CALCULATION")
			
			#Begin by pulling centroid location as function of time, Rc and Zc
			##pulse 	= pulse
			Rc,x,tc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
			ppfdata(pulse,'EFIT','RC',seq=0,uid='jetppf',device='JET',\
			fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);
			
			# Get vertical position of current centroid
			Zc,x,tc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
			ppfdata(pulse,'EFIT','ZC',seq=0,uid='jetppf',device='JET',\
			fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);
			
			try:
				Zc 			= np.interp(time,tc,Zc)
				Rc 			= np.interp(time,tc,Rc) #Interpolate for time resolution
			except:
				Zc 			= np.zeros(time.size)+Z0;
				Rc 			= np.zeros(time.size)+R0;

			
			Zpos 			= np.meshgrid(Z_p,np.linspace(0,1,Zc.shape[0]),indexing="ij")[0]
			Rpos 			= np.meshgrid(R_p,np.linspace(0,1,Rc.shape[0]),indexing="ij")[0]
			phi0 			= np.arctan2(Zpos-Zc,Rpos-Rc)    #Poloidal probe angles vs time
			phi0 			= np.transpose(phi0)  #time vs probes. ijth is angle of ith probe in jth time bin
			
			#var_t 			= (np.power(Zpos-Zc,2)+np.power(Rpos-Rc,2))**m_var_power;
			var_t 			= (((Zpos-Zc)**2+(Rpos-Rc)**2)**0.5)**m_var_power;
			var_t 			= np.transpose(var_t);   #Variance as a time vs probes matrix. Same shape as phi0
				
			shape 			= avgCOH.shape
			avgCOH 			= np.ravel(avgCOH)
			indices 		= np.nonzero(avgCOH)[0]
			avgCOH 			= avgCOH.reshape(shape)
			indices 		= indices.astype(int)       #Nonzero indices of avgCOH found
			pol_spect 		= spectograms[-len(R_p):]   #Spectograms created by poloidal probes
				
			PHASE0 			= np.zeros((len(indices),len(R_p)))
			for j in range(len(R_p)):
				pol_sp		= np.ravel(pol_spect[j])
				Fks 		= pol_sp[indices] #Unravel spect and consider nonzero Fks only
				PHASE0[:,j]	= np.angle(Fks)
			PHASE0 			= PHASE0[:,np.newaxis,:]
			
			#Find phi_nz, probe angles vs time for all nonzero indice times
			indices 		= np.nonzero(avgCOH)
			time_indices 	= np.sort((np.array(indices))[1])
			time_indices 	= time_indices.astype(int)
	
			#Create phi_nz and var_nz, these are probe angles and variance vs times for non-zero time indices
			phi0_nz 		= np.zeros(shape=PHASE0.shape)
			VAR 			= np.zeros(shape=PHASE0.shape)
			for i in range(PHASE0.shape[0]):
				tb_cur 		= time_indices[i]
				phi0_nz[i]	= phi0[tb_cur]
				VAR[i]		= var_t[tb_cur]
			time_indices 	= None
			
			#MMax			= NEST*MMaxFactor;
			MMax			= abs(NEST)*MMaxFactor; # make sure this can be positive for NEST<0
			MEst 			= nmnc.calculate_mode_num(phi0_nz,PHASE0,nzs,MMax,NMinAllowed=NEST,var=VAR)[2]
			
			#Broadcast back to avgCOH.shape
			MEST 			= np.zeros([avgCOH.size,1]);MEST[:]=np.nan
			MEST[nzs]		= MEst
			MEST 			= np.reshape(MEST,avgCOH.shape)
			print("\nPoloidal mode number calculation complete")	
			MEST_sp 		= csr_matrix(MEST)
	else:
		MEST 				= None
		MEST_sp 			= None
		
	##-----------------------------------------------------------------------------
	
		
#-------------------Averaged Fourier amplitude matrix-----------------------------			
			
	if amplitude_spect	== True:
		print("\n\nCREATING AVERAGED AMPLITUDE MATRIX")
		absolute_vals 	= np.abs(spectograms)  #Returns absolute values of Fourier coefficients
		avgAMPL 		= (np.sum(absolute_vals,axis=0))/len(spectograms)  #Averaged Fourier amplitude spectogram		
		avgAMPL  		= np.multiply(avgAMPL,avgCOH)
		AMPL_spect_sp 	= csr_matrix(avgAMPL)
		print("Averaged amplitude matrix found")
	else:
		avgAMPL  		= None
		AMPL_spect_sp 	= None
			
	##---------------------------------------------------------------------


	##----------------------ACQUIRE PLASMA PARAMETERS----------------------------
	
	if get_parameters 	== True:
		parameters_dict = get_params.get_parameters(params,pulse,time) #Dataframe for ONE PULSE. See get_parameters in Auxscripts
		
	##---------------------------------------------------------------------------


	##----------------------REMOVE MIN RESOLVABLE FREQUENCIES----------------------------
	
	#Try/excepts to deal with SKIP or if TMN/PMN/avgAMPL flags set to off
	
	try:
		min_fresbin 			= int(round(min_fres/df))
		avgCOH[:min_fresbin]	= 0	
		freqs[:min_fresbin]		= 0
		NEST[:min_fresbin]		= 0
		MEST[:min_fresbin]		= 0
		avgAMPL[:min_fresbin]	= 0
	except:						pass
		
	##---------------------------------------------------------------------		
		
		
	##----------------------SAMPLING FREQUENCY----------------------------
	
	#If some pulses have f_Nyquist = 500kHz, then we must restrict all spectrograms to this range. 
	#Perhaps a better way to do this, i.e. when spectrograms are intially compiled.
	
	try:
		print("\n\n Some probes sampling at 1000 kHz. Restricting spectrograms to 0-500kHz range.")
		
		sampling_freqs 	= np.array(sampling_freqs)
		if np.any(sampling_freqs == int(1e6)):
			f_nyq 		= 500000; 
			fnyqbin		= int(np.ceil(f_nyq/fnyqbin))
			avgCOH 		= avgCOH[:fnyqbin]
			avgAMPL		= avgAMPL[:fnyqbin]
			NEST 		= NEST[:fnyqbin]
			MEST 		= MEST[:fnyqbin]
			freqs 		= freqs[:fnyqbin]
	except:				pass

	##---------------------------------------------------------------------
	

	##----------------------SAVING SPECTROGRAM IMAGES----------------------------
	if image_save==True:
	#Saves avgCOH,avgAMPL,NEST, MEST spectrograms to image_folder given in main setup
	
		if num_pulses%image_save_freq==0:
			
			print("\n\n SAVING FIGURES")
			
			i=0
			
			spectrogram_list 	= [avgCOH,avgAMPL,NEST,MEST]
			savename_list 		= ["_averagedCoherence","_averagedAmplitude","_TMN","_PMN"]
			title_list 			= ["{} Averaged Coherence spectrogram ","{} Averaged Amplitude spectrogram ",
								   "{} Toroidal mode number ","{} Poloidal mode number "]
			cbarlabel_list 		= ["Coherence","Amplitude $Ts^{-1}$","",""];
			cbarlog_list 		= [False,True,False,False];
			max_list 			= [0,0,NMax-1,MMax-1];
			min_list 			= [0,0,NMin,MMin];
			
			for i in range(len(spectrogram_list)):
				plt.clf();
				savename 	= image_folder+"JPN_"+str(pulse)+savename_list[i];
				'''
				if i in [0,1]:
					impl.image(		spectrogram_list[i],time,freqs,
									title=title_list[i].format(pulse),
									min_f=im_minf,max_f=im_maxf,
									cbarlabel=cbarlabel_list[i],
									cbarlog=cbarlog_list[i],
									savepath=savename)
				elif i in [2,3]:
					impl.image_n(	spectrogram_list[i],time,freqs,
									title=title_list[i].format(pulse),
									min_f=im_minf,max_f=im_maxf,
									NMax=max_list[i],NMin=min_list[i],
									savepath=savename)
				'''
				if i==0:
					plt.clf()
					savename 	= image_folder+"JPN_"+str(pulse)+"_averagedCoherence"
					
					impl.image(	spectrogram_list[i],time,freqs,
								title="{} Averaged Coherence spectrogram ".format(pulse),
								min_f=im_minf,max_f=im_maxf,
								cbarlabel="Coherence",
								cbarlog=False,
								savepath=savename)
				elif i==1:
					
					plt.clf()
					savename 	= image_folder+"JPN_"+str(pulse)+"_averagedAmplitude"
					
					impl.image( spectrogram_list[i],time,freqs,
								title="{} Averaged Amplitude spectrogram ".format(pulse),
								min_f=im_minf,max_f=im_maxf,
								cbarlabel="Amplitude $Ts^{-1}$",
								cbarlog=True,
								savepath=savename)
				elif i==2:
					plt.clf()
					savename 	= image_folder+"JPN_"+str(pulse)+"_TMN"
					
					impl.image_n(	spectrogram_list[i],time,freqs,
									title="{} Toroidal mode number ".format(pulse),
									min_f=im_minf,max_f=im_maxf,
									NMax=NMax-1,NMin=NMin,
									savepath=savename)
				elif i==3:
					plt.clf()
					savename = image_folder+"JPN_"+str(pulse)+"_PMN"
					
					impl.image_n(	spectrogram_list[i],time,freqs,
									title="{} Poloidal mode number ".format(pulse),
									min_f=im_minf,max_f=im_maxf,
									NMax=3*NMax-1,NMin=MMin,
									savepath=savename)
				
		else:	pass
	##---------------------------------------------------------------------
	
		
	##----------------------SAVE ONLY NONZEROS----------------------------
	#Saving whole matrices requires large memory usage, so save only nonzeros which can then be broadcasted back
	
	try:
		Shape 		= avgCOH.shape
		NZS 		= np.nonzero(avgCOH) 	#Take nonzeros. NZS[0] ordered frequency indices, NZS[1] time indices (unordered)
		avgCOHnzs 	= avgCOH[NZS] 			#Row ordered nonzeros
		TMNnzs 		= NEST[NZS]
		PMNnzs 		= MEST[NZS]
		avgAMPLnzs = avgAMPL[NZS]
		NZS 		= np.vstack((NZS[0],NZS[1]))
	except:			pass
	
	##---------------------------------------------------------------------	

	
	##----------------------EXTRACT TO DATBASE ROW----------------------------
	#Consider using pd.categorical for mode number storage
	
	##Free parameters----------

	#--------------------
	if extract==True:
			
		if 1:# ALL == True:		
			new = {
				"Time array (s)": time,
				"Frequency array (Hz)": freqs,
				"NZ indices":NZS,
				"Averaged Coherence":avgCOHnzs,
				"TMN spectogram": TMNnzs,
				"PMN spectogram": PMNnzs,
				"Averaged Amplitudes":avgAMPLnzs,
				"Shape":Shape
				}
		else:	#If some flags off
			new = {
				"Time array (s)": time,
				"Frequency array (Hz)": freqs,
				"Averaged Coherence":avgCOHnzs,
				"Shape":np.array(avgCOH.shape)
				}
			
			
		if get_parameters==True:	new.update(parameters_dict)

		if num_pulses==1:			data_store={} #Initialise store
			
		data_store[str(pulse)]		= new #Add new KV pair
		
		print("\n\nSpectogram Database with parameters updated")

		if storing==True:
			#if periodic save true, saves periodically to mat file (and prelim DB if selected).
			if periodic_save == True:
				
				'''
				if num_pulses%save_freq==0 or pulse_counter==len(pulses): 
					
					keys 	= list(data_store.keys())
					keys.sort()
					print(keys)
					store 	= h5.File(storing_file,"a") #Open for appending
					print(str(len(list(store.keys())))+" keys in DATABASE file")
						
					for key in keys:
						f			= store.create_group(key)
						pulse_data 	= data_store[key]
						pulse_keys 	= list(pulse_data.keys())
						for dkey in pulse_keys:		f.create_dataset(dkey,data=pulse_data[dkey])
					
					store.close()
					
					with open(logs_file,"a") as f:
						f.write("\n\n"+str(save_freq)+" pulses starting from "+str(keys[0])+" and ending at "+str(keys[-1]))
						if num_pulses%save_freq==0: 	f.write("\n"+str(pulse_counter))
					print("\n\nSaved to storage file") 
					print("\n\nSaved "+str(len(keys))+" to storage file") 
					
					data_store.clear()
				'''
				
				
				if pulse_counter==len(pulses):   #FINAL PULSE CASE
					
					keys 	= list(data_store.keys())
					keys.sort()
					print(keys)
					store 	= h5.File(storing_file,"a") #Open for appending
					print(str(len(list(store.keys())))+" keys in DATABASE file")
						
					for key in keys:
						f			= store.create_group(key)
						pulse_data 	= data_store[key]
						pulse_keys 	= list(pulse_data.keys())
						for dkey in pulse_keys:		f.create_dataset(dkey,data=pulse_data[dkey])
					
					store.close()
					
					with open(logs_file,"a") as f:
						f.write("\n\n"+str(save_freq)+" pulses starting from "+str(keys[0])+" and ending at "+str(keys[-1]))

					print("\n\nSaved to storage file") 
					print("\n\nSaved "+str(len(keys))+" to storage file") 
					
					data_store.clear()
				
				
				
				if num_pulses%save_freq==0:
					keys 	= list(data_store.keys())
					keys.sort()
					print(keys)
					store 	= h5.File(storing_file,"a") #Open for appending
					#print(str(len(list(store.keys())))+" keys in DATABASE file")

					for key in keys:
						f 			=	store.create_group(key)
						pulse_data 	= data_store[key]
						pulse_keys 	= list(pulse_data.keys())
						for dkey in pulse_keys:
							f.create_dataset(dkey,data=pulse_data[dkey])
					
					store.close()
					
					with open(logs_file,"a") as f:
						f.write("\n\n"+str(save_freq)+" pulses starting from "+str(keys[0])+" and ending at "+str(keys[-1]))
						f.write("\n"+str(pulse_counter))
					
				
					print("\n\nSaved to storage file") 
					print("\n\nSaved "+str(len(keys))+" to storage file") 
					
					data_store.clear()
				
			
			
	#----------------------------------------------------------------------------
	
	if 1: #ADDITIONAL
		
		pulse_counter+=1  #Move onto next pulse
		
		if SKIP==True:			pass
		else:
			print(num_pulses)
			num_pulses+=1 #Number of pulses iterated through
			pulses_complete.append(pulse)
		
		#Break at certain pulse number
		if pulse_counter==2000:
			print('FINISHED')
			break	
			
		if test_pulse==True:	break #Break for one pulse when testing
		
		#RESET
		#SKIP = False
		#spect_AC=True;antenna_remove=True;ICRH_remove=True;tmn=True;pmn=True;amplitude_spect=True;get_parameters=True;image_save=True;extract=True;storing=True;periodic_save=True
				
			
	#--------------------------END----------------------------------			
			
print("\n\nDATABASE COMPILATION FINISHED")		
			
