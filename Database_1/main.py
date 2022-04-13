"""
The following script compiles Database I - discussed further in attached report.
This file mainly comprises of setting various parameters for compilation

General structure:
    - Imports
    - Set node and probe parameters to call data from JET servers 
    - Set parameters to obtain relevant plasma data (e.g heating power)
    - Set parameters for various script functionalities
    - Set pulse parameters and saving files
    - Main loop


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
	#import pdb #debugging
	#pdb.set_trace() #debugging


	from Auxscripts import preprocessing as prep
	from Auxscripts import plotting as impl
	from Auxscripts import new_mode_number_calculation as nmnc
	from Auxscripts import get_ICRH_signal as get_ICRH
	from Auxscripts import get_parameters as get_params
 
  #General algo

#-----------------------NODE AND PROBE SETTINGS---------------------------
if 1:
	node 	= 'DA/C1M-';    #Node from which to collect data from



	test_pr = False; k=8


	#Toroidal probes and angles (radians)
	toroidal_probes_main = ['T001', 'T002', 'T006', 'T007', 'T008', 'T009',
			   'H301', 'H302', 'H303', 'H304', 'H305', 'H306', 'H307', 'H308'];

	phi_t = [3.0,42.2,182.9,222.2,257.1,290.4,77.0,92.9,103.1,108.4,110.4]; 
	phi_t_main = np.radians(phi_t) #radians
	n_torprobes_main = len(toroidal_probes_main)


	#Poloidal probes and location (R,Z).
	poloidal_probes_main = ["PP801","PP802","PP803","PP804","PP805","PP806","PP807"]
	R_p_main = np.array([3.766,3.843,3.952,3.983,3.863,3.673,3.450]) #   Radial location of probes. 
	Z_p_main = np.array([1.249,1.075,0.694,0.184,-0.383,-0.770,-1.057]) #Vertical axis location of probes
	
	
	#Angle measurement uncertainty
	thetaUNC = np.radians(30)

	n_polprobes_main = len(poloidal_probes_main)

	probes_main = toroidal_probes_main+poloidal_probes_main

##------------------------------------------------------------------------


#----------------PLASMA AND OPERATIONAL PARAMETERS------------------------
if 1:
	#To add param, follow structure below
	params = [['magn','bvac','jetppf'  ], 								# Vacuum toroidal field at R=2.96 m (T)
				['magn','ipla','jetppf'  ], 								# Plasma current (A)
				['hrtx','ne0' ,'jetppf'  ], 								# Thomson scattering, electron density on axis, m-3 - chosen by Alex
				['HRTX','TE0' ,'jetppf'  ], 								# Thomson scattering, electron temperature on axis, eV?
				['efit','q95' ,'jetppf'  ], 								# q at 95% flux surface - chosen by Alex
				['efit','s95' ,'jetppf'  ], 								# shear (r/q)*(dq/dr) at 95% flux surface
				['efit','elon','jetppf'  ], 								# elongation - chosen by Alex
				#['efit','xpfl','jetppf'  ], 								# limiter/diverted (x-point) flag - chosen by Alex 
				['NBI' ,'NBLM','jetppf'  ], 								# Neutral beam power, W
				['ICRH','PTOT','jetppf'  ], 								# The total ICRH power coupled to the plasma, W
				['efit','triu','jetppf'  ], 								# upper triangularity
				['efit','tril','jetppf'  ], 								# lower triangularity
				['efit','q'   ,'jetppf'  ,'q0'],							# q profile, but just get q0
				#['efit','xlid','jetppf'  ], 								# li - diamagnetic
				#['efit','xlim','jetppf'  ], 								# li - MHD
				['efit','btnm','jetppf'  ], 								# beta normalized (1)
				#['elma','freq','chain1'  ,'felm'], 						# ELM frequency (Hz)
				['ks3b','hthd','jetppf'  ,'fH_ks3b'], 					# this is H/(H+D+T) so ~ n_H/n_e
				['kt5p','hthd','jetppf'  ,'fH_kt5p'], 					# this is H/(H+D+T) so ~ n_H/n_e
				['kt5p','dthd','jetppf'  ,'fD_kt5p'], 					# this is D/(H+D+T) so ~ n_D/n_e
				['kt5p','tttd','jetppf'  ,'fT_kt5p'], 					# this is T/(H+D+T) so ~ n_T/n_e
				['kt5b','che' ,'edelabie','fHe3_kt5b'], 					# this is n_He/n_e, so max is 0.5
				['kt5b','che' ,'edelabie','fHe4_kt5b'], 					# this is n_He/n_e, so max is 0.5
				#['gash','c34c','jetppf'  ,'eflow_He3'], 					# flow of He3 electrons
				#['gash','c35c','jetppf'  ,'eflow_He4'], 					# flow of He4 electrons
				['hrts','n95' ,'jetppf'  ], 								# Thomson scattering, electron density at q95, m-3 - chosen by Alex
				['hrts','gradn95','jetppf'], 								# Thomson scattering, gradient of n95
				['cxd6','angf','jetppf','angf_cxd6'], 					# Rotation angular frequency (rad/s)
				['cxg6','angf','jetppf','angf_cxg6'], 					# Rotation angular frequency (rad/s)
				['cxh4','angf','jetppf','angf_cxh4'], 					# Rotation angular frequency (rad/s)
				['cxhm','angf','jetppf','angf_cxhm'], 					# Rotation angular frequency (rad/s)
				['cxkm','angf','jetppf','angf_cxkm'], 					# Rotation angular frequency (rad/s)
				];
	nParams = len(params);
	param_cols = [p[1] for p in params]

	parameters_df = pd.DataFrame(columns=param_cols)

##------------------------------------------------------------------------


#----------------SCRIPT FUNCTIONALITIES & FREE PARAMETERS----------------

if 1:

	#tStart = float(input("Enter start time: "))
	#tStop = float(input("Enter end time: "))
	#pulses = list(input("Enter Pulse list: "))


	#---------------Generate list of spectograms and averaged Coherence spectogram. Required true for all other functionalities

	#Free parameters (NOTE: some FPs dependent on JPN, so defined in loop. They have been commented out here).--------

	threshold = 0.6     #avgCOH coherence threshold
	take_freq = 4     #In list of NC2 pairings of N probe spectrograms, how often to choose pairing to include in avgCOH calc. Low take_freq --> stronger Coherence
	#noverlap = N_bin//4   #Overlap quarter of previous bin for spectogram creation
	coh_n = 5   		#Time bin overlap for coherence matrix creation
	T_bin = 2.048e-3
	df = 1//T_bin
	min_fres = 1/T_bin #Minimum resolvable frequency from Fourier Transform
	#-----------------------------



	#------------Antenna signal remove 
	#Free parameters (NOTE: some FPs dependent on JPN, so defined in loop. They have been commented out here)---------

	#fb_width = int(1500//df)   #Frequency width of antenna signal
	#tb_width = int(0.01//T_bin)  #Time bin width of antenna signal
	Ithreshold 	= 2;   #Antenna current threshold (we only consider input currents of 2A or greater)
	fpres = 1         #Used to retain signal that intersects with TAEs
	tpres = 25		  #Used to retain signal that intersects with TAEs
	af_threshold = 0.5 #Removes bins that have not been completely cut out by removing reconstructed antenna signal
	#------------------


	#-----------ICRH signal remove

	#---------------------

	#------------Toroidal mode number spectogram

	NMax = 8; NMin = -1 
	#-----------------------

	#--------------Poloidal mode number

	MMax = 3*NMax;MMin=NMin
	#-------------------

	#----------------Averaged Amplitude spectgram

	#------------------------
	
	#----------------Saving spectrograms

	
	#------------------------

	#----------------Acquire plasma parameters

	#------------------------

	#-------------------Extract to Database
	#extract is general extraction (can extract individual modes and/or spectograms),storing is storing data in hdf5 file, periodic save - Save every save_freq pulses.

	#------------------------
	#(NOTE: some FPs dependent on JPN, so defined in loop. They have been commented out here).




	#To skip code. This is used when most probes do not return signal. KEEP AS FALSE
	SKIP = False

##------------------------------------------------------------------------
	

##--------------------------MAIN SETUP------------------------------------

if 1: #Pulse setup prior to main loop. Option for one test pulse here
	
	"""
	IMPORTANT:
		#To alter the fast magnetics probe config, see -NODE AND PROBE SETTINGS
		#To alter the parameter time-series incorporated into the DB, see -PLASMA AND OPERATIONAL PARAMETERS
		#To alter the main script functionalities and associated free parameters, see -SCRIPT FUNCTIONALITIES & FREE PARAMETERS
	"""
	
	##FUNCTIONALITIES
	spect_AC=True;antenna_remove=True;ICRH_remove=True;tmn=True;
	pmn=True;fixed_centroid=False;amplitude_spect=True;get_parameters=True;image_save=True;
	extract=True;storing=True;periodic_save=True
	
	
	#Saving settings
	#HD5 file path
	storing_file = "DATABASE_testing.h5"
	save_freq = 10 #Save to hdf every `save_freq` pulses
	#Logs path
	logs_file = "/home/ht2059/Documents/pulses_completed_testing.txt"
	#Image folder
	image_folder = "/home/ht2059/Documents/plots/"
	image_save_freq = 10;im_minf=0;im_maxf=500  #Save every image_save_freq pulse. Choose frequency range for PLOTS ONLY
	
	
	
	
	#Pulses to analyse here. 
	pulse_data = sio.loadmat("/home/atinguel/work/forIyngkarranK/pulses_w_ICRH_93000-98006_Pthreshold500kW_Tthreshold500ms")
	pulses = pulse_data["pulses"][0]
	tStarts = pulse_data["tStart"][0]
	tStops = pulse_data["tEnd"][0]
	
	
	
	
	
	pulse_counter = 0  #Set to k-1 to start at kth pulse
	
	#OR
	
	#pulse_counter = np.where(pulses==94700)[0][0] 
	
	

	
	test_pulse = False
	if test_pulse==True:
		pulse = 96851
		pulse_counter = np.where(pulses==pulse)[0][0] #1 specific pulse
		
	test_timings = True #True if want to compile for test_time_length interval (for debugging mainly)
	test_time_length = 0.5


##------------------------------------------------------------------------	
		

##--------------------------MAIN LOOP------------------------------------

num_pulses=1 #Number pulses iterated through
pulses_complete = [] #List of completed pulses

for pulse in pulses[pulse_counter:]:
    

    if 1:
            
        toroidal_probes = deepcopy(toroidal_probes_main)
        phi_t = deepcopy(phi_t_main)
        ntorprobes = len(toroidal_probes)


        poloidal_probes = deepcopy(poloidal_probes_main)
        R_p = deepcopy(R_p_main)
        Z_p = deepcopy(Z_p_main)
        npolprobes = len(poloidal_probes)

        probes = toroidal_probes+poloidal_probes

	
	
    print("\n\n\n\n\n\n\n\n\nPulse "+str(pulse_counter)+" of "+str(len(pulses))+" pulses.")
    print("\n\nCURRENT JPN: "+str(pulse))
    #if pulse_counter==10:
    #	break

    T_bin = T_bin
	
	
	#Defining pulse dependent parameters (tStart, tSTop etc)
    if 1:

        tStart = tStarts[pulse_counter]
        tStop = tStops[pulse_counter]

        #dt=1/2000000
        #tStart = round(tStarts[pulse_counter]+10*dt,2)
        #tStop = round(tStops[pulse_counter]-10*dt,2)


        if test_timings==True:
        #TEST TIMES
            tStop_orig = tStop
            tStop = tStart+test_time_length
   
    
    if spect_AC==True:
        
        spectrograms = get_spectrograms(pulse)
        
        avgCOH,time = get_avgCOH_time(pulse)
        
        t_bin = t_bin = ((time[-1]-time[0])/(avgCOH.shape[1])) #Effective T_bin.
        
        
    
    if antenna_remove==True:
        
        
        
        
    
        
        
         