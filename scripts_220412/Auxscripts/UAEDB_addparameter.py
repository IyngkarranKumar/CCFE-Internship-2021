"""
The purpose of this script is to add parameters to the UAE_DB. 
"""

if 1: ### Import #######################################################

	import numpy as np
	import matplotlib
	import matplotlib.pyplot as plt
	import matplotlib.mlab as mlab
	import scipy.signal as ssi
	from scipy.optimize import curve_fit
	from array import *
	import datetime
	import sys
	import pwd
	import os
	sys.path[:0]=['/jet/share/lib/python']
	from ppf import *
	import getdat
	sys.path[:0]=['/home/pgpuglia/AEAD/GUI/class']
	from AEAD_class import *
	from INCAA_functions_all import INCAA_Iant_Vant
	from scipy.linalg import solve
	from scipy.optimize import minimize
	import scipy.io as sio
	sys.path[:0] = ['/work/atinguel/PPF'];
	from ppfFunctions import *
	sys.path[:0] = ['/home/ht2059/Documents'];
	import h5py as h5
	from Auxscripts import readUAEDB as rDB
	from Auxscripts import get_parameters as get_params
	
	
if 1: #SETUPS
	pulses = rDB.pulse_list()
	#pulses = pulses[:1]
	nPulses = len(pulses)
	DB_path = "/home/ht2059/Documents/DATABASE_MAIN.h5"
	DATABASE = h5.File(DB_path,"a")
	
if 1: #PARAMETERS TO ADD
	#eg:
	params = [['ka2' ,'ka2-1','jetppf','ka2-1'], 						# Fast ion loss detector signal 1
			  ['ka2' ,'ka2-2','jetppf','ka2-2'], 						# Fast ion loss detector signal 2
			  ['ka2' ,'ka2-3','jetppf','ka2-3'], 						# Fast ion loss detector signal 3
			  ['ka2' ,'ka2-4','jetppf','ka2-4'], 						# Fast ion loss detector signal 4
			  ['ka2' ,'ka2-5','jetppf','ka2-5'], 						# Fast ion loss detector signal 5
			  ]
	

if 1: #UPDATE DATABASE WITH NEW PARAMETER DATA
	
	for pulse in pulses[1994:1995]:
		

		
		
		sub_db = rDB.find_pulse(pulse)[0]
		sub_database = DATABASE[sub_db]
		pulse_dset = sub_database[str(pulse)]
		
		timebase = pulse_dset["Time array (s)"]
		timebase = np.asarray(timebase)
		param_dset = get_params.get_parameters(params,pulse,timebase) #Return is dictionary, label: data, where labels are labels of new params, data is associated data
		
		for new_param in list(param_dset.keys()):
			print(len(list(pulse_dset.keys())))
			pulse_dset[new_param]=param_dset[new_param]
			
		print("\n PARAMETERS UPDATED FOR JPN "+str(pulse))
		print(len(list(pulse_dset.keys())))
		
if 1: #FINISH
	DATABASE.close()	
	print("\n\n DATABASE PARAMETERS UPDATED")	
		
