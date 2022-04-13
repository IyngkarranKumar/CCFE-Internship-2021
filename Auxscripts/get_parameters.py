"""The purpose of this script is to get a range of plasma parameters for a number of pulses and timebases that will be
entered into the Unstable Alfven Eigenmode database. For the parameters obtained, see below. Note that parameter array
MUST BE in the form ['dda', 'parameter', 'userid',('altname')]
It is heavily based off the following script: /home/atinguel/work/forIyngkarranK/save_parameters.py
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.signal as ssi
from scipy.optimize import curve_fit
import pandas as pd
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


# I will take as many signals as possible from https://users.euro-fusion.org/openwiki/index.php/List_of_Recommended_Signals
# In a list of ['dda', 'parameter', 'userid',('altname')]
"""params = [['magn','bvac','jetppf'  ], 								# Vacuum toroidal field at R=2.96 m (T)
		  ['magn','ipla','jetppf'  ], 								# Plasma current (A)
		  ['hrtx','ne0' ,'jetppf'  ], 								# Thomson scattering, electron density on axis, m-3 - chosen by Alex
		  ['HRTX','TE0' ,'jetppf'  ], 								# Thomson scattering, electron temperature on axis, eV?
		  ['efit','q95' ,'jetppf'  ], 								# q at 95% flux surface - chosen by Alex
		  ['efit','s95' ,'jetppf'  ], 								# shear (r/q)*(dq/dr) at 95% flux surface
		  ['efit','elon','jetppf'  ], 								# elongation - chosen by Alex
		  ['efit','xpfl','jetppf'  ], 								# limiter/diverted (x-point) flag - chosen by Alex 
		  ['NBI' ,'NBLM','jetppf'  ], 								# Neutral beam power, W
		  ['ICRH','PTOT','jetppf'  ], 								# The total ICRH power coupled to the plasma, W
		  ['efit','triu','jetppf'  ], 								# upper triangularity
		  ['efit','tril','jetppf'  ], 								# lower triangularity
		  ['efit','q'   ,'jetppf'  ,'q0'],							# q profile, but just get q0
		  ['efit','xlid','jetppf'  ], 								# li - diamagnetic
		  ['efit','xlim','jetppf'  ], 								# li - MHD
		  ['efit','btnm','jetppf'  ], 								# beta normalized (1)
		  ['elma','freq','chain1'  ,'felm'], 						# ELM frequency (Hz)
		  ['ks3b','hthd','jetppf'  ,'fH_ks3b'], 					# this is H/(H+D+T) so ~ n_H/n_e
		  ['kt5p','hthd','jetppf'  ,'fH_kt5p'], 					# this is H/(H+D+T) so ~ n_H/n_e
		  ['kt5p','dthd','jetppf'  ,'fD_kt5p'], 					# this is D/(H+D+T) so ~ n_D/n_e
		  ['kt5p','tttd','jetppf'  ,'fT_kt5p'], 					# this is T/(H+D+T) so ~ n_T/n_e
		  ['kt5b','che' ,'edelabie','fHe3_kt5b'], 					# this is n_He/n_e, so max is 0.5
		  ['kt5b','che' ,'edelabie','fHe4_kt5b'], 					# this is n_He/n_e, so max is 0.5
		  ['gash','c34c','jetppf'  ,'eflow_He3'], 					# flow of He3 electrons
		  ['gash','c35c','jetppf'  ,'eflow_He4'], 					# flow of He4 electrons
		  ['hrts','n95' ,'jetppf'  ], 								# Thomson scattering, electron density at q95, m-3 - chosen by Alex
		  ['hrts','gradn95','jetppf'], 								# Thomson scattering, gradient of n95
		  ];
nParams = len(params);
param_cols = [p[1] for p in params]

parameters_df = pd.DataFrame(columns=param_cols)"""


def get_parameters(params,pulse,timebase):
	
	#Params is array of arrays containing relevant parameter calling data. See above
	#Returns a dataframe of parameters in param_cols interpolated on given timebase for pulse.
	
	
	nParams = len(params);
	param_cols = [p[1] for p in params]
	#params_df = pd.DataFrame(columns=param_cols)
	params_df = {}
	
	t0 = round(timebase[0],2)
	tn1 = round(timebase[-1],2)
	print("\n\nGETTING PARAMETERS FOR PULSE "+str(pulse)+" for times "+str(t0)+" - "+str(tn1))
	

	for i in range(nParams): # for each parameter
			dda 	= params[i][0];
			param 	= params[i][1];
			userid 	= params[i][2];

			try: # to look for alternative label/name
				label = params[i][3];
			except: # if nothing, use param
				label = param;
				
			try:
				# Call specific functions if applicable
				if label 			== 's95':
					temp,t 			= gets95(pulse);
				elif label 			== 'q0':
					temp,t 			= getq0(pulse);
				elif label 			== 'fH_ks3b':
					temp,t,_,__ 	= getfH(pulse);
				elif label 			== 'fH_kt5p':
					_,__,temp,t 	= getfH(pulse);
				elif label 			== 'fD_kt5p':
					temp,t 			= getfD(pulse);
				elif label 			== 'fT_kt5p':
					temp,t 			= getfT(pulse);
				elif label 			== 'fHe3_kt5b':
					temp,_,t 		= getfHe(pulse);
				elif label 			== 'fHe4_kt5b':
					_,temp,t 		= getfHe(pulse);
				elif label 			=='NBLM':
					temp,t			= getNBI(pulse);
				elif label 			== 'PTOT':
					temp,t			= getICH(pulse);
				elif label 			== 'n95':
					temp,_,t		= getn95(pulse);
				elif label 			== 'gradn95':
					_,temp,t		= getn95(pulse);
				elif 'angf' in label:
					temp,t			= getangf(pulse,dda);
					
				elif 'ka2' in label:
					temp,t 			= getka2(pulse,param);
				else:
					temp,x,t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = \
					ppfdata(pulse,dda,param,seq=0,uid=userid,device="JET",fix0=0,reshape=1,no_x=0,no_t=0,no_data=0);
			except:
				temp = [];t=[]
		
			#Check if array returned is emtpy
			if len(temp) == 0 or len(t) == 0:
				t 	 = np.array([np.NaN]);
				temp = np.array([np.NaN]);	
				print(" NOTE: No data returned for "+label+" parameter.")
				
			#Interpolation:
			temp_time = np.interp(timebase,t,temp,left=np.nan,right=np.nan);
			
			
			params_df[label]=temp_time
			
			


	return params_df

	
	
	
	
