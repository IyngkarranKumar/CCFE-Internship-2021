''' The purpose of this script is to define functions which calculate the
	toroidal mode number for Iyngkarran. Copied from ppfFunctions.py.
	Alex Tinguely 210615
'''

if 1: ### Import ###
	import numpy as np
	import matplotlib.pyplot as plt
	from copy import deepcopy
		
# Calculate chi-square for mode number spectrum
def calculate_mode_num(phi0,PHASE0,nzs,NMaxAllowed,NMinAllowed=[],var=[],phaseUnc=30*np.pi/180,Ip=-1):
	# phi0 			= angle of probes (rad), can be an array or matrix
	# PHASE0 		= phase of each probe (rad). shape=(l,1,phi). l = len(nonzero_indices)
	# NMaxAllowed 	= maximum N that we want to fit (even if higher mode numbers are resolvable), can be matrix
	# NMinAllowed 	= minimum N that we want to fit, can be empty, can be matrix
	# var 			= variance of data for weighting, can be empty, can be a matrix
	# phaseUnc 		= uncertainty in phase, only needed for resolvable mode numbers (not used here)
	# Ip 			= plasma current, really only need sign to determine direction of n>0

	# Get shape of phase matrix
	shape0 	= np.shape(PHASE0)[0]; # number of non-zero indices, l
	shape1 	= np.shape(PHASE0)[1]; # 1
	shape2 	= np.shape(PHASE0)[2]; # should be same size as phi, number of probes

	if len(var)==0:				var = np.ones(shape=np.shape(phi0));
	if type(NMinAllowed)==list: NMinAllowed=-NMaxAllowed;	
	
	pi 		= np.pi;
	val 	= -1; # 'nonsensical' value for not fittable mode numbers
	if len(np.shape(NMaxAllowed))==0:
		nmax = NMaxAllowed
	else:
		
		nmax = np.amax(NMaxAllowed[np.logical_not(np.isnan(NMaxAllowed))])
		NMaxAllowed_RS = np.zeros(len(nzs))
		NMaxAllowed_RS[:] = np.ravel(NMaxAllowed)[nzs] #Reshape the NMaxAllowed input matrix to dimensions lxn
	if len(np.shape(NMinAllowed))==0:
		nmin = NMinAllowed
	else:
		nmin = np.amin(NMinAllowed[np.logical_not(np.isnan(NMinAllowed))])
		NMinAllowed_RS = np.zeros(len(nzs))
		NMinAllowed_RS[:] = np.ravel(NMinAllowed)[nzs] #Reshape the NMaxAllowed input matrix to dimensions lxn


	NRange 	= np.arange(nmin,nmax+1);
	
	# Initialize X2 matrix
	X2 	 = np.zeros(shape=(shape0,shape1,len(NRange))); #Shape1=1
	
	### Make grids the size of phase ###
	# for phi
	Nphi 	= np.shape(phi0)[-1];
	dimphi 	= len(np.shape(phi0));
	if dimphi < 3:	_,__,PHI0	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phi0,indexing='ij');
	else:			PHI0		= deepcopy(phi0);
	
	# for var
	dimvar 	= len(np.shape(var));		
	if dimvar < 3:	_,__,VAR	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),var,indexing='ij');
	else:			VAR			= deepcopy(var);
	
	
		
	# for nrange
	NN 		= len(NRange);
	_,__,NRANGE = np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),NRange,indexing='ij');
		
	# for nmin
	dimmin 	= len(np.shape(NMinAllowed));
	if dimmin == 0: 					NMIN = np.ones(np.shape(NRANGE))*NMinAllowed;
	elif dimmin == 2:	
		"""NMIN = np.stack((NMinAllowed,NMinAllowed),axis=2);	
		for iN in range(len(NRange)-2): NMIN = np.concatenate((NMIN,NMinAllowed[:,:,np.newaxis]),axis=2);"""
		NMIN = NMinAllowed_RS[:,np.newaxis,np.newaxis]

	else:								NMIN = deepcopy(NMinAllowed);
		
	# for nmax
	dimmax 	= len(np.shape(NMaxAllowed));
	if dimmax == 0: 					NMAX = np.ones(np.shape(NRANGE))*NMaxAllowed;
	elif dimmax == 2: 	
		"""NMAX = np.stack((NMaxAllowed,NMaxAllowed),axis=2);
		for iN in range(len(NRange)-2): NMAX = np.concatenate((NMAX,NMaxAllowed[:,:,np.newaxis]),axis=2);"""
		NMAX = NMaxAllowed_RS[:,np.newaxis,np.newaxis]
	
	else:								
		NMAX = deepcopy(NMaxAllowed);


	# Sort the data by toroidal angle
	iSort 	= np.argsort(PHI0,axis=2);
	PHI0 	= np.take_along_axis(PHI0,iSort,2);
	PHASE0 	= np.take_along_axis(PHASE0,iSort,2);
	VAR		= np.take_along_axis(VAR,iSort,2);

	for j in range(Nphi): # try fitting with each probe shifted to "zero"
		# Shift the data so that phi[j] = 0 and phase[j] = 0
		#print('...shift '+str(j)+' of '+str(Nphi));
		PHI 	= deepcopy(PHI0);
		PHASE 	= deepcopy(PHASE0);
		for k in range(Nphi): 
			PHI[:,:,k] 	-= PHI0[:,:,j];
			PHASE[:,:,k]-= PHASE0[:,:,j];

		# And put phase in range [0, 2pi)
		PHI		= PHI % (2*pi);
		PHASE 	= PHASE % (2*pi); 

		for i in range(len(NRange)):
			NTemp 		= NRange[i];
			PHASETemp 	= NTemp*PHI;
			PHASETemp 	= PHASETemp % (2*pi);
			
			# Take the difference squared, but consider the minimum distance when adding/subtracting 2pi
			dPHASE2 = np.stack([ (PHASETemp        - PHASE)**2,
								 (PHASETemp + 2*pi - PHASE)**2,
								 (PHASETemp - 2*pi - PHASE)**2],
								 axis=3);
			dPHASE2 = np.min(dPHASE2,axis=3);
			X2[:,:,i]	= X2[:,:,i] + np.sum(dPHASE2/VAR,axis=2)/np.sum(1/VAR,axis=2)/Nphi; # chi-square normalized by number of probes

	X2 = X2/Nphi; # normalize again by number of probes (since we iterated for different "zeros")

	if Ip < 0: # N defined in the direction of Ip (which is usually negative in JET)
		X2 = np.flip(X2,axis=2);

	boolMin 			= NRANGE >= NMIN;
	boolMax 			= NRANGE <= NMAX;
	boolResolvable 		= np.logical_and(boolMin,boolMax);
	X2[~boolResolvable] = 1e10;	# set X2 values to huge number so that they can't be the min
	
	# Get estimate of toroidal mode number
	iMin 		= np.argmin(X2,axis=2);
	NRESOLVABLE = NRANGE;
	NEst 		= np.take(NRESOLVABLE,iMin);

	return X2, NRange, NEst, NRESOLVABLE, iMin
	
"""
if 1: ### Test ###
	shape0 			= 2;
	shape1 			= 3;
	shape2 			= 10;

	#nMax			= 10; # maximum tor mode number to consider
	#nMin 			= -10; # minimum tor mode number to consider

	n 				= -5; # toroidal mode number
	phit 			= np.random.rand(shape2)*2*np.pi;
	var 			= np.random.rand(len(phit))*0.1;
	phase 			= (n*phit + np.random.rand(len(phit))*1) % (2*np.pi); # rad, phase values with random noise in [0,2pi)
	_,__,PHI0 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phit,indexing='ij');
	_,__,VAR0 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),var,indexing='ij');
	_,__,PHASE0 	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phase,indexing='ij');

	n 				= 4; # toroidal mode number
	phit 			= np.random.rand(shape2)*2*np.pi;
	var 			= np.random.rand(len(phit))*0.1;
	phase 			= (n*phit + np.random.rand(len(phit))*1) % (2*np.pi); # rad, phase values with random noise in [0,2pi)
	_,__,PHI1 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phit,indexing='ij');
	_,__,VAR1		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),var,indexing='ij');
	_,__,PHASE1 	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phase,indexing='ij');

	#phit 			= np.array([2,68,119,334,277,191,15,200])*np.pi/180;
	phit 			= np.random.rand(len(phit))*2*np.pi;
	var 			= np.random.rand(len(phit))*0.9;
	n 				= 1; # toroidal mode number
	phase 			= (n*phit + np.random.rand(len(phit))*1) % (2*np.pi); # rad, phase values with random noise in [0,2pi)
	_,__,PHI2 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phit,indexing='ij');
	_,__,VAR2 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),var,indexing='ij');
	_,__,PHASE2 	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phase,indexing='ij');

	n 				= -2; # toroidal mode number
	phase 			= (n*phit + np.random.rand(len(phit))*1) % (2*np.pi); # rad, phase values with random noise in [0,2pi)
	_,__,PHI3 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phit,indexing='ij');
	_,__,VAR3 		= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),var,indexing='ij');
	_,__,PHASE3 	= np.meshgrid(np.linspace(0,1,shape0),np.linspace(0,1,shape1),phase,indexing='ij');

	PHI 			= np.hstack((np.vstack((PHI0,PHI1)),np.vstack((PHI2,PHI3))));
	VAR 			= np.hstack((np.vstack((VAR0,VAR1)),np.vstack((VAR2,VAR3))));
	PHASE 			= np.hstack((np.vstack((PHASE0,PHASE1)),np.vstack((PHASE2,PHASE3))));
	
	NMIN0 			= -10;
	NMAX0 			= 10;
	#NMIN0 			= np.random.randint(-10,-1,size=(np.shape(PHI)[0],np.shape(PHI)[1]));
	#NMAX0 			= np.random.randint(1,10,size=(np.shape(PHI)[0],np.shape(PHI)[1]));

	# Calculate estimated mode number from minimum chi-square (X2)
	# Let plasma current be in positive direction Ip > 0
	X2,nRange,NEST,NRESOLVABLE,iMin = calculateN_X2(PHI,PHASE,NMAX0,Ip=1,NMinAllowed=NMIN0,var=VAR)
"""
