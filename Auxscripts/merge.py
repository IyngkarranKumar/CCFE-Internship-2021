"""
This script defines a function that merges significant structures identified by the mode extraction algorithm. 

"""
if 1: #IMPORT
	import numpy as np
	from copy import deepcopy
	import itertools as it



def merge(SIG_STRS_ALL,t_o,f_o):
    
    



	num = len(SIG_STRS_ALL)
	el = np.arange(0,num,1)
	arr = np.arange(num)
	pairs = list(it.combinations(el,2))

	connected = []
	SIG_STRS_ALL_NEW=[]
	i=0
	for pair in pairs:
		SS1 = SIG_STRS_ALL[pair[0]]
		SS2 = SIG_STRS_ALL[pair[1]]
		SS1_time = set(SS1[:,0]) ; SS1_freq = set(SS1[:,1])
		SS2_time = set(SS2[:,0]) ; SS2_freq = set(SS2[:,1])
		
		S1_time_arr = SS1[:,0] ; S1_freq_arr = SS1[:,1]
		S2_time_arr = SS2[:,0] ; S2_freq_arr = SS2[:,1]
		
		S1_mt = np.mean(S1_time_arr) ; S1_mf = np.mean(S1_freq_arr)
		S2_mt = np.mean(S2_time_arr) ; S2_mf = np.mean(S2_freq_arr)
		t_diff = np.abs(S1_mt-S2_mt); f_diff = np.abs(S1_mf-S2_mf)
		t_mn = np.mean([(S1_time_arr[-1]-S1_time_arr[0]),(S2_time_arr[-1]-S2_time_arr[0])]); f_mn = np.mean([S1_mf,S2_mf]);
		f_mn = np.mean([(S1_freq_arr[-1]-S1_freq_arr[0]),(S2_freq_arr[-1]-S2_freq_arr[0])])
		
		nt_overlap = len(SS1_time.intersection(SS2_time))
		nf_overlap = len(SS1_freq.intersection(SS2_freq))
		
		sz = np.min([len(SS1),len(SS2)])

		

		
		t_t= t_o*(t_mn/2e-3); f_t = f_o*(f_mn/488)
		if t_diff<t_t and f_diff<10:#Those which overlap `enough' are added to connected to be stitched together
			connected.append(pair)
			
	#Connected is list of all connected pairings
	connected_list = [list(pair) for pair in connected]
	connected_list = np.ravel(np.array(connected_list))
	connected_list = np.unique(connected_list)

	set_diff = np.setdiff1d(arr,connected_list)


	#Add those with no/little connection to new SIG_STRS_ALL
	for j in set_diff:
		SIG_STRS_ALL_NEW.append(SIG_STRS_ALL[int(j)])
		

	def network(L): #Finds connected structures
		
		olength = len(L)
		if type(L)!=list:
			raise TypeError("Input must be list")
			
		web=[]
		for link in L:
			for pair in connected:
				if link in pair:
					web=web+list(pair)
		
		web = list(set(web))
		
		if olength-len(web)==0:
			return list(web)
			
		
		
		else:
			return network(web)
	  
	  
	  
	#Finally,fill SIG_STRS_ALL_NEW with merged arrays
	done = []
	for i in connected_list:
		if i in done:
			continue
		else:
			net = network([i])  #Find all connected arrays
			done=done+net
			arr_0 = SIG_STRS_ALL[net[0]]
			for i in net[1:]:
				arr_0 = np.vstack((arr_0,SIG_STRS_ALL[i]))
			SIG_STRS_ALL_NEW.append(arr_0)
			


	SIG_STRS_ALL = SIG_STRS_ALL_NEW

	return SIG_STRS_ALL



def overlap(SIG_STRS_ALL):
	print("in")
	SIG_STRS_ONE =[]
	for nstruct in SIG_STRS_ALL:
		for struct in nstruct:
			SIG_STR_ONE = SIG_STRS_ONE.append(np.array(struct))
	SIG_STRS_ALL = deepcopy(SIG_STRS_ONE)
	SIG_STRS_ONE=None



	num = len(SIG_STRS_ALL)
	print(num)
	el = np.arange(0,num,1)
	arr = np.arange(num)
	pairs = list(it.combinations(el,2))

	connected = []
	SIG_STRS_ALL_NEW=[]
	i=0
	for pair in pairs:
		SS1 = SIG_STRS_ALL[pair[0]]
		SS2 = SIG_STRS_ALL[pair[1]]
		SS1_time = set(SS1[:,0]) ; SS1_freq = set(SS1[:,1])
		SS2_time = set(SS2[:,0]) ; SS2_freq = set(SS2[:,1])
		
		nt_overlap = len(SS1_time.intersection(SS2_time))
		nf_overlap = len(SS1_freq.intersection(SS2_freq))
		overlap = nt_overlap+nf_overlap
		
		print(overlap)
        
        
