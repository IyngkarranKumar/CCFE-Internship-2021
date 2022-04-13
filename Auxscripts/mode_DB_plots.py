"""

"""


savefile = "/home/ht2059/Documents/IDENTIFIED_MODES_50khz_500khz.csv"
save_to = "/home/ht2059/Documents/analysis_plots/mode growth time"
parameters = ["bvac","ne0","n95","ipla","TE0","s95","NBLM","PTOT","q0","q95","btnm","ka3_ccds","gradn95"]
pkl = "/home/ht2059/Documents/IM_DATABASE.pkl"

if 1: #IMPORT
	import numpy as np
	import sys
	import matplotlib
	import matplotlib.pyplot as plt
	import matplotlib
	from matplotlib.pyplot import cm
	import pandas as pd
	from random import randint
	import h5py as h5
	from scipy import io as sio
	from numpy.polynomial import polynomial as poly
	import matplotlib.font_manager as font_manager



	from Auxscripts import readUAEDB as rDB




def num_vs(param,DB_path):
	num_modes=[]
	data = []
	main_DB_path = "/home/ht2059/Documents/IM_DATABASE_0_500khz.pkl"
	mdb = pd.read_pickle(DB_path)


	
	for row in mdb:
		par = row[param]
		
		data.append(par)

	
		
		
	data = np.array(data)
	data=data[~np.isnan(data)]
	
	if 0: #Linear fit
		m,c = poly.polyfit(pulse_param,num_modes,1)
		
		xr = np.arange(np.amin(pulse_param),np.amax(pulse_param),100)
		yr = m*xr+c
		
		plt.plot(xr,yr,"r")
		
	plt.plot(pulse_param,num_modes,"bo")
	plt.xlabel(param)
	plt.ylabel("Number of modes identified")
	plt.title("Number of modes vs "+str(param))
	plt.show()





def power_ratio(DB_path,hist=False,min_f=0,max_f=500):
	
	modedb = pd.read_pickle(DB_path)
	boolf = np.zeros(len(modedb),dtype=bool)
	
	#Setting boolean array and filtering database
	if 1:
		for i in range(len(modedb)):
			f_list = modedb.iloc[i]["F"]
			f_med = f_list[1]
			if f_med >= min_f and f_med <=max_f:
				boolf[i]=True
			else:
				pass
				
		modedb = modedb[boolf]
			
		
		
	ICRH = modedb["PTOT"]; NBI = modedb["NBLM"]
	ICRH = [i[2] for i in ICRH]; NBI = [N[2] for N in NBI]
	ICRH = np.array(ICRH) ; NBI = np.array(NBI)
	ICRH_nan = ICRH==np.nan; NBI_nan = NBI==np.nan
	
	if len(NBI_nan)>len(ICRH_nan):
		inds=NBI_nan
		NBI = NBI[~NBI_nan]
		ICRH = ICRH[~NBI_nan]
	else:
		inds=ICRH_nan
		NBI = NBI[~ICRH_nan]
		ICRH = ICRH[~ICRH_nan]
		
	
	p_ratio = ICRH/NBI
	
	if 0:  #Not relevant anymore
		##Linear Fit
		c,m = poly.polyfit(power_ratio,num_modes,1)
		xr = np.arange(np.amin(power_ratio),np.amax(power_ratio),100)
		yr = m*xr+c
	
	
	plt.clf()
	
	if hist==True:
		bins=np.logspace(-1,np.log10(4),21)
		N, bin_edges,patches = plt.hist(param_data,bins=bins,log=False)
		bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

		if 0:
			plt.errorbar(
				bin_centers,
				N,
				yerr = N**0.5,
				marker = '.',
				drawstyle="steps-mid",
				elinewidth=0.10,
				ecolor="#682117",
				barsabove=True
			)
			
			plt.title("Identified modes 0-1000khz, ICRH/NBI distribution")
			plt.text(bin_centers[-4],np.amax(N),"N = "+str(len(param_data)))
			plt.grid(True)
			plt.xscale("log")
			plt.xlabel("ICRH/NBI")
			plt.ylabel("Number of modes")
		
	else:
		A = modedb["A"]
		AMPL = A[~inds]
		
		if 1: #Colour
			def colour(c):
				c = np.array(c); c = c/255; c = tuple(c)
				return c
		
		
		plt.clf()
		plt.scatter(p_ratio,AMPL,alpha=0.5,s=5)
		#plt.plot(xr,yr,"r")
		plt.ylabel("Mode amplitude (T/s)")
		plt.yscale("log")
		plt.ylim([10e-8,10e-1])
		plt.xlabel("ICRH/NBI power")
		plt.xscale("log")
		plt.xlim([10e-3,10e2])
		plt.title("Mode amplitude vs ICRH/NBI 100-300khz")
		plt.grid(True)
		plt.show()
		#plt.savefig("/home/ht2059/Documents/analysis_plots/AMPLITUDEdist/p_ratio")
	

		
	
	plt.show()
	
	return param_data
	
	
def num_hist(DB_path,param,title,xlabel,ylabel,savepath=None,min_f=0,max_f=500):
	
	modedb = pd.read_pickle(DB_path)
	boolf = np.zeros(len(modedb),dtype=bool)
	
	#Setting boolean array and filtering database
	if 1:
		for i in range(len(modedb)):
			f_list = modedb.iloc[i]["F"]
			f_med = f_list[1]
			if f_med >= min_f and f_med <=max_f:
				boolf[i]=True
			else:
				pass
				
		modedb = modedb[boolf]
		
	if param=="n" or param=="m":

		param_data= modedb[param]
		mn = []
		for i in range(len(modedb)):
			p = param_data.iloc[i]
			p = p[1]
			mn.append(p)
		param_data_100ms=np.array(mn)
		bins = np.arange(np.amin(param_data)-0.5,np.amax(param_data)+1.5,1);
			
		
			
	else:
		bins = 30;	
		param_data = modedb[param]
		param_data_100ms = np.array([p[2] for p in param_data]) #Params with t_prio=100ms
		param_data_100ms = param_data_100ms[~np.isnan(param_data_100ms)]
	
	
	units = rDB.units()
	
	def colour(c):
		c = np.array(c); c = c/255; c = tuple(c)
		return c

	param_DATA = param_data_100ms
	fig,ax = plt.subplots()
	param_DATA = np.round(param_DATA,2)
	#param_DATA = np.random.rand(10000)
	A,__,__ = ax.hist(param_DATA,bins=bins,log=True);ax.clear()
	N, bin_edges,patches = ax.hist(param_DATA,bins=bins,
									log=True,
									density=True,
									color=colour((102,0,128)),
									edgecolor="black",
									alpha=0.8)
	
	bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
	
	A = np.array(A); 
	N = np.array(N)
	

	font_prop = font_manager.FontProperties(size=12)

	
	if 1:
		ax.errorbar(
			bin_centers,
			N,
			yerr = N**0.5,
			marker = '.',
			drawstyle="steps-mid",
			elinewidth=1.5,
			capsize=3,
			ecolor=(0,0,0),
			barsabove=True,
		)
	
	
	
	ax.set_xlabel(xlabel,fontproperties=font_prop)
	ax.set_xticks(np.round(bin_edges,2))
	#ax.xaxis.set_major_locator(plt.MaxNLocator(5))
	#ax.xaxis.set_major_formatter('{x:9<3.1f}')
	
	ax.set_ylabel(ylabel,fontproperties=font_prop)

	fig.suptitle(title)	

	#ax.text(0.8,
	#		0.9,
	#		"N = "+str(len(param_DATA)),
	#		bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
	#		boxstyle='square,pad=1'),
	#		transform=ax.transAxes)
	
	ax.grid(True)
	#plt.savefig("/home/ht2059/Documents/Hist_plots/"+str(param))
	
	fig.show()
	




def RF_vs_NBI():
	mode_db = pd.read_csv(savefile)
	ICRH = mode_db["PTOT"]
	NBI = mode_db["NBLM"]
	ICRH = ICRH.dropna()
	NBI = NBI.dropna()
	ICRH = ICRH[NBI.index]
	ICRH = ICRH/1e6 ; NBI = NBI/1e6
	plt.scatter(ICRH,NBI,alpha=0.2)
	plt.xlabel("Mode ICRH power (MW)")
	plt.ylabel("Mode NBI power (MW)")
	plt.show()
	
	
	
	
def amp_vs(DB_path,param,xlabel,ylabel,title,savepath=None,min_f=0,max_f=500):
	
	
	if 1: #Load and filter db based on frequencies
		modedb = pd.read_pickle(DB_path)
		
		boolf = np.zeros(len(modedb),dtype=bool)


		for i in range(len(modedb)):
			f_list = modedb.iloc[i]["F"]
			f_med = f_list[1]
			if f_med >= min_f and f_med <=max_f:
				boolf[i]=True
			else:
				pass
					
		modedb = modedb[boolf]
		
	
			
	if 1: #Getting and setting data structures
		if param!="prat":
		
			params=np.zeros([len(modedb),3])
			params_data = modedb[param]
			for i in range(len(modedb)):
				params[i]=params_data.iloc[i]
				
			params = params[:,2]
			params = np.array(params)
			boolnan = np.isnan(params)
			params = params[~boolnan]
			
		else:
			ICRH = np.zeros([len(modedb),3]); NBI = ([len(modedb),3])
			ICRH_data = modedb["PTOT"]; NBI_data = modedb["NBLM"]
			for i in range(len(modedb)):
				row = modedb.iloc[i]
				ICRH_r = row[ICRH];NBI_r = row["NBLM"]
				ICRH_data[i]=ICRH_r;NBI_data[i]=NBI_r
				
			ICRH = ICRH[:,2]; NBI = NBI[:,2]
			ICRH = np.array(ICRH); NBI = np.array(NBI)

			prat = ICRH/NBI
			prat = prat[~np.isnan(prat)];prat = prat[~np.isinf(prat)]
			params = prat
			param = "ICRH/NBI"
		
		AMPL = modedb["A"]
		AMPL = AMPL[~boolnan]

		nfactor = np.amax(np.abs(params))
		nfactor = 1.5e20
		
		AMPL_norm = AMPL/np.amax(AMPL)
		params_norm=params/nfactor
		
		
		
		
		
		
		
		
		
	
	if 1: #Linear fit
		m,c= poly.polyfit(params_norm,AMPL_norm,1)
		
		xr = np.linspace(np.amin(params),np.amax(params),100)
		yr = m*xr+c
	
	
	if 1: #plot
		plt.clf()
		plt.scatter(params,AMPL,alpha=0.5,s=5)
		
		plt.errorbar(params,
					AMPL,
					yerr = N/A**0.5,
					marker = '.',
					drawstyle="steps-mid",
					elinewidth=1.5,
					capsize=3,
					ecolor=(0,0,0),
					barsabove=True,
					)
		
		plt.plot(xr,yr,"r")
		
		plt.xlabel(xlabel)
		#plt.xscale("log")
		#plt.xlim([10e3,10e6])
		
		
		plt.ylabel(ylabel)
		plt.yscale("log")
		plt.ylim([10e-5,10e0])
		

		
		plt.title(title)
		plt.grid(True)
	
	
	if 1: #Save or show
		if savepath!=None:
			plt.savefig(savepath+"/"+str(param))
		else:
			plt.show()














def t_vs(DB_path,param,savepath=None,min_f=0,max_f=500):
	
	
	if 1: #Load and filter db based on frequencies
		modedb = pd.read_pickle(DB_path)
		
		boolf = np.zeros(len(modedb),dtype=bool)


		for i in range(len(modedb)):
			f_list = modedb.iloc[i]["F"]
			f_med = f_list[1]
			if f_med >= min_f and f_med <=max_f:
				boolf[i]=True
			else:
				pass
					
		modedb = modedb[boolf]
	
	
	
	if 1 : #Getting and filtering data
		time_arrs = modedb["T"]
		times = np.zeros([len(modedb),5])
		for i in range(len(modedb)):
			times[i]=time_arrs.iloc[i]
			
			
			#tmodestart,tmodemed,tmodeend,tStart
		T_START = times[:,0]; 
		T_STOP = times[:,4]; #Change 2 to 4 for t_A, time to reach max amplitude
		T_MODE = T_STOP
		
		params=np.zeros([len(modedb),3])
		params_data = modedb[param]
		for i in range(len(modedb)):
			params[i]=params_data.iloc[i]
			
			
		params = params[:,2]
		params = np.array(params)
		boolnan = np.isnan(params)
		params = params[~boolnan]
		
		#inds = np.array()

		T_MODE = T_MODE[~boolnan]
		
		inds = T_MODE < 75e-3
		T_MODE = T_MODE[~inds]
		params = params[~inds]
	
	
	
	
	if 1: #Linear fit
		params_fit = params/np.amax(np.abs(params))
		T_FIT = T_MODE/np.amax(np.abs(T_MODE))
		m,c= poly.polyfit(params_fit,T_FIT,1)
		
		xr = np.linspace(np.amin(params),np.amax(params),len(T_MODE))
		yr = m*xr+c
		
		
		
	if 1: #plot
		plt.clf()
		plt.scatter(params,T_MODE,alpha=0.5,s=20)
		plt.plot(xr,yr,"r","-")
		plt.ylabel("Mode time length (s)")
		plt.yscale("log")
		plt.ylim([10e-4,2*10e1])
		
		plt.xlabel(param)
		
		
		plt.title("Identified modes amplitude {0}kHz - {1}kHz {2} distribution | N= {3}".format(min_f,max_f,param,len(params)))

		plt.grid(True)
		
	return m,c,xr,yr,params,T_MODE,params_fit,T_FIT



	if 1: #Save or show
		if savepath!=None:
			plt.savefig(savepath+"/"+str(param))
		else:
			plt.show()
	















def AvsB(DB_path,A,B,xlabel=None,ylabel=None,title=None,xlog=False,ylog=False):
	mode_db = pd.read_pickle(DB_path)
	A = mode_db[A]; A = A.dropna()
	B = mode_db[B]; B = B.dropna()
	A = A.loc[~(A==0)]
	B = B.loc[~(B==0)]
	
	
	low = np.min([len(A),len(B)])
	if low==len(A):
		B=B[A.index]
		N = len(A)
	else:
		A=A[B.index]
		N = len(B)
	
	plt.clf()
	
	plt.scatter(A,B,alpha=0.3)
	
	if xlog==True:
		plt.xscale("log")
		
	if ylog==True:
		plt.yscale("log")
		
	if xlabel!=None:
		plt.xlabel(xlabel)
		
	if ylabel!=None:
		plt.ylabel(ylabel)
		
	if title!=None:
		plt.title(title)
		
	plt.grid(True)
	plt.xlim([np.min(A),np.max(A)])
	plt.ylim([np.min(B),np.max(B)])
	plt.text(0.8,0.1,"N= "+str(N))
	plt.show()
	














def find_prat(DB_path):
	#Find ICRH/NBI power for all modes identified
	mode_db = pd.read_pickle(DB_path)
	ICRH = mode_db["PTOT"]; NBI = mode_db["NBLM"]
	ICRH = mode_db
	ICRH = ICRH.dropna();NBI = NBI.dropna()
	ICRH = ICRH.loc[~(ICRH==0)]; NBI = NBI.loc[~(NBI==0)];
	ICRH = ICRH[NBI.index]
	param_data = ICRH/NBI
	
	return param_data
"""
def n_vs(DB_path,param,ylabel=None,title=None):
	mode_db = pd.read_pickle(DB_path)
	ns=[]
	empty=[]
	for i in mode_db.index:
		try:
			row = mode_db.iloc[i]
			n_arr = row["n"]
			n = n_arr[0]
			ns.append(n)
		except:
			emtpy.append(i)
	for i in empty:
		
"""



def weighted_q(n,m):
	pmns = m[::2]
	weights = m[1::2]
	try:
		pmns = pmns[:2]
		weights = weights[:2]; norm = np.sum(weights)
		weights_= weights/norm
	except:
		pass

	weighted_prod = np.multiply(pmns,weights)
	est_q = (np.sum(weighted_prod))/(n*len(pmns))
	return est_q
		
def q_surf_dist(DB_path,min_f=0,max_f=100):
	
	modedb = pd.read_pickle(DB_path)
		
	boolf = np.zeros(len(modedb),dtype=bool)


	for i in range(len(modedb)):
		f_list = modedb.iloc[i]["F"]
		f_med = f_list[1]
		if f_med >= min_f and f_med <=max_f:
			boolf[i]=True
		else:
			pass
					
	modedb = modedb[boolf]
	
	ns = modedb["n"]
	ms = modedb["m"]
	
	units = rDB.units()
	

	bins = 30;
	
	q_data=[]
	for i in range(len(modedb)):
		n = ns.iloc[i]
		#print(n)
		try:
			n=n[0] #Take first 
		except:
			continue
		m = ms.iloc[i]
		est_q = weighted_q(n,m)
		q_data.append(est_q)
	q_data = np.array(q_data)
	q_data=q_data[~np.isnan(q_data)]
	q_data=q_data[~np.isinf(q_data)]
	
		
		
	
	plt.clf()
	N, bin_edges,patches = plt.hist(q_data,bins=bins,log=False)
	bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

	if 1:
		plt.errorbar(
			bin_centers,
			N,
			yerr = N**0.5,
			marker = '.',
			drawstyle="steps-mid",
			elinewidth=0.10,
			ecolor="#682117",
			barsabove=True
		)
	
	

		
	plt.title("Identified modes 0-500khz, q surface distribution")
	plt.text(bin_centers[-4],np.amax(N),"N = "+str(len(q_data)))
	plt.grid(True)
	plt.xlabel("q surface")
	plt.xlim([0,6])
	
	plt.show()
	
	
	
def pulseparam_vs_t(pulse,param):
	MAIN_DB_PATH = "/home/ht2059/Documents/DATABASE_MAIN.h5"
	DB = h5.File(MAIN_DB_PATH,"r")
	sub_db=rDB.find_pulse(pulse)[0]
	pulse_data = DB[sub_db][str(pulse)]
	
	param_data = pulse_data[param]
	param_data = param_data[~np.isnan(param_data)]
	
	if len(param_data)==0:
		print("No data recorded for "+str(param)+" for this pulse and time")
	
	time = pulse_data["Time array (s)"]
	
	return param_data,time
	

def get_data(param,min_f=0,max_f=500):
	
	modedb = pd.read_pickle("/home/ht2059/Documents/mdb.pkl")
	boolf = np.zeros(len(modedb),dtype=bool)
	
	#Setting boolean array and filtering database
	if 1:
		for i in range(len(modedb)):
			f_list = modedb.iloc[i]["F"]
			f_med = f_list[1]
			if f_med >= min_f and f_med <=max_f:
				boolf[i]=True
			else:
				pass
				
	modedb = modedb[boolf] 

	modedb_param = modedb[param]
	param_data = np.zeros([len(modedb),10])
	
	for i in range(len(modedb)):
		
		param_data_r=np.asarray(modedb_param.iloc[i])
		
		try:
			pad = np.zeros(10-len(param_data_r))
			hold = np.hstack((param_data_r,pad))
			param_data[i]=hold
		except:
			param_data[i,0]=param_data_r
		
	
	param_data = param_data[:,~np.all(param_data==0,axis=0)]
	
	return param_data
		


def sc_hist(p1,p2):
	fig 	= plt.figure();
	alpha 	= 0.5;
	bins 	= 20;

	scatter_axes 	= plt.subplot2grid((3,3),(1,0),rowspan=2,colspan=2);
	x_hist_axes 	= plt.subplot2grid((3,3),(0,0),colspan=2,sharex=scatter_axes);
	y_hist_axes 	= plt.subplot2grid((3,3),(1,2),rowspan=2,sharey=scatter_axes);

	""""
	x1 = dIpdt_Ip_C; 	x2 = dIpdt_Ip_ILW; 		xbins = np.linspace(0,250,11);		xlabel = 'dIp/dt/Ip (1/s)';
	x1 = Btor_C; 		x2 = Btor_ILW; 			xbins = np.linspace(1,4,11); 		xlabel = 'Btor (T)';
	x1 = Ip_C*1e-6;		x2 = Ip_ILW*1e-6;		xbins = np.arange(0.75,7,0.5); 		xlabel = 'Ip (MA)';
	x1 = q95_C*Ip_C*1e-6/Btor_C;x2=q95_ILW*Ip_ILW*1e-6/Btor_ILW;xbins=np.linspace(1,6,11);xlabel='q95 Ip/Bt (MA/T)';
	#x1 = q95_C; 		x2 = q95_ILW; 			xbins = np.linspace(2,10,17); 		xlabel = 'q95';
	#x1 = idx_C; 		x2 = idx_ILW; 			xbins = np.arange(0.75,3.75,0.5); 	xlabel = 'outcome';
	#y1 = Btor_C; 		y2 = Btor_ILW; 			ybins = np.linspace(1,4,11); 		ylabel = 'Btor (T)';
	y1 = Ir_Ip_C; 	 	y2 = Ir_Ip_ILW; 		ybins = np.linspace(0,1,21);		ylabel = 'Ir/Ip';
	#y1 = idx_C; 		y2 = idx_ILW; 			ybins = np.arange(0.75,3.75,0.5); 	ylabel = 'outcome';"""
	
	x = get_data(p1) 
	y = get_data(p2)
	
	
	
	#If x is param, y is A
	x = x[:,2]
	y = y[:,0]
	xbins = np.linspace(np.amin(x[~np.isnan(x)]),np.amax(x[~np.isnan(x)]),20)
	ybins = np.linspace(np.amin(y[~np.isnan(y)]),np.amax(y[~np.isnan(y)]),20)


	scatter_axes.plot(x,y,'o',alpha=alpha,label='C');
	x_hist_axes.hist(x[~np.isnan(x)],bins=xbins,alpha=alpha);
	y_hist_axes.hist(y[~np.isnan(y)],bins=ybins,alpha=alpha,orientation='horizontal');

	#scatter_axes.plot(x2,y2,'s',alpha=alpha,label='ILW');
	#x_hist_axes.hist(x2[~np.isnan(x2)],bins=xbins,alpha=alpha,normed=True);
	#y_hist_axes.hist(y2[~np.isnan(y2)],bins=ybins,alpha=alpha,normed=True,orientation='horizontal');

	scatter_axes.set_xlim((min(xbins),max(xbins)));
	scatter_axes.set_ylim((min(ybins),max(ybins)));
	scatter_axes.set_xlabel(p1);
	scatter_axes.set_ylabel(p2);



	fig.tight_layout();
	fig.show();

	if 0: ### Save fig ###
		filename = ylabel.replace('/','_')+'_v_'+xlabel.replace('/','_');
		fig.savefig(filename+'.png',dpi=250);
