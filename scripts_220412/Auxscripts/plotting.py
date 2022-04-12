"""
The purpose of this script is to provide functions that create various plots. 
Examples are: 
image - Plots spectogram amplitudes raised to power p
image_n - PLots spectogram but for mode numbers as supposed to Fourier coefficient amplitudes

More will be added 
"""
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import cmath
from random import randint
import h5py as h5
from scipy import io as sio
from numpy.polynomial import polynomial as poly



from Auxscripts import readUAEDB as rDB



Fs = 2e6
dt = 1/Fs
T = 20
N = int(T*Fs)
df = 488



def image(S,T,F,p=0.2,min_f=170,max_f=250,title="Spectogram",cbarlabel=None,cbarlog=False,cmap='plasma_r',savepath=None):

	#Plots/saves spectogram

	#plot image prepared as imagesc from matlab, with limits and all
	extent=[T[0],T[len(T)-1],F[0]/1000,F[len(F)-1]/1000]   #x axis, y axis min and max
	if cbarlog==False:
		norm= matplotlib.colors.Normalize()
	else:
		norm = matplotlib.colors.LogNorm()
	
	
	im=plt.imshow(np.abs(S)**p,
				  extent=extent,
				  aspect='auto',
				  norm=norm,
				  interpolation='none',
				  origin='lower',
				  cmap=cmap)
	plt.ylim([min_f,max_f])
	plt.ylabel("Frequency[kHz]")
	plt.xlabel("Time[s]")
	plt.title(title)
	cbar=plt.colorbar()
	plt.tight_layout();
	
	if cbarlabel!=None:
		cbar.set_label(cbarlabel)
		
	if savepath==None:
		plt.show()
	
	else:
		plt.savefig(savepath)
		plt.close()
	
		
		
		
		
	


def image_n(S,T,F,NMax,min_f=170,max_f=250,NMin=[],title="Mode number spectogram",cbarlabel="",savepath=None):

	"""Plots/saves mode number spectrogram"""



	#Colour specified by two letters, first is 'main' colour, second is light,dark or medium. E.g: RL is red,light

	RL = np.array([255,179,179])/255;RM = np.array([255,0,0])/255;RD = np.array([128,0,0])/255;
	BL = np.array([179,179,255])/255;BM = np.array([0,0,255])/255;BD = np.array([0,0,128])/255;
	GL = np.array([179,255,179])/255;GM = np.array([0,255,0])/255;GD = np.array([0,128,0])/255;
	YL = np.array([255,255,179])/255;YM = np.array([255,255,0])/255;YD = np.array([128,128,0])/255;
	CL = np.array([179,255,255])/255;CM = np.array([0,255,255])/255;CD = np.array([0,128,128])/255;
	ML = np.array([255,179,255])/255;MM = np.array([255,0,255])/255;MD = np.array([128,50,128])/255;

	"""
	red     = [1,0,0]; green     = [0,1,0];blue     = [0,0,1]; 
	orange     = np.array([255,140,0])/255;purple = np.array([128,0,128])/255;cyan     = [0,1,1]; 
	magenta = [1,0,1]; yellow     = [1,1,0];grey     = [0.5,0.5,0.5];black     = [0,0,0];
	"""

	
	
	
	if type(NMin)==list:
		levels = np.arange(-NMax-0.5,NMax+0.5,1)
	else:
		levels = np.arange(NMin-0.5,NMax+1.5,1)
	
	colors     = np.array([RM,RL,RD,BM,BL,BD,YM,YL,YD,CM,CL,CD,GM,GL,GD,MM,ML,MD])
	#colors     = np.array([RM,BM,YM,GM,CM,MM,RD,BL,YD,GL,CD,ML,RL,BD,YL,GD,CL,MD])
	colors     = np.array([RM,BM,YM,CM,GM,MM,RL,BL,YL,CL,GL,ML,RD,BD,YD,CD,GD,MD]);#BM,BL,BD,YM,YL,YD,CM,CL,CD,GM,GL,GD,MM,ML,MD])
	#colors 	   = np.array([R,C,B,M,G,Y,R,C,B,M,G,Y,R,C,B,M,G,Y]);	
	colors 	   = np.array([RM,CM,BM,MM,GM,YM,RD,CD,BD,MD,GD,YD,RL,CL,BL,ML,GL,YL]);	
	colors 	   = np.array([RM,CD,BL,MM,GD,YL,RD,CL,BM,MD,GL,YM,RL,CM,BD,ML,GM,YD]);	
	#np.random.shuffle(colors)
	colors     = colors[:len(levels)-1]
	
	cmap     = matplotlib.colors.LinearSegmentedColormap.from_list('test',colors,N=len(levels))
	print(levels);print("\n"+str(len(levels)))
	

	

	#plot image prepared as imagesc from matlab, with limits and all
	extent=[T[0],T[len(T)-1],F[0]/1000,F[len(F)-1]/1000]   #x axis, y axis min and max
	#S[0,0]=levels[0];S[0,1]=levels[-1]
	im=plt.imshow(S,
				  extent=extent,
				  aspect='auto',
				  interpolation='none',
				  origin='lower',
				  cmap=cmap)

	plt.ylim([min_f,max_f])
	plt.ylabel("Frequency[kHz]")
	plt.xlabel("Time[s]")
	plt.title(title)
	cbar=plt.colorbar(ticks=levels+0.5)
	cbar.set_label(cbarlabel)
	plt.clim(levels[0],levels[-1]+1)
	plt.tight_layout();
	
	if savepath==None:
		plt.show()
	
	else:
		plt.savefig(savepath)
		plt.close();


def mode_plot(NEST,MEST,T,F,NMax,MMax,NMin=[],MMin=[],min_f=170,max_f=250):
	fig = plt.figure()
	ax1 = fig.add_axes([0.15,0.55,0.6,0.4])
	ax2 = fig.add_axes([0.15,0.03,0.6,0.4])
	fig.tight_layout(h_pad=5)
	
	"""Plots spectogram. Imshow is significantly faster than ax.plot"""



	#Colour specified by two letters, first is 'main' colour, second is light,dark or medium. E.g: RL is red,light

	RL = np.array([255,179,179])/255;RM = np.array([255,0,0])/255;RD = np.array([128,0,0])/255;
	BL = np.array([179,179,255])/255;BM = np.array([0,0,255])/255;BD = np.array([0,0,128])/255;
	GL = np.array([179,255,179])/255;GM = np.array([0,255,0])/255;GD = np.array([0,128,0])/255;
	YL = np.array([255,255,179])/255;YM = np.array([255,255,0])/255;YD = np.array([128,128,0])/255;
	CL = np.array([179,255,255])/255;CM = np.array([0,255,255])/255;CD = np.array([0,128,128])/255;
	ML = np.array([255,179,255])/255;MM = np.array([255,0,255])/255;MD = np.array([79,0,79])/255;
	
	colors     = np.array([RM,RL,RD,CM,CL,CD,YM,YL,YD,GM,GL,GD,MM,ML,MD,BM,BL,BD])
	colors = np.array([MD,BD,YD,GD,RD,CD,RL,BL,YL,GL,ML,CL,RM,BM,GM,MM,YM,CM])

	"""
	red     = [1,0,0]; green     = [0,1,0];blue     = [0,0,1]; 
	orange     = np.array([255,140,0])/255;purple = np.array([128,0,128])/255;cyan     = [0,1,1]; 
	magenta = [1,0,1]; yellow     = [1,1,0];grey     = [0.5,0.5,0.5];black     = [0,0,0];
	"""

	
	#Define maps for toroidal and poloidal plots
	
	if type(NMin)==list:
		nlevels = np.arange(-NMax-0.5,NMax+1.5,1)
	else:
		nlevels = np.arange(NMin-0.5,NMax+1.5,1)
		
		
	if type(MMin)==list:
		mlevels = np.arange(-MMax-0.5,MMax+1.5,1)
	else:
		mlevels = np.arange(MMin-0.5,MMax+1.5,1)
		
	ncolors     = colors[:len(nlevels)+1]
	mcolors     = colors[:len(mlevels)+1]
	
	n_cmap     = matplotlib.colors.LinearSegmentedColormap.from_list('test',ncolors,N=len(nlevels)-1)
	m_cmap     = matplotlib.colors.LinearSegmentedColormap.from_list('test',mcolors,N=len(mlevels)-1)
	#print(levels);print("\n"+str(len(levels)))

	#plot image prepared as imagesc from matlab, with limits and all
	extent=[T[0],T[len(T)-1],F[0]/1000,F[len(F)-1]/1000]   #x axis, y axis min and max. Same for both toroidal and poloidal
	

	im_n=ax1.imshow(NEST,
				  extent=extent,
				  aspect='auto',
				  interpolation='none',
				  origin='lower',
				  cmap=n_cmap,
				  vmax = nlevels[-1],
				  vmin = nlevels[0])
				 

	ax1.set_ylim(min_f,max_f)
	ax1.set_ylabel("Frequency[kHz]")
	ax1.set_xlabel("Time[s]")
	ax1.set_title("Toroidal number spectogram")
	divider = make_axes_locatable(ax1)
	cax1 = divider.new_vertical(size='10%', pad=0.6, pack_start=True)
	fig.add_axes(cax1)
	cb2 = fig.colorbar(im_n,ticks=mlevels+0.5,cax=cax1,orientation="horizontal",aspect=30,shrink=0.7)
	cb2.set_clim(mlevels[0],nlevels[-1])

		
		

	im_m=ax2.imshow(MEST,
				  extent=extent,
				  aspect='auto',
				  interpolation='none',
				  origin='lower',
				  cmap=m_cmap,
				  vmax = mlevels[-1],
				  vmin = mlevels[0])

	ax2.set_ylim(min_f,max_f)
	ax2.set_ylabel("Frequency[kHz]")
	ax2.set_xlabel("Time[s]")
	ax2.set_title("Poloidal number spectogram")
	#cb2 = fig.colorbar(im_m,ticks=mlevels+0.5,ax=ax2,orientation="horizontal")
	#dummy_m = plt.imshow(np.meshgrid(np.linspace(0,6,5),mlevels)[1],cmap=m_cmap)
	divider = make_axes_locatable(ax2)
	cax2 = divider.new_vertical(size='10%', pad=0.6, pack_start=True)
	fig.add_axes(cax2)
	cb2 = fig.colorbar(im_m,ticks=mlevels+0.5,cax=cax2,orientation="horizontal",aspect=30,shrink=0.7)
	cb2.set_clim(mlevels[0],mlevels[-1])


	print(nlevels);print(mlevels)
	plt.show()
	
	
#CREATE NEW FUNCTION. INPUT ARE INDICES FROM SOME MODE NUMBER, OUTPUT IS SPECTOGRAM PLOT
def check(Shape,n_structures,min_f,max_f,tStart,tStop,pulse):
	MODE = np.zeros(shape=Shape)
	MODE[:]=np.nan
	i=0
	for struct in n_structures:
		#i+=1
		#print("\n ")
		#print(i)
		#print(len(str))
		arr = np.array(struct)
		#print(np.amin(arr[:,0]),np.amax(arr[:,0]))
		#print(np.amin(arr[0,:]),np.amax(arr[0,:]))
		for index in struct:
			#index in form (t,f)
			MODE[int(index[1]),int(index[0])]=i
		i+=1
	#fig,ax = plt.subplots()
	extent=[tStart,tStop,0,1e3]   #x axis, y axis min and max
	plt.imshow(MODE,
				  extent=extent,
				  aspect='auto',
				  interpolation='none',
				  origin='lower',
				  cmap="tab10")
	min_fb = min_f*1000//df
	max_fb = max_f*1000//df
	plt.ylim(min_f,max_f)
	plt.ylabel("Frequency (kHz)")
	plt.xlabel("Time (s)")
	plt.title("Mode identification indices JPN "+str(pulse))
	#fig.show()
	#return MODE
	



def line_plot(Shape,struct,min_f,max_f,colour,tStart,tStop,pulse):
	

		
	#print("\n ")
	#print(len(struct))

	arr = np.array(struct)
	
	#print(np.amin(arr[:,0]),np.amax(arr[:,0]))    #Time min and max
	#print(np.amin(arr[:,1]),np.amax(arr[:,1]))    #Freq min and max
	
	colour = np.array((np.random.choice(range(256), size=3)));colour = colour/256
	time = tStart+(arr[:,0]*((tStop-tStart)/Shape[1]))
	freq = arr[:,1]*(1e6/Shape[0])
	freq = freq/1000
	plt.plot(time,freq,color=colour,linewidth=4) #Freq vs time
		
	plt.xlabel("Time (s)")
	plt.ylabel("Frequency (kHz)")
	plt.title("Identified Modes JPN "+str(pulse))
	
	plt.xlim([tStart,tStop])
	min_fb=(min_f//1000)//df; max_fb=(max_f//1000)//df
	plt.ylim([min_f//1000,max_f//1000])
	

		
		
def plot_histogram(param,xlabel,ylabel,title,All=False,pulse_rng=(),freq_range=None):
	
	pulse_data = sio.loadmat("/home/atinguel/work/forIyngkarranK/pulses_w_ICRH_93000-98006_Pthreshold500kW_Tthreshold500ms")
	all_pulses = pulse_data["pulses"][0]
	units = rDB.units()
	
	DB_path = "/home/ht2059/Documents/DATABASE_MAIN.h5"
	DATABASE = h5.File(DB_path,"r")
	
	
	if 1: #Set pulse range
		pulses=[]
		#print(All);print(len(pulse_rng)
		
		
		if len(pulses)==0:
			print("here")
			pulses=all_pulses
		
		
			

		
		
		

			
	print(len(pulses))
	
	if 1:	#Get data
		print("here")
		param_DATA=[]
		
		i=0
		
		DATABASE = h5.File(DB_path,"r")
		
		print("\nGETTING DATA..")
		
		for pulse in pulses:
			#print(pulse)
			try:
				sub_db = rDB.find_pulse(pulse)[0]
			except:
				print(str(pulse)+" not found")
				continue
			
			if "TMN" in param or "PMN" in param:
				data = DATABASE[sub_db][str(pulse)][param]
				
				
					
				data = data[::10]
				print(len(data))
				
			else:
				data = DATABASE[sub_db][str(pulse)][param]

			
			if freq_range!=None:
				if len(freq_range)!=2:
					raise ValueError("Frequency range specified as tuple (fmin,fmax)")
				else:
					freq_indices = np.array(DATABASE[sub_db][str(pulse)]["NZ indices"])
					startf = freq_range[0];stopf=freq_range[-1]
					boolLf = freq_indices < startf
					boolRf = freq_indices > stopf
					boolAllf = np.logical_or(boolLf,boolRf)
					print(len(data));print(len(boolAllf))
					data = data[boolAllf]
					
					
			data = data[~np.isnan(data)];
			data = data[~np.isinf(data)];
			
			data = list(data)
			param_DATA= param_DATA+data

					
			
			
			
		print(len(param_DATA))
		param_DATA = np.array(param_DATA)
		scale = 1e6
		param_DATA = np.round(param_DATA/scale,2)
		
		
		
	if 1: #Plotting histogram

		
		if 'PMN' in param or 'TMN' in param:
			bins = np.arange(min(param_DATA)-0.5,max(param_DATA)+1.5,1);
			
		else:
			bins = 30;
		
		if 1: #Colour
			def colour(c):
				c = np.array(c); c = c/255; c = tuple(c)
				return c
			
			
		
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
		
		A = np.array(A); N = np.array(N)
		
		if 1:
			ax.errorbar(
				bin_centers,
				N,
				yerr = N/A**0.5,
				marker = '.',
				drawstyle="steps-mid",
				elinewidth=1.5,
				capsize=3,
				ecolor=(0,0,0),
				barsabove=True,
			)
		
		
		
		ax.set_xlabel(xlabel)
		ax.set_xticks(bin_edges)
		ax.xaxis.set_major_locator(plt.MaxNLocator(5))
		#ax.xaxis.set_major_formatter('{x:9<3.1f}')
		
		ax.set_ylabel(ylabel)

		fig.suptitle(title)	

		ax.text(0.8,
				0.9,
				"N = "+str(len(param_DATA)),
				bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
				boxstyle='square,pad=1'),
				transform=ax.transAxes)
		
		ax.grid(True)
		#plt.savefig("/home/ht2059/Documents/Hist_plots/"+str(param))
		
		fig.show()
	
#####################################################		





def num_vs(param,DB):
	num_modes=[]
	pulse_param=[]
	main_DB_path = "/home/ht2059/Documents/DATABASE_MAIN.h5"
	MAIN_DB = h5.File(main_DB_path,"r")
	mode_DB_path = "/home/ht2059/Documents/MODE_DATABASE.h5"
	MODE_DB = h5.File(mode_DB_path,"r")
	pulses = rDB.pulse_list_mdb()
	
	for pulse in pulses:
		sub_db = rDB.find_pulse(pulse)[0]
		
		#GET NUMBER OF MODES IDENTIFIED
		pulse_data = MODE_DB[str(pulse)]
		num = len(pulse_data)
		
		#GET AVG VALUE OF PARAMETER OVER PULSE
		sub_db = rDB.find_pulse(pulse)[0]
		param_arr = MAIN_DB[sub_db][str(pulse)][param]
		param_arr = param_arr[~np.isnan(param_arr)]
		param_avg = np.mean(param_arr)
		
		
		if len(param)==0:    #In case that no data returned for whole timebase
			continue
		else:
			num_modes.append(num)
			pulse_param.append(param_avg)
		
		
	num_modes = np.asarray(num_modes)
	pulse_param = np.asarray(pulse_param)
	
	c,m = poly.polyfit(pulse_param,num_modes,1)
	
	xr = np.arange(pulse_param[0],pulse_param[-1],100)
	yr = m*xr+c
	
	plt.plot(pulse_param,num_modes,"bo")
	plt.plot(xr,yr,"r")
	plt.xlabel(param)
	plt.ylabel("Number of modes identified")
	plt.title("Number of modes vs "+str(param))
	plt.show()
	
	
		
	
	
	
	

	
