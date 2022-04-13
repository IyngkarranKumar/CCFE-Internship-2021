"""
The purpose of the following script is to define some functions for the visualisation of the processed data contained in 
Database I and Database II. 


Andrew Tafurkumar - 20/08/21
"""

if 1: ##IMPORT
	
	import matplotlib.pyplot as plt
	import numpy as np
	import pandas as pd
	from copy import deepcopy
	import h5py as h5
	
	import readUAEDB as rDB
	
	


DATABASE_I_path = "/home/ht2059/Documents/Final Scripts and Databases/DATABASE_MAIN.h5"
DATABASE_II_path = "/home/ht2059/Documents/Final Scripts and Databases/Identified_Unstable_modes.pkl"

plt.ioff()

def frequency_distribution():
	
	xlabel = "Frequency (kHz)"
	ylabel = ""
	title  = "Mode frequency distribution"
	titlefontsize=20
	labelfontsize=15
	ticksize=20
	
	
	
	
	if 1: #GET DATA
		modedb = pd.read_pickle(DATABASE_II_path)
		F = modedb["F"]
		
		param_data = []
		for i in range(len(modedb)):
			row = F.iloc[i]
			f_med = row[1] #Median freq
			f_med = f_med.astype("float64")
			param_data.append(f_med)
			
		F = np.array(F)
		F = param_data
	
	
	
	if 1: ##PLOT
		def colour(c): #Colour 
			c = np.array(c); c = c/255; c = tuple(c)
			return c
		
		fig,ax = plt.subplots()	
		
		
		if 1: #HISTOGRAM PLOT
			bins = 30		
			A,__,__ = ax.hist(param_data,bins=bins,log=True);ax.clear()
			N, bin_edges,patches = ax.hist(param_data,bins=bins,
											log=True,
											density=True,
											color=colour((102,0,128)),
											edgecolor="black",
											alpha=0.8)
			
			bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
			
			A = np.array(A); N = np.array(N)
		
			if 1: #ERRORBARS
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
				
				
			bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])	
		
		
		if 1: #AESTHETICS

			#X
			ax.set_xlabel(xlabel,fontsize=labelfontsize)
			ax.tick_params(axis="x", labelsize=ticksize)
			plt.axvline(x=100,color="r")
			plt.axvline(x=300,color="r")
		
			
			
			
			#Y
			ax.set_ylabel(ylabel,fontsize=labelfontsize)
			ax.set_ylim([0,10])
			ax.tick_params(axis="y", labelsize=ticksize)

			
			
			
			#OTHER	
			text = "N = {0}".format(len(param_data))
			ax.text(0.8,
					0.7,
					text,
					bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
					boxstyle='square,pad=1'),
					transform=ax.transAxes)
			ax.grid(True)
			fig.suptitle(title,fontsize=titlefontsize)
			
			
		#SAVE/SHOW	
		fig.show()
	
	
	
	
	
def parameter_distribution():
	
	
	if 1: ##MAIN SETUP
		fig,axs = plt.subplots(3,1,sharex=True)
		param="TE0"; scalingf = 1e-3 #Use scalingf if want to change units (e.g. Use scalingf=1e-6 to plot ICRH power in MW) 
		xlabel = r"Central electron temperature , $T_{e0}, (keV)$"
		ylabel = ""
		title  = r""
		titlefontsize=20
		labelfontsize=15
		ticksize=20
	
	
	
	
	#MAIN LOOP
	for i in range(len(axs)):
		ax = axs[i]
		
		if i==0: ##DATABASE I GET DATA
			
			param_data=[]
			
			j=0
				
			DATABASE_I = h5.File(DATABASE_I_path,"r")
			
			print("\nGETTING DATA..")
			
			pulses = rDB.pulse_list()
			
			for pulse in pulses:
				try:
					sub_db = rDB.find_pulse(pulse)[0]
				except:
					print(str(pulse)+" not found")
					continue
				
				if "TMN" in param or "PMN" in param:
					data = DATABASE_I[sub_db][str(pulse)][param]
					data = data[::10]
					
				else:
					data = DATABASE_I[sub_db][str(pulse)][param]
					

				data = data[~np.isnan(data)];
				data = data[~np.isinf(data)];
				
				data = list(data)
				param_data= param_data+data

		
		else: ##DATABASE II GET DATA
			num_modes=[]
			param_data = []

			if i==1:
				f_min = 0
				f_max = 100
				
			elif i==2:
				f_min = 100
				f_max = 300
				
			modedb = pd.read_pickle(DATABASE_II_path)
			
			
			for j in range(len(modedb)):
				row = modedb.iloc[j]
				frange = row["F"]
				
				if frange[0]>=f_min and frange[2]<f_max:
					pass
				else:
					continue #skip
				
				par = row[param]
				par = par[2]
				param_data.append(par)
				
				
			param_data = np.array(param_data)
			param_data=param_data[~np.isnan(param_data)]
			param_data = param_data[~np.isinf(param_data)]
			
			
		if 1: ##HISTOGRAM PLOT
			param_data = np.array(param_data)*scalingf
			
			if i==0:
				minr = np.amin(param_data)
				maxr = np.amax(param_data)
				bins = np.linspace(minr,maxr,30)
			def colour(c): #Colour
				c = np.array(c); c = c/255; c = tuple(c)
				return c
					
			A,__,__ = ax.hist(param_data,bins=bins,log=True);ax.clear()
			N, bin_edges,patches = ax.hist(param_data,bins=bins,
											log=False,
											density=True,
											color=colour((102,0,128)),
											edgecolor="black",
											alpha=0.8)
			
			
			bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
			
			A = np.array(A); N = np.array(N)
		
			if 1: #ERROR BARS
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
				
			
		if 1: ##AESTHETICS
			
			
			if 1: #AXES SPECIFIC
				
				if i==0:
					axtitle = "DATABASE I distribution"
				elif i==1:
					axtitle = "DATABASE II 0-100 kHz distribution"
				elif i==2:
					axtitle = "DATABASE II 100-300 kHz distribution"
			
			#X
			if i==2:
				ax.set_xlabel(xlabel,fontsize=labelfontsize)
				ax.tick_params(axis="x", labelsize=ticksize)

		
			
			
			
			#Y
			ax.set_ylabel(ylabel,fontsize=labelfontsize)
			ax.set_ylim([0,3])
			ax.tick_params(axis="y", labelsize=ticksize)

			
			
			
			#OTHER	
			text = "N = {0}".format(len(param_data))
			ax.text(0.8,
					0.7,
					text,
					bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
					boxstyle='square,pad=1'),
					transform=ax.transAxes)
			ax.grid(True)
			ax.set_title(axtitle)
			
			
	fig.suptitle(title,fontsize=titlefontsize)
	
	#SAVE/SHOW
	
	fig.show()
			
	
	
	
	
def amplitude_scatter():
	
	
	if 1: ## MAIN SETUP
		param = "btnm"
		xlabel = r"Normalized plasma beta, $\beta_{N}$"
		ylabel = r"Mode amplitude, $Ts^{-1}$"
		title  = r"Mode amplitude vs $\beta_{N}$"
		titlefontsize=20
		labelfontsize=15
		ticksize=20
	
	
	if 1: ##GET DATA
		
		modedb = pd.read_pickle(DATABASE_II_path)
		
		amp_low=[]
		param_low=[]
		amp_interm=[]
		param_interm=[]
		amp_high=[]
		param_high=[]
		
		fig,ax = plt.subplots()
		
		for i in range(len(modedb)):
			row = modedb.iloc[i]
			flist = row["F"]
			
			amp = row["A"]
				
			param_list = row[param]
			p = param_list[2]
			#print(param_list)
			#print("param is"+str(p))
			
			if np.isnan(p)==True:
				#print(param)
				continue
			
			
			if flist[1]>0 and flist[1]<100:
				amp_low.append(amp)
				param_low.append(p)
			
			elif flist[1]>100 and flist[1]<300:
				amp_interm.append(amp)
				param_interm.append(p)
				
			elif flist[1]>300 and flist[1]<500:
				amp_high.append(amp)
				param_high.append(p)
				
				
			else:
				continue
				
	
	if 1: ##PLOT
		
		"""
		if 0: #Linear regression
			amp_whole = amp_low+amp_interm
			param_whole = param_low+param_interm
			amp_whole_norm = amp_whole/np.amax(amp_whole)
			param_whole_norm = param_whole/np.amax(param_whole)
			m,c = poly.polyfit(param_whole_norm,amp_whole_norm,1)
			xr = np.linspace(np.amin(param_whole),np.amax(param_whole),1000)
			yr = m*xr+c
		"""
		
		def colour(c):
			c = np.array(c); c = c/255; c = tuple(c)
			return c
			

		if 1: ##Scatter
		
			ax.scatter(param_low,amp_low,color="blue",label="0-100kHz band",s=10,alpha=0.5)
			ax.scatter(param_interm,amp_interm,color="red",label="100-300kHz band",s=10,alpha=0.3)
			#ax.scatter(param_high,amp_high,color="green",label="300-500kHz band",s=200,alpha=0.5)
			#ax.plot(xr,yr,color="green",label = "Linear regression")
			
			#ax.scatter(amp_low,param_low,color="blue",label="0-100kHz band",s=10,alpha=0.5)
			#ax.scatter(amp_interm,param_interm,color="blue",label="100-300kHz band",s=10,alpha=0.5)

			ax.set_xlabel(xlabel,fontsize=labelfontsize)
			ax.tick_params(axis="x", labelsize=ticksize)
			#ax.set_xscale("log")
			ax.set_xlim([0,2.5])

			
			ax.set_ylabel(ylabel,fontsize=labelfontsize)
			ax.set_yscale("log")
			ax.tick_params(axis="y", labelsize=ticksize)
			ax.set_ylim([10e-5,10e-1])
			
			"""
			ax.text(0.8,
				0.8,
				text,
				bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
				boxstyle='square,pad=1'),
				transform=ax.transAxes)
			"""			
			
			ax.legend(fontsize=12)
			ax.grid()
	
	fig.suptitle(title,fontsize=titlefontsize)
	

	
def ka3_scatter():
	
	if 1: ## MAIN SETUP
		param = "ka3_ccds"
		xlabel = r"Mode amplitude ($Ts^{-1}$)"
		ylabel = r"Fast ion losses"
		title  = r"Amplitude vs FI losses"
		titlefontsize=20
		labelfontsize=15
		ticksize=15

	
	if 1: ##GET DATA
		
		modedb = pd.read_pickle(DATABASE_II_path)
		
		amp_low=[]
		param_low=[]
		amp_interm=[]
		param_interm=[]
		amp_high=[]
		param_high=[]
		
		fig,ax = plt.subplots()
		
		for i in range(len(modedb)):
			row = modedb.iloc[i]
			flist = row["F"]
			
			amp = row["A"]
				
			param_list = row[param]
			p = param_list[2]
			#print(param_list)
			#print("param is"+str(p))
			
			if np.isnan(p)==True:
				#print(param)
				continue
			
			
			if flist[1]>0 and flist[1]<100:
				amp_low.append(amp)
				param_low.append(p)
			
			elif flist[1]>100 and flist[1]<300:
				amp_interm.append(amp)
				param_interm.append(p)
				
			elif flist[1]>300 and flist[1]<500:
				amp_high.append(amp)
				param_high.append(p)
				
				
			else:
				continue
				
	
	if 1: ##PLOT
		
		"""
		if 0: #Linear regression
			amp_whole = amp_low+amp_interm
			param_whole = param_low+param_interm
			amp_whole_norm = amp_whole/np.amax(amp_whole)
			param_whole_norm = param_whole/np.amax(param_whole)
			m,c = poly.polyfit(param_whole_norm,amp_whole_norm,1)
			xr = np.linspace(np.amin(param_whole),np.amax(param_whole),1000)
			yr = m*xr+c
		"""
		
		def colour(c):
			c = np.array(c); c = c/255; c = tuple(c)
			return c
			

		if 1: ##Scatter
		
			#ax.scatter(param_low,amp_low,color="blue",label="0-100kHz band",s=10,alpha=0.5)
			#ax.scatter(param_interm,amp_interm,color="red",label="100-300kHz band",s=10,alpha=0.3)
			#ax.scatter(param_high,amp_high,color="green",label="300-500kHz band",s=200,alpha=0.5)
			#ax.plot(xr,yr,color="green",label = "Linear regression")
			
			ax.scatter(amp_low,param_low,color="blue",label="0-100kHz band",s=10,alpha=0.5)
			ax.scatter(amp_interm,param_interm,color="red",label="100-300kHz band",s=10,alpha=0.2)

			ax.set_xlabel(xlabel,fontsize=labelfontsize)
			ax.tick_params(axis="x", labelsize=ticksize)
			ax.set_xscale("log")
			ax.set_xlim([10e-5,10e-1])

			
			ax.set_ylabel(ylabel,fontsize=labelfontsize)
			ax.set_yscale("log")
			ax.tick_params(axis="y", labelsize=ticksize)
			ax.set_ylim([10e1,10e5])
			
			"""
			ax.text(0.8,
				0.8,
				text,
				bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
				boxstyle='square,pad=1'),
				transform=ax.transAxes)
			"""			
			
			ax.legend(fontsize=12)
			ax.grid()
	
	fig.suptitle(title,fontsize=titlefontsize)
	
def ndist():
	
	if 1: ## MAIN SETUP
		xlabel = r"Toroidal mode number, $n$"
		ylabel = r""
		title  = r"Toroidal mode number distribution"
		titlefontsize=20
		labelfontsize=15
		ticksize=20
	
	if 1: ##MAIN LOOP
		if 1: ##GET DATA
			fig,axs = plt.subplots(2,1,sharex=True)
		
			for i in range(2):
				ax = axs[i]
				
				if i==0:
					f_min = 0
					f_max = 100
						
				else:
					f_min = 100
					f_max = 300
					

				modedb = pd.read_pickle(DATABASE_II_path)

				
				param_data=[]
				for j in range(len(modedb)):

					
					row = modedb.iloc[j]
					frange = row["F"]
					
					if frange[0]>=f_min and frange[2]<f_max:
						pass
					else:
						continue
						
					
					par = row["n"]
					
					try:
						par = par[0] #Get n
						
						param_data.append(par)
					except:
						pass

				
				

				param_data = np.array(param_data)
				print(len(param_data))




				if 1: ##PLOT
					bins = np.arange(np.amin(param_data)-0.5,np.amax(param_data)+1.5,1)
					def colour(c):
						c = np.array(c); c = c/255; c = tuple(c)
						return c
							
					A,__,__ = ax.hist(param_data,bins=bins,log=True);ax.clear()
					N, bin_edges,patches = ax.hist(param_data,bins=bins,
													log=True,
													density=True,
													color=colour((102,0,128)),
													edgecolor="black",
													alpha=0.8)
					
					
					bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
					
					mean = np.average(bin_centers,weights=N)
					print(mean)
					
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
					
					
					
					
					
					

					#ax.set_ylim([0,0.1])


				
				
				if 1: #AXES SPECIFIC

					if i==0:
						axtitle = "DATABASE II 0-100 kHz distribution"
					else:
						axtitle = "DATABASE II 100-300 kHz distribution"
				
					#X
					if i==1:
						ax.set_xlabel(xlabel,fontsize=labelfontsize)
						ticks = np.arange(-1,9,1)
						ax.set_xticks(ticks)
						ax.tick_params(axis="x", labelsize=ticksize)

				
					
					
					
					#Y
					ax.set_ylabel(ylabel,fontsize=labelfontsize)
					ax.set_ylim([0,5])
					ax.tick_params(axis="y", labelsize=ticksize)

					
					
					
					#OTHER	
					text = "N = {0}".format(len(param_data))
					ax.text(0.8,
							0.7,
							text,
							bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
							boxstyle='square,pad=1'),
							transform=ax.transAxes)
					ax.grid(True)
					ax.set_title(axtitle)
			
			
	fig.suptitle(title,fontsize=titlefontsize)


def qdist():
	
	if 1: ## MAIN SETUP
		xlabel = r"Estimated q (=m/n) surface"
		ylabel = r""
		title  = r"Mode distribution across q (=m/n) surfaces"
		titlefontsize=20
		labelfontsize=15
		ticksize=20
		
	if 1: ##MAIN LOOP
		if 1: ##GET DATA
			fig,axs = plt.subplots(2,1,sharex=True)
		
			for i in range(2):
				ax = axs[i]
				
				if i==0:
					f_min = 0
					f_max = 100
						
				else:
					f_min = 100
					f_max = 300
					

				modedb = pd.read_pickle(DATABASE_II_path)

				
				param_data=[]
				for j in range(len(modedb)):

					
					row = modedb.iloc[j]
					frange = row["F"]
					
					if frange[0]>=f_min and frange[2]<f_max:
						pass
					else:
						continue
						
					
					nrow = row["n"]
					
					if len(nrow)!=2:
						continue
					
					try:
						n=nrow[0]
					except:
						continue
					
					mrow = row["m"]
					
					if len(mrow)==0:
						continue
					elif len(mrow)>4:
						continue
					else:
						ms = np.array(mrow[::2])
						
						for m in ms:
							q = m/n
							
							if q<1 or q>3:
								continue
								
							param_data.append(q)

						

			
				
				
			
				param_data = np.array(param_data)
				param_data = param_data[~np.isnan(param_data)]
				param_data = param_data[~np.isinf(param_data)]





				if 1: ##PLOT
					bins = 20
					bins = np.arange(0.9,3.2,0.2)
					def colour(c):
						c = np.array(c); c = c/255; c = tuple(c)
						return c
							
					A,__,__ = ax.hist(param_data,bins=bins,log=True);ax.clear()
					N, bin_edges,patches = ax.hist(param_data,bins=bins,
													log=True,
													density=True,
													color=colour((102,0,128)),
													edgecolor="black",
													alpha=0.8)
					
					
					bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
					
					mean = np.average(bin_centers,weights=N)
					print(mean)
					
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
					
					
					
					
					
					

					#ax.set_ylim([0,0.1])


				
				
				if 1: #AXES SPECIFIC

					if i==0:
						axtitle = "DATABASE II 0-100 kHz distribution"
					else:
						axtitle = "DATABASE II 100-300 kHz distribution"
				
					#X
					if i==1:
						ax.set_xlabel(xlabel,fontsize=labelfontsize)
						#ticks = bins
						#ax.set_xticks(bins)
						ax.tick_params(axis="x", labelsize=ticksize)
						print(np.amax(param_data),np.amin(param_data))

				
					
					
					
					#Y
					ax.set_ylabel(ylabel,fontsize=labelfontsize)
					ax.set_ylim([0,5])
					ax.tick_params(axis="y", labelsize=ticksize)

					
					
					print(bins,bin_edges,bin_centers)
					#OTHER	
					text = "N = {0}".format(len(param_data))
					ax.text(0.8,
							0.7,
							text,
							bbox=dict(facecolor=colour((204, 204, 255)), alpha=0.5,
							boxstyle='square,pad=1'),
							transform=ax.transAxes)
					ax.grid(True)
					ax.set_title(axtitle)
			
			
	fig.suptitle(title,fontsize=titlefontsize)
		
