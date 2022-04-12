

Hi Alex. Here's a quick overview of this folder, which contains the final work from the internship. 
As mentioned further on, I recommend checking out test_environment.py to see how some of the scripts function



UAE_DB1_COMPILE - A script for the compilation of 'Database I' (4 spectrograms & 31 parameter time arrays). 

MEA.py - Mode extraction algorithm code. For experimentation with this, I recommend using the image_save flag to save extraction results (line plots),
then tune the threshold until satisfied with results. 

------------------------------------------------------------------------------------------------------------------------

DB_MAIN - HDF5 file that allows access to the 5 sub Databases that are also stored in this file. 


Interfacing scripts provided by readUAEDB in Auxscripts:
	- Use the "find_pulse(JPN)" function to locate the sub database that a pulse is stored, and index relative to that sub database
	- The spect_reconstruct(JPN) function returns the 4 spectrogram from the nonzero indices stored in DATABASE I
	- pulse_list() and param_list() functions return lists of all pulses stored in Database I, and the 31 numpy arrays associated with each one, respectively

------------------------------------------------------------------------------------------------------------------------

Auxscripts - A number of supporting scripts, mostly, but not all, for Database I compilation. Some important files are:
 
	Plotting script provides functions for visualisation of pulse spectrograms. This is what I used over the course of the internship. 

	Analysis plotting holds some functions for statistical analysis of Databases 

	new_mode_number_calculation is the script used for n & m calculation for Database I spectrograms

	UAEDB_addparameter can be run to add parameter time array to Database I. 

------------------------------------------------------------------------------------------------------------------------

Identified_Unstable_modes.pkl contains Database II (~7179 unstables modes stored as pandas Dataframe)
The structure of the dataframe is as follows:
18 columns. Associated JPN are dataframe indexes

n - toroidal mode number array. If n comprises p% of indices of mode, with p>40, n included in array, as well as p/100. Most only have one n though.
m - Same as above, but threshold to be incorporated is 25 percent only
t - tuple of times. t_start, t_med, t_stop, t_Amax, tStart. t_Amax time to reach max A, tStart is PULSE start time.
f - f_min,f_med,f_max,f_Amax
A - Averaged amplitude of mode

13 parameters arrays make up remaining columns. Each array is 3 elements, parameter averaged over mode time AND 500,250 and 100ms prior to first appearance of mode

------------------------------------------------------------------------------------------------------------------------

test_environment.py is a good place to start to see how some of the scripts above work. 

------------------------------------------------------------------------------------------------------------------------

analysis_plots folder contains the plots compiled over weeks of analysis. 

------------------------------------------------------------------------------------------------------------------------


IMPORTANT!
In all of the scripts I specified Database locations as absolute paths, so functions will not work as expected until these paths are updated.
This must be done in readUAEDB.py, where the DATABASE_I location is defined as a global variable in line 36 
Also defined in testing_environment.py









