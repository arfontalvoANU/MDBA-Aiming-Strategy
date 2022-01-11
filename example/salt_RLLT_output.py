import os
import sys
import inspect
import datetime
from solsticepy.cal_sun import *
import colorama
colorama.init()

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import salt_aiming as aiming
from simulation import optical
from run_solstice import *
from python_postprocessing import *
from flux_reader import *

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

if __name__=='__main__':

#	ratio = 3
#	DNI_ratio=[1.24,1.00,0.76,0.52]

	# define a unique case folder for the user
	basefolder = os.path.join(currentdir,'case-a230salt-RLLT')

	if not os.path.exists(basefolder):
		os.makedirs(basefolder)

	folder=basefolder
	source_path=parentdir

	# Inputs
	num_bundle=18         # Number of panels
	r_height=10.5         # Receiver height [m]
	r_diameter=8.5        # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=114.75        # Tower height [m]
	azimuth=0.0           # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=950.0             # Beam irradiance [W/m2]
	num_fp=2              # Two panels per flow path (num_bundle/2)
	D0=42.16              # Panel tube OD [mm]
	num_rays=int(5e6)     # Number of rays

	# Creating receiver model
	Model=aiming.one_key_start(
		folder,
		source_path,
		num_bundle,
		num_fp,
		r_height,
		r_diameter,
		bins,
		tower_h,
		azimuth,
		elevation,
		DNI,
		D0)
	#Model.New_search_algorithm()

	# Creating RLLT text file
#	f = open('%s/RLLT.csv'%basefolder,'a')
#	header = \
#		'h_conv_ext (W/m2/K), T_avg (K), (T**4_avg)**0.25 (K),'+\
#		'conv_loss (W), rad_loss (W), receiver_loss (W),'+\
#		'receiver_output (W), eta_th\n'
#	f.write(header)

	# Importing RLLT inputs
	data = np.genfromtxt('%s/RLLT_input.csv'%currentdir, delimiter=',', skip_header=1)
	nrows = data.shape[0]
	for i in range(1):
		azimuth = data[i,4]
		elevation = data[i,5]
		dni = data[i,6]
		Tamb = data[i,7]
		Wspd = data[i,8]
		print(yellow('row: %s\tazi: %s [deg]\tele: %s [deg]\tdni: %s [W/m2]\tTamb: %s [deg C]\tWspd: %s [m/s]'%(i,azimuth,elevation,dni,Tamb,Wspd)))
		casefolder = '%s/row_%s'%(basefolder,i)
		if not os.path.exists(casefolder):
			os.makedirs(casefolder)

		# Re-aiming
		Model=aiming.one_key_start(
			casefolder,
			source_path,
			num_bundle,
			num_fp,
			r_height,
			r_diameter,
			bins,
			tower_h,
			azimuth,
			elevation,
			dni,
			D0)
		C_start = 0.5
		E_start = 2.0
		A_start = 0.5
		# Running aiming for design point
		Model.sweeping_algorithm(C_start,E_start,A_start)

		# Optical postprocessing
		eta,q_results,eta_exc_intec=proces_raw_results(
				'%s/vtk/simul'%casefolder,
				casefolder)

		results,aiming_results,Strt=Model.HT_model(Tamb,Wspd)
		print results
#		f.write('%s,'%results[5])                       #h_conv_ext
#		f.write('%s,'%results[3])                       #T_avg
#		f.write('%s,'%results[4])                       #(T**4_avg)**0.25
#		f.write('%s,'%(results[8]*1e6))                 #conv_loss
#		f.write('%s,'%(results[7]*1e6))                 #rad_loss
#		f.write('%s,'%((results[7]+results[8])*1e6))    #receiver_loss
#		f.write('%s,'%(results[0]*results[9]*1e6))      #receiver_output
#		f.write('%s,'%results[9])                       #receiver_output
#		f.write('\n')
#	f.close()

