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

import onekey_aiming as aiming
from simulation import optical
from run_solstice import *
from python_postprocessing import *
from flux_reader import *

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def runSOLSTICE(azimuth,elevation,yaml_folder,simul_folder,num_rays):
	'''
	The scrip to run SOLSTICE with the user-defined sun positions:
	Inputs:
		azimuth:                  Solar azimuth
		elevation:                Solar elevation
		yaml_folder:              Directory to the YAML input for SOLSTICE
		simul_folder:             Directory to save the simul file
		num_rays:                 Number of rays
	Outputs:
		simul:                    SOLSTICE text file saved in simul_folder
	'''

	print(yellow('solstice -D%s,%s -v -n %s -R %s/demo-rcv.yaml -fo %s/simul %s/demo.yaml'%(
		azimuth,
		elevation,
		int(num_rays),
		yaml_folder,
		simul_folder,
		yaml_folder
		))
		)
	os.system('solstice -D%s,%s -v -n %s -R %s/demo-rcv.yaml -fo %s/simul %s/demo.yaml'%(
		azimuth,
		elevation,
		int(num_rays),
		yaml_folder,
		simul_folder,
		yaml_folder
		)
		)

if __name__=='__main__':
	# define a unique case folder for the user
	snum = 0
	suffix = ""
	while 1:
		dt = datetime.datetime.now()
		ds = dt.strftime("%a-%H-%M")
		basefolder = os.path.join(currentdir,'case-%s%s'%(ds,suffix))
		if os.path.exists(basefolder):
			snum+=1
			suffix = "-%d"%(snum,)
			if snum > 200:
				raise RuntimeError("Some problem with creating basefolder")
		else:
			# good, we have a new case dir
			break

	if not os.path.exists(basefolder):
		os.makedirs(basefolder)

	folder=basefolder
	os.system("cp pos_and_aiming.csv %s"%(basefolder))
	source_path=parentdir

	# Inputs
	num_bundle=16         # Number of panels
	r_height=24.0         # Receiver height [m]
	r_diameter=16.0       # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=175.0         # Tower height [m]
	azimuth=0.0           # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=980.0             # Beam irradiance [W/m2]
	num_fp=8              # Two panels per flow path (num_bundle/2)
	D0=60.33              # Panel tube OD [mm]
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
	Model.New_search_algorithm()

	# Creating RLLT text file
	f = open('%s/RLLT.csv'%basefolder,'a')
	header = \
		'h_conv_ext (W/m2/K), T_avg (K), (T**4_avg)**0.25 (K),'+\
		'conv_loss (W), rad_loss (W), receiver_loss (W),'+\
		'receiver_output (W), eta_th\n'
	f.write(header)

	# creating a case folder for each new simul file
	designfolder = '%s/vtk'%basefolder
	casefolder = '%s/pos'%basefolder
	if not os.path.exists(casefolder):
		os.makedirs(casefolder)

	# Importing RLLT inputs
	data = np.genfromtxt('%s/RLLT_input.csv'%currentdir, delimiter=',', skip_header=1)
	nrows = data.shape[0]
	for i in range(nrows):
		azimuth = data[i,4]
		elevation = data[i,5]
		dni = data[i,6]
		Tamb = data[i,7]
		Wspd = data[i,8]

		# Creating solstice scene
		scene=SolsticeScene(
			mainfolder=basefolder,
			num_rays=Model.num_rays,
			dni=DNI,
			azimuth=azimuth,
			zenith=elevation,
			att_factor=Model.att_factor,
			csv='%s/pos_and_aiming_new.csv'%basefolder,
			tower_h=Model.tower_h,
			r_cyl=Model.r_diameter/2.,
			h_cyl=Model.r_height,
			num_bundle=Model.num_bundle
		)

		# Creating YAML inputs
		scene.gen_YAML()

		# Running solstice with new DNI
		scene.runSOLSTICE(savefile=casefolder,view=True)

		# Optical postprocessing
		eta,q_results,eta_exc_intec=proces_raw_results(
				'%s/simul'%casefolder,
				casefolder)

		# Read flux map
		read_data(
				casefolder,
				Model.r_height,
				Model.r_diameter,
				Model.num_bundle,
				Model.bins,
				flux_file=True
				)

		print(yellow('azi: %s [deg]\tele: %s [deg]\tdni: %s [W/m2]\tTamb: %s [deg C]\tWspd: %s [m/s]'%(azimuth,elevation,dni,Tamb,Wspd)))
		results=Model.simple_HT_model(Tamb,Wspd,casefolder)
		f.write('%s,'%results[5])                       #h_conv_ext
		f.write('%s,'%results[3])                       #T_avg
		f.write('%s,'%results[4])                       #(T**4_avg)**0.25
		f.write('%s,'%(results[8]*1e6))                 #conv_loss
		f.write('%s,'%(results[7]*1e6))                 #rad_loss
		f.write('%s,'%((results[7]+results[8])*1e6))    #receiver_loss
		f.write('%s,'%(results[0]*results[9]*1e6))      #receiver_output
		f.write('%s,'%results[9])                       #receiver_output
		f.write('\n')
	f.close()

