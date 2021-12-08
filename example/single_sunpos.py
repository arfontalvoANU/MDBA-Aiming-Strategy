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

	azimuth=270.0 - Model.phi
	if (azimuth>=360.0 or azimuth<0.0):
		azimuth = (azimuth+360.0)%(360.0)

	# creating a case folder for each new simul file
	designfolder = '%s/vtk'%basefolder
	casefolder = '%s/pos'%basefolder
	if not os.path.exists(casefolder):
		os.makedirs(casefolder)

	runSOLSTICE(
		azimuth,
		Model.elevation,
		designfolder,
		casefolder,
		num_rays
		)

	# Optical postprocessing
	eta,q_results,eta_exc_intec=proces_raw_results(
			'%s/simul'%casefolder,
			casefolder)
	eff_interception=eta/eta_exc_intec
	print 'Optical efficiency: %s'%(eta)

	# Read flux map
	read_data(
			casefolder,
			Model.r_height,
			Model.r_diameter,
			Model.num_bundle,
			Model.bins,
			flux_file=True
			)

	Model.simple_HT_model(20.0,0.0,casefolder)

