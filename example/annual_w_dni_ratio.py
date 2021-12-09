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

	ratio = 3
	DNI_ratio=[1.24,1.00,0.76,0.52]

	# define a unique case folder for the user
	snum = 0
	suffix = ""
	while 1:
		dt = datetime.datetime.now()
		ds = dt.strftime("%a-%H-%M")
#		basefolder = os.path.join(currentdir,'case-%s%s'%(ds,suffix))
		basefolder = os.path.join(currentdir,'Haynes230sodium-DNI_ratio_%s'%ratio)
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

	# Annual discretisation
	ndec=5
	nhra=25
	E = np.zeros((ndec, nhra))
	# Declination angle
	Dec=np.linspace(-23.45, 23.45, num=ndec)
	# Solar hour angle
	Hra=np.linspace(-180.0, 180.0, num=nhra)
	sun=SunPosition()
	irow = 0
	icol = 0

	# creating a case folder for each new simul file
	designfolder = '%s/vtk'%basefolder

	# Estimating the azimuth and elevation vectors
	res=0
	f = open('%s/OELT_verification.csv'%basefolder,'w+')
	f.write('dec,hra,azi,ele,eta\n')
	for irow,dec in enumerate(Dec):
		for icol,hra in enumerate(Hra):
			# Getting the azimuth angle, elevation angle and DNI
			daytime,sunrise=sun.solarhour(dec, Model.latitude)
			zen=sun.zenith(Model.latitude, dec, hra)
			azi=sun.azimuth(Model.latitude, zen, dec, hra)
			if zen > 90.0:
				zen = 90.0
			ele=90.0-zen
			# Converting azimuth into SOLSTICE convention
			azimuth = azi
			elevation = ele
			azi = -(90 + azi)
			if (azi>=360.0 or azi<0.0):
				azi = (azi+360.0)%(360.0)

			if ele>8.0:
				casefolder = '%s/sunpos_%s'%(basefolder,res)
				if not os.path.exists(casefolder):
					os.makedirs(casefolder)

				DNI = DNI_ratio[ratio]*Model.get_I_Meinel(ele)
				print('DNI ratio: %s'%DNI_ratio[ratio])
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
					DNI,
					D0)
				Model.New_search_algorithm()

				# Optical postprocessing
				eta,q_results,eta_exc_intec=proces_raw_results(
						'%s/vtk/simul'%casefolder,
						casefolder)
				E[irow,icol] = eta
				f.write('%s,%s,%s,%s,%s\n'%(dec,hra,azi,ele,eta))
				res += 1

				# Read flux map
				read_data(
						'%s/vtk/'%casefolder,
						Model.r_height,
						Model.r_diameter,
						Model.num_bundle,
						Model.bins,
						flux_file=True
						)

	# Writting outputs to OELT file
	f.close()
	np.savetxt('%s/OELT.txt'%basefolder,E,fmt='%s', delimiter=',')

