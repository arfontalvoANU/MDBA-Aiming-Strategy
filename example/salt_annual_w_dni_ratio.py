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

	ratio = 3
	DNI_ratio=[1.24,1.00,0.76,0.52]

	# define a unique case folder for the user
	basefolder = os.path.join(currentdir,'case-test-a230salt-DNI_ratio_%s'%ratio)

	if not os.path.exists(basefolder):
		os.makedirs(basefolder)

	folder=basefolder
	source_path=parentdir

	# Inputs (Gemasolar)
	num_bundle=18         # Number of panels
	r_height=10.5         # Receiver height [m]
	r_diameter=8.5        # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=114.75        # Tower height [m] (Floor to receiver bottom)
	azimuth=0.0           # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=950.0             # Beam irradiance [W/m2]
	num_fp=2              # Number of flow path
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
	f = open('%s/OELT_verification.csv'%basefolder,'a')
	f.write('res,dec,hra,azi,ele,eta,sucess\n')
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
				DNI = DNI_ratio[ratio]*Model.get_I_Meinel(ele)
				if res == 0:
					print(yellow('sunpos: %s\t DNI: %s'%(res,DNI)))
					casefolder = '%s/sunpos_%s'%(basefolder,res)
					if not os.path.exists(casefolder):
						os.makedirs(casefolder)
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
					C_start = 0.7
					E_start = 2.0
					A_start = 0.5
					# Running aiming for design point
					Model.sweeping_algorithm(C_start,E_start,A_start)

					# Optical postprocessing
					eta,q_results,eta_exc_intec=proces_raw_results(
							'%s/vtk/simul'%casefolder,
							casefolder)
					E[irow,icol] = eta
					f.write('%s,%s,%s,%s,%s,%s,%s\n'%(res,dec,hra,azi,ele,eta,str(Model.sucess)))
					print(yellow('Sunpos %s -- Sucess %s'%(res,str(Model.sucess))))

					# Read flux map
					read_data(
							'%s/vtk/'%casefolder,
							Model.r_height,
							Model.r_diameter,
							Model.num_bundle,
							Model.bins,
							flux_file=True
							)
				#endif
				res += 1

	# Writting outputs to OELT file
	f.close()
	np.savetxt('%s/OELT.txt'%basefolder,E,fmt='%s', delimiter=',')

