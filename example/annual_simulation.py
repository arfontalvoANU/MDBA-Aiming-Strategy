import os
import sys
import inspect
from solsticepy.cal_sun import *

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import onekey_aiming as aiming
from simulation import optical
from run_solstice import *
from python_postprocessing import *
from flux_reader import *

if __name__=='__main__':
	# define a unique case folder for the user
	snum = 0
	suffix = ""
	while 1:
		import datetime
		dt = datetime.datetime.now()
		ds = dt.strftime("%a-%H-%M")
		basefolder = os.path.join(os.getcwd(),'case-%s%s'%(ds,suffix))
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
	num_bundle=16         # Number of panels
	r_height=24.0         # Receiver height [m]
	r_diameter=16.0       # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=175.0         # Tower height [m]
	phi=0.0               # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=980.0             # Beam irradiance [W/m2]
	num_fp=8              # Two panels per flow path (num_bundle/2)
	D0=60.33              # Panel tube OD [mm]
	num_rays_1=int(5e6)
	mainfolder_1=basefolder
	csv_1='%s/pos_and_aiming_new.csv'%basefolder

	Model=aiming.one_key_start(
		folder,
		source_path,
		num_bundle,
		num_fp,
		r_height,
		r_diameter,
		bins,
		tower_h,
		phi,
		elevation,
		DNI,
		D0)
	Model.New_search_algorithm()

	ndec=5
	nhra=13
	# Declination angle
	Dec=np.linspace(-23.45, 23.45, num=ndec)
	# Solar hour angle
	Hra=np.linspace(-180., 0.0, num=nhra)
	azi = []
	ele = []
	sun=SunPosition()
	for dec in Dec:
		for hra in Hra:
			# Getting the azimuth angle, elevation angle and DNI
			daytime,sunrise=sun.solarhour(dec, Model.latitude)
			theta=sun.zenith(Model.latitude, dec, hra)
			phi=sun.azimuth(Model.latitude, theta, dec, hra)
			if theta > 90.0:
				theta = 90.0
			elevation=90.0-theta
			phi = -(90 + phi)
			if (phi>=360.0 or phi<0.0):
				phi = (phi+360.0)%(360.0)
			azi.append(phi)
			ele.append(elevation)
	azi = np.array(azi)
	ele = np.array(ele)

	azimuth=270.0 - Model.phi
	if (azimuth>=360.0 or azimuth<0.0):
		azimuth = (azimuth+360.0)%(360.0)

	att_factor=Model.attenuation(Model.csv_trimmed)

	scene=SolsticeScene(
		mainfolder=basefolder,
		num_rays=num_rays_1,
		dni=DNI,
		azimuth=azimuth,
		zenith=elevation,
		att_factor=att_factor,
		csv=csv_1,
		tower_h=tower_h,
		r_cyl=num_fp,
		h_cyl=r_height,
		num_bundle=num_bundle
	)

	casefolder = basefolder + '/pos'
	if not os.path.exists(casefolder):
		os.makedirs(casefolder)

	scene.gen_YAML()
	scene.runSOLSTICE(
		savefile=casefolder,
		view=True)

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

