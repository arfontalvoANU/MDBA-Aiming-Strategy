import os
import sys
import inspect
from solsticepy.cal_sun import *

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import onekey_aiming as aiming
from simulation import optical

if __name__=='__main__':
	folder=currentdir
	source_path=parentdir
	num_bundle=16		#Number of panels
	r_height=24.		#
	r_diameter=16.		#
	bins=50				#Vertical bins
	tower_h=175.		#tower height
	phi=0.0				#solar azimuth angle
	elevation=55.15	#solar elevation angle
	DNI=980.0			#Beam irradiance
	num_fp=num_bundle/2	#Two panels per flow path
	D0=60.33				#Panel tube OD
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
	#Model.simple_HT_model(10,0.7)

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
#			print(phi,elevation)
			azi.append(phi)
			ele.append(elevation)
	azi = np.array(azi)
	ele = np.array(ele)

	sim = optical(
					Dec,
					Hra,
					azi,
					ele,
					rec_slices=num_bundle,
					rec_stacks=bins,
					rec_diameter=r_diameter,
					rec_height=r_height,
					DNI=DNI,
					sunshape='pillbox',
					half_angle_deg = np.degrees(4.65e-3),
					csr = 0.000002,			#Same as Shuang
					num_rays=5e6,				#Same as Shuang
					layoutfile = Model.csv_aiming,
					rec_pos_x=0.,
					rec_pos_y=0.,
					rec_pos_z=tower_h+0.5*r_height,
					rec_abs=0.98,
					receiver='cylinder',
					hemisphere='North',
					hst_w=12.2,             #Same as Shuang
					hst_h=12.2,             #Same as Shuang
					rho_refl=0.90,          #Same as Shuang
					slope_error=1.5e-3,     #Same as Shuang
					tower_h=tower_h,
					tower_r=0.01)
	#sim.sim_group(userdefinedfolder=True,folder='OELT')
	sim.sim_single(
					azimuth=Model.phi,                 #solar azimuth angle
					elevation=Model.elevation,         #solar elevation angle
					userdefinedfolder=False,
					folder='OELT')
