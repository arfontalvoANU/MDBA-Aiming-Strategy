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
import pickle
import argparse
import matplotlib.pyplot as plt

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

class gemasolar:
	def __init__(self,
		num_bundle=18, r_height=10.5, r_diameter=8.5, bins=50,
		tower_h=114.75, azimuth=0.0, elevation=52.38, DNI=950.0,
		num_fp=2, D0=22.40, num_rays=int(5e6), latitude=37.56,
		testcase='gemasolar',new_algorithm=False):

		self.num_bundle=num_bundle         # Number of panels
		self.r_height=r_height             # Receiver height [m]
		self.r_diameter=r_diameter         # Receiver diameter [m]
		self.bins=bins                     # Vertical bins
		self.tower_h=tower_h               # Tower height [m] (Floor to receiver bottom)
		self.azimuth=azimuth               # Solar azimuth angle [degrees]
		self.elevation=elevation           # Solar elevation angle [degrees]
		self.DNI=DNI                       # Beam irradiance [W/m2]
		self.num_fp=num_fp                 # Number of flow path
		self.D0=D0                         # Panel tube OD [mm]
		self.num_rays=num_rays             # Number of rays
		self.latitude=latitude             # Plant latitude (Seville, Spain)
		self.testcase=testcase             # Name of the case study
		self.new_algorithm = new_algorithm

	def design_point(self,ratio,C_start,E_start,A_start):

		# define a unique case folder for the user
		self.basefolder = os.path.join(currentdir,'%s_%s'%(self.testcase,ratio))

		if not os.path.exists(self.basefolder):
			os.makedirs(self.basefolder)

		self.folder=self.basefolder
		self.source_path=parentdir

		Model=aiming.one_key_start(
			self.folder, self.source_path, self.num_bundle, self.num_fp,
			self.r_height, self.r_diameter, self.bins, self.tower_h,
			self.azimuth, self.elevation, self.DNI, self.D0,
			lat=self.latitude)

		if self.new_algorithm:
			Model.New_search_algorithm()
		else:
			Model.sweeping_algorithm(C_start,E_start,A_start)

	def simulate_dni_ratio(self,ratio,C_start,E_start,A_start):

		DNI_ratio=[1.39,1.00,0.87,0.56]

		# define a unique case folder for the user
		self.basefolder = os.path.join(currentdir,'%s_%s'%(self.testcase,ratio))

		if not os.path.exists(self.basefolder):
			os.makedirs(self.basefolder)

		self.folder=self.basefolder
		self.source_path=parentdir

		# Creating receiver model
		Model=aiming.one_key_start(
			self.folder, self.source_path, self.num_bundle, self.num_fp,
			self.r_height, self.r_diameter, self.bins, self.tower_h,
			self.azimuth, self.elevation, self.DNI, self.D0,
			lat=self.latitude)

		Model.sweeping_algorithm(C_start,E_start,A_start)

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
		designfolder = '%s/vtk'%self.basefolder

		# Estimating the azimuth and elevation vectors
		res=0
		f = open('%s/OELT_verification.csv'%self.basefolder,'a')
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
					if res >= 0:
						print(yellow('sunpos: %s\t DNI: %s'%(res,DNI)))
						casefolder = '%s/sunpos_%s'%(self.basefolder,res)
						if not os.path.exists(casefolder):
							os.makedirs(casefolder)
						Model=aiming.one_key_start(
							casefolder,
							self.source_path,
							self.num_bundle,
							self.num_fp,
							self.r_height,
							self.r_diameter,
							self.bins,
							self.tower_h,
							azimuth,
							elevation,
							DNI,
							self.D0)
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
		np.savetxt('%s/OELT.txt'%self.basefolder,E,fmt='%s', delimiter=',')

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Run representative sun positions of annual simulation for a specific DNI ratio')
	parser.add_argument('--ratio', type=int, default=1, help='DNI ratio to be simulated. Default=1')
	parser.add_argument('--fpath', type=int, default=1, help='Flowpath data to be compiled. Default=1')
	parser.add_argument('--C_start', type=float, default=0.0, help='The starting value of the aiming extent')
	parser.add_argument('--E_start', type=float, default=2.0, help='The starting value of the aiming exponent')
	parser.add_argument('--A_start', type=float, default=0.5, help='The starting value of the aiming asymetry factor')
	parser.add_argument('--design', type=bool, default=False, help='Run only the design point')
	parser.add_argument('--new_algorithm', type=bool, default=False, help='Run the new search algorithm')
	args = parser.parse_args()

	cyl_receiver = gemasolar(new_algorithm=args.new_algorithm)
	if args.design:
		cyl_receiver.design_point(args.ratio, args.C_start, args.E_start, args.A_start)
	else:
		cyl_receiver.simulate_dni_ratio(args.ratio, args.C_start, args.E_start, args.A_start)
	#get_m_flow_lookup_tables(args.fpath)
