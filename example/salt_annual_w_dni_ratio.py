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

def simulate_dni_ratio(ratio,C_start,E_start,A_start):

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

def get_flux_lookup_tables(fpath):

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

	N = int(num_bundle/num_fp*bins)
	basefolder = os.path.join(currentdir,'case-a230salt')
	source_path=parentdir

	# Creating receiver model
	Model=aiming.one_key_start(
		basefolder,
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

	# Annual discretisation
	ndec=5
	nhra=25
	# Declination angle
	Dec=np.linspace(-23.45, 23.45, num=ndec)
	# Solar hour angle
	Hra=np.linspace(-180.0, 180.0, num=nhra)
	sun=SunPosition()

	for r in range(4):
		f = open('flux_a230_salt_FP%s_DNIr%s.motab'%(fpath,r),'w+')
		f.write('#1\n')
		E = np.zeros((ndec, nhra))
		F = []
		for k in range(N):
			f.write('double flux_%s(%d,%d)\n'%(k+1, ndec+1, nhra+1))
			f.write('0.0,')
			res = 0
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
					azi = -(90 + azi)
					if (azi>=360.0 or azi<0.0):
						azi = (azi+360.0)%(360.0)
					if ele>8.0:
						filename = 'case-a230salt/DNI_ratio_%s/sunpos_%s/flux-table'%(r, res)
						fileo = open(filename,'r')
						data = pickle.load(fileo)
						fileo.close()
						areas = data['areas']
						q_net = data['q_net']
						fp = data['fp']
						flux = q_net[fp[fpath-1]]/areas[fp[fpath-1]]/1e3
						E[irow,icol] = flux[k]
						res += 1
			x = np.linspace(-180.0, 180.0, nhra)
			y = np.linspace(-23.45, 23.45, ndec)
			for j in range(E.shape[1]):
				f.write('%s'%(x[j]))
				if j == E.shape[1]-1:
					f.write('\n')
				else:
					f.write(',')
			for i in range(E.shape[0]):
				f.write('%s,'%y[i])
				for j in range(E.shape[1]):
						f.write('%s'%(E[i,j]))
						if j == E.shape[1]-1:
							f.write('\n')
						else:
							f.write(',')
			f.write('\n')
		#
		f.close()

def get_m_flow_lookup_tables(fpath):

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

	N = int(num_bundle/num_fp*bins)
	basefolder = os.path.join(currentdir,'case-a230salt')
	source_path=parentdir

	# Creating receiver model
	Model=aiming.one_key_start(
		basefolder,
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

	# Annual discretisation
	ndec=5
	nhra=25
	# Declination angle
	Dec=np.linspace(-23.45, 23.45, num=ndec)
	# Solar hour angle
	Hra=np.linspace(-180.0, 180.0, num=nhra)
	sun=SunPosition()

	f = open('mflow_a230_salt_FP%s.motab'%(fpath),'w+')
	f.write('#1\n')
	for r in range(4):
		E = np.zeros((ndec, nhra))
		F = []
		f.write('double mflow_%s(%d,%d)\n'%(r, ndec+1, nhra+1))
		f.write('0.0,')
		res = 0
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
				azi = -(90 + azi)
				if (azi>=360.0 or azi<0.0):
					azi = (azi+360.0)%(360.0)
				if ele>8.0:
					filename = 'case-a230salt/DNI_ratio_%s/sunpos_%s/flux-table'%(r, res)
					print 'case-a230salt/DNI_ratio_%s/sunpos_%s/flux-table'%(r, res)
					fileo = open(filename,'r')
					data = pickle.load(fileo)
					fileo.close()
					mflow = data['m']
					n_tubes = data['n_tubes']
					E[irow,icol] = mflow[int(fpath)-1]/n_tubes[0]
					res += 1
		x = np.linspace(-180.0, 180.0, nhra)
		y = np.linspace(-23.45, 23.45, ndec)
		for j in range(E.shape[1]):
			f.write('%s'%(x[j]))
			if j == E.shape[1]-1:
				f.write('\n')
			else:
				f.write(',')
		for i in range(E.shape[0]):
			f.write('%s,'%y[i])
			for j in range(E.shape[1]):
					f.write('%s'%(E[i,j]))
					if j == E.shape[1]-1:
						f.write('\n')
					else:
						f.write(',')
		f.write('\n')
	
	f.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Run representative sun positions of annual simulation for a specific DNI ratio')
	parser.add_argument('--ratio', type=int, default=1, help='DNI ratio to be simulated. Default=1')
	parser.add_argument('--fpath', type=int, default=1, help='Flowpath data to be compiled. Default=1')
	parser.add_argument('--C_start', type=float, default=0.5, help='The starting value of the aiming extent')
	parser.add_argument('--E_start', type=float, default=2.0, help='The starting value of the aiming exponent')
	parser.add_argument('--A_start', type=float, default=0.5, help='The starting value of the aiming asymetry factor')
	args = parser.parse_args()

	#simulate_dni_ratio(args.ratio, args.C_start, args.E_start, args.A_start)
	get_m_flow_lookup_tables(args.fpath)
