import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import sys
import re
import shutil
import datetime
from Deviation_aiming_new3 import *
from scipy.optimize import curve_fit
from sys import path
from python_postprocessing import *
from Open_CSPERB import *
from Open_CSPERB_plots import *
from HC import *
from Tube_materials import *
from Flux_reader import *
from cal_sun import *
from scipy import interpolate
from run_solstice import *

class one_key_start:
	def __init__(self, folder, source, num_bundle, num_fp, r_height,
			r_diameter, bins, tower_h, phi, elevation, DNI, D0,
			num_rays=5e6, lat=34.85, ndec=5, nhra=25,abs_t=0.94,ems_t=0.88):
		"""
		Instantiation of a receiver object to adjust the aiming points
		using the MDBA method from Wang et al. (2021)
		https://doi.org/10.1016/j.solener.2021.07.059
		
		Inputs:
			folder:              The case folder
			source:              Source code folder
			num_bundle:          Number of panels [slices]
			num_fp:              Number of flow paths in the receiver
			r_height:            Receiver height [m]
			r_diameter:          Receiver diameter [m]
			bins:                Number of vertical elements [stacks]
			tower_h:             Tower height (Base to receiver bottom) [m]
			phi:                 Solar azimuth angle using Duffie and Beckman's covention [degrees]
			elevation:           Solar elevation [degrees]
			DNI:                 Beam irradiance [W/m2]
			D0:                  Pipe outer diameter [mm]
			num_rays:            Number of rays
			lat:                 Latitude [degrees]
			ndec:                Number of declination angle discretisations (annual simulation only)
			nhra:                Number of hour angle discretisations (annual simulation only)
		Outputs:
			self:                Object
		"""

		# Assigning the values
		self.folder=folder
		self.source=source
		self.r_diameter=r_diameter
		self.r_height=r_height
		self.num_bundle=num_bundle
		self.num_fp=num_fp
		self.bins=bins
		self.tower_h=tower_h
		self.phi=phi
		self.elevation=elevation
		self.DNI=DNI
		self.D0=D0
		self.csv_aiming='%s/pos_and_aiming_new.csv' % self.folder
		self.csv_trimmed='%s/pos_and_aiming_salt.csv'%self.source
		self.latitude=lat
		self.dis_delta=ndec
		self.dis_omega=nhra
		self.num_rays=int(num_rays)
		self.abs_t=abs_t
		self.ems_t=ems_t

	def run_SOLSTICE(self, dni, phi, elevation, att_factor, num_rays,csv, user_folder=False, ufolder='vtk'):
		"""
		New function to create a YAML file and run SOLSTICE.
		Modified from Wang's original
		Inputs:
			dni [W/m2]:           Beam irradiance
			phi [degrees]:        Solar azimuth (SOLSTICE convention)
			elevation [degrees]:  Solar elevation
			att_factor [-]:       Attenuation factor
			num_rays [int]:       Number of rays
			csv [-]:              Field coordinates (X,Y)
		Outputs:
			demo.yaml
			demo_rcv.yaml
			simul
		"""

		# Converting azimuth into the SOLSTICE convention
		phi=270.0-phi
		if (phi>=360.0 or phi<0.0):
			phi = (phi+360.0)%(360.0)
		print('azi: %s [deg]\tele: %s [deg]\tdni: %s [W/m2]'%(phi, elevation, dni))

		# Defining folder to save yaml and simul files
		if ufolder:
			vtk_path='%s/%s'%(self.folder,ufolder)
		else:
			vtk_path='%s/vtk'%self.folder

		if os.path.exists(vtk_path):
			shutil.rmtree(vtk_path)

		# Creating SOLSTICE scene
		scene=SolsticeScene(
			mainfolder=self.folder,         #
			num_rays=num_rays,              #
			dni=dni,                        #
			azimuth=phi,                    #
			zenith=elevation,               #
			att_factor=att_factor,          #
			csv=csv,                        #
			tower_h=self.tower_h,
			r_cyl=self.r_diameter/2.,
			h_cyl=self.r_height,
			num_bundle=self.num_bundle
		)

		if not os.path.exists(vtk_path):
			os.makedirs(vtk_path)

		# Generate YAML file
		scene.gen_YAML()

		# Run SOLSTICE
		scene.runSOLSTICE(
			savefile=vtk_path,
			view=True
			)

	def attenuation(self, csv):
		"""
		Funtion to calculate the attenuation factor of a given field.
		Inputs:
			csv:                  Field coordinates (X,Y)
		Outputs:
			att_factor:           The attenuation factor [-]
		"""
		hst_info=np.loadtxt(csv,delimiter=',', skiprows=2)
		foc=hst_info[:,3]

		# Get the attenuation factor
		def func(x, b):
			return np.exp(-b * x)
		def fun_two(x):
			return 0.99321-0.0001176*x+1.97e-8*x**2
		xdata = np.linspace(0, np.max(foc), np.max(foc)*100)
		y = fun_two(xdata)
		ydata = y
		popt, pcov = curve_fit(func, xdata, ydata)
		y2 = [func(i, popt[0]) for i in xdata]
		att_factor =popt[0]
		self.att_factor = att_factor
		return att_factor

	def HT_model(self, T_amb, V_wind):
		"""
		Function to run the receiver thermal simulation
		at design point (defined by self.phi and self.elevation)
		Inputs:
			T_amb:                Ambient temperature [deg C]
			V_wind:               Wind velocity [m/s]
		Output:
			results:              
			aiming_results:       
			Strt:                 
		"""
		flux_folder = '201015_N06230_thermoElasticPeakFlux_velocity_salt'
		rec = Cyl_receiver(
				radius=0.5*self.r_diameter, 
				height=self.r_height,
				n_banks=self.num_bundle,
				n_elems=self.bins,
				D_tubes_o=self.D0/1000.,
				D_tubes_i=self.D0/1000.-2.*1.2e-3, 
				abs_t=self.abs_t, 
				ems_t=self.ems_t, 
				k_coating=1.2, 
				D_coating_o=self.D0/1000.+45e-6)
		Strt = rec.flow_path(
				option='NES-NWS',
				fluxmap_file=self.folder+'/flux-table.csv')
		rec.balance(
				HC=Solar_salt(),
				material=Haynes230(),
				T_in=290+273.15,
				T_out=565+273.15,
				T_amb=T_amb+273.15,
				h_conv_ext='SK',
				filesave=self.folder+'/flux-table',
				air_velocity=V_wind)
		flux_limits_file = \
				'%s/%s/N06230_OD%s_WT1.20_peakFlux.csv'%(
				self.source,
				flux_folder,
				round(self.D0,2))
		results,aiming_results,vel_max = tower_receiver_plots(
				files=self.folder+'/flux-table',
				efficiency=False,
				maps_3D=False,
				flux_map=False,
				flow_paths=True,
				saveloc=None,
				billboard=False,
				flux_limits_file=flux_limits_file,
				C_aiming=self.C_aiming)

		return results,aiming_results,Strt

	def aiming_loop(self,C_aiming,Exp,A_f):
		"""
		The aiming strategy loop: optical + thermal
		"""

		# The input for optical modelling
		self.C_aiming=C_aiming
		print C_aiming
		print Exp
		print A_f
		att_factor=self.attenuation(self.csv_trimmed)

		# Change of aiming points 
		aiming(
				self.folder,
				self.r_height,
				self.r_diameter,
				C_aiming,
				self.csv_trimmed,
				self.tower_h,
				self.num_bundle,
				Exp,
				A_f)

		# Optical simulation
		self.run_SOLSTICE(
				dni=self.DNI,
				phi=self.phi,
				elevation=self.elevation,
				att_factor=att_factor,
				num_rays=self.num_rays,
				csv=self.csv_aiming)

		# Optical postprocessing
		eta,q_results,eta_exc_intec=proces_raw_results(
				'%s/vtk/simul'% self.folder,
				'%s/vtk'% self.folder)
		eff_interception=eta/eta_exc_intec
		print 'Interception efficiency: ' + str(eff_interception)

		# Read flux map
		read_data(
				self.folder,
				self.r_height,
				self.r_diameter,
				self.num_bundle,
				self.bins,
				flux_file=True,
				flux_map=True)

		# Thermal simulation
		results,aiming_results,Strt=self.HT_model(20.,0.)

		# Print aiming_results
		print aiming_results[1]
		return aiming_results,eff_interception,Strt

	def get_I_Meinel(self,elevation):
		"""
		Meinel clear-sky model
		"""
		I0=1363.
		zenith=90.-elevation
		AM=1./np.cos(zenith/180.*np.pi)
		I=I0*0.7**(AM**0.678)
		return I

	def New_search_algorithm(self):
		"""
		The net one-key algorithm to get the optimised aiming points
		"""
		# Aiming extent
		C_aiming=np.zeros(self.num_bundle)
		# Equatorial aiming
		C_aiming[:]=0.0
		# Shape exponent
		Exp=np.zeros(self.num_bundle)
		# Initialising exponent
		Exp[:]=2.0
		# Asymmetry factor
		A_f=np.zeros(self.num_bundle)
		if self.num_bundle/self.num_fp == 1:
			A_f[:]=0.75
		elif self.num_bundle/self.num_fp == 2:
			A_f[:int(0.25*self.num_bundle)]=A_f[int(0.75*self.num_bundle):]=0.67
			A_f[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=0.33

		# New search algorithm
		C_aiming[:]=0.5
		aiming_results,eff_interception,Strt=\
				self.aiming_loop(C_aiming,Exp,A_f)
		gap=0.05
		while np.all(aiming_results[1])==False and np.all(C_aiming<1.):
			# Searching for E
			C_aiming_old=np.ones(self.num_bundle)
			C_aiming_old[:]=C_aiming[:]
			for i in range(self.num_bundle):
				if aiming_results[1][i]==False:
					C_aiming[Strt[i]]+=gap
					if Strt[i]==self.num_bundle-1:
						C_aiming[0]+=gap
					else:
						C_aiming[Strt[i]+1]+=gap
					C_aiming[Strt[i]-1]+=gap
				# Searching for A
				if A_f[Strt[i]]>0.5:
					if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
						A_f[Strt[i]]+=0.02
					elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
						A_f[Strt[i]]-=0.02
				else:
					if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
						A_f[Strt[i]]-=0.02
					elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
						A_f[Strt[i]]+=0.02
				# Searching for S
				if aiming_results[5][i]>0.55:
					Exp[Strt[i]]-=0.2
				elif aiming_results[5][i]<0.45:
					Exp[Strt[i]]+=0.2
			C_aiming[C_aiming-C_aiming_old>gap]=\
				C_aiming_old[C_aiming-C_aiming_old>gap]+gap
			aiming_results,eff_interception,Strt=\
				self.aiming_loop(C_aiming,Exp,A_f)

	def sweeping_algorithm(self, C_start, E_start, A_start):
		"""
		An algorithm to sweep the aiming extent (C_aiming) from a start value to 2.0
		
		Parameters:
		- C_start: Starting value for the aiming extent
		- E_start: Starting value for the shape exponent
		- A_start: Starting value for the asymmetry factor
		
		For equatorial aiming set C_start to 0.0
		"""
		# Initialising values
		C_aiming=C_start*np.ones(self.num_bundle)                         # Aiming extent
		Exp=E_start*np.ones(self.num_bundle)                       # Shape exponent
		A_f=A_start*np.ones(self.num_bundle)                       # Asymmetry factor

		# Trying the starting values
		isuccessful = False
		while np.all(C_aiming<2.) and not isuccessful:
			try:
				aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
				isuccessful = True
				print('Initialization successful')
			except ValueError:
				C_aiming[:] += 0.05
		gap=0.05
		# Sweeping algorithm for the aiming extent
		while np.all(aiming_results[1])==False and np.all(C_aiming<2.):
			for i in range(self.num_bundle):
				if aiming_results[1][i]==False:
					C_aiming[Strt[i]]+=gap
					if Strt[i]==self.num_bundle-1:
						C_aiming[0]+=gap
					else:
						C_aiming[Strt[i]+1]+=gap
					C_aiming[Strt[i]-1]+=gap
			try:
				aiming_results,eff_interception,Strt=self.aiming_loop(C_aiming,Exp,A_f)
			except ValueError:
				pass
		self.sucess = np.all(aiming_results[1])

if __name__=='__main__':
	# define a unique case folder for the user
	basefolder = os.path.join(os.getcwd(),'salt-case')
	if not os.path.exists(basefolder):
		os.makedirs(basefolder)

	folder=basefolder
	#os.system("cp example/pos_and_aiming.csv %s"%(basefolder))
	source_path=os.getcwd()

	# Inputs
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
	Model=one_key_start(
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

	# Getting the azimuth angle, elevation angle and DNI
	sun=SunPosition()
	dec=0
	hra=0
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

	DNI=Model.get_I_Meinel(ele)
	Model=one_key_start(
		folder,
		source_path,
		num_bundle,
		num_fp,
		r_height,
		r_diameter,
		bins,
		tower_h,
		azi,
		ele,
		DNI,
		D0)
	C_start = 0.5
	E_start = 2.0
	A_start = 0.5
	# Running aiming for design point
	Model.sweeping_algorithm(C_start,E_start,A_start)

