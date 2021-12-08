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
			num_rays=5e6, lat=34.85, ndec=5, nhra=25):
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
		self.csv_trimmed='%s/pos_and_aiming.csv'%self.folder
		self.latitude=lat
		self.dis_delta=ndec
		self.dis_omega=nhra
		self.num_rays=int(num_rays)

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
		flux_folder = '201015_N06230_thermoElasticPeakFlux_velocity'
		rec = Cyl_receiver(
				radius=0.5*self.r_diameter, 
				height=self.r_height,
				n_banks=self.num_bundle,
				n_elems=50,
				D_tubes_o=self.D0/1000.,
				D_tubes_i=self.D0/1000.-2.*1.2e-3, 
				abs_t=0.94, 
				ems_t=0.88, 
				k_coating=1.2, 
				D_coating_o=self.D0/1000.+45e-6)
		Strt = rec.flow_path(
				option='cmvNib%s'%self.num_fp,
				fluxmap_file=self.folder+'/flux-table.csv')
		rec.balance(
				HC=Na(),
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
#		results=[
#			N.sum(fluxmap*areas[ahr_map])/1e6,
#			eff_abs,
#			eff_ems,
#			N.average(T_ext),
#			N.sqrt(N.sqrt(N.average(T_ext**4))),
#			h_conv_ext,
#			N.sum(q_ref)/1e6,
#			N.sum(q_rad)/1e6,
#			N.sum(q_conv)/1e6,
#			eff_rec
#			]
		print('results:\n%s'%results)
#		aiming_results=[
#			Success,
#			Positive,
#			A_over,
#			C_safe,
#			C_net,
#			S_ratio
#			]
		for i,res in enumerate(aiming_results):
			print('aiming_results[%d]:\n%s'%(i,res))
		print('Strt:\n%s'%Strt)
#		Strt=[
#			8, 0, 7, 15, 9, 1, 6, 14, 10, 2, 5, 13, 11, 3, 4, 12
#			]

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
				flux_file=True)

		# Thermal simulation
		results,aiming_results,Strt=self.HT_model(20.,0.)

		# Print aiming_results
		print aiming_results[1]
		return aiming_results,eff_interception,Strt

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

if __name__=='__main__':
	# define a unique case folder for the user
	snum = 0
	suffix = ""
	while 1:
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
	os.system("cp example/pos_and_aiming.csv %s"%(basefolder))
	source_path=os.getcwd()

	# Inputs
	num_bundle=16         # Number of panels
	r_height=24.0         # Receiver height [m]
	r_diameter=16.0       # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=175.0         # Tower height [m]
	azimuth=0.0           # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=980.0             # Beam irradiance [W/m2]
	num_fp=8              # Number of flow path
	D0=60.33              # Panel tube OD [mm]
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

	# Running aiming for design point
	Model.New_search_algorithm()

