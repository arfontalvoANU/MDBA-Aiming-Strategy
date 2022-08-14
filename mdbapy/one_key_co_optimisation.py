import os
import re
import shutil

from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

import random
import numpy as np
from sys import path
from scipy.optimize import curve_fit
from scipy import interpolate

from .cal_sun import SunPosition
from .cal_layout_r import radial_stagger, aiming_cylinder
from .Deviation_aiming_new3 import aiming
from .Open_CSPERB import eval_v_max, Cyl_receiver
from .Open_CSPERB_plots import tower_receiver_plots
from .HC import Na, Solar_salt
from .Tube_materials import Inconel740H, Haynes230, Incoloy800H
from .Flux_reader import read_data
from .Loss_analysis import receiver_correlation
from .output_motab import output_motab, output_matadata_motab
from .python_postprocessing import proces_raw_results, get_heliostat_to_receiver_data
from .SOLSTICE import SolsticeScene


class one_key_start:
	def __init__(self, casedir, tower_h, Q_rec, T_in, T_out, HTF, rec_material, r_diameter, r_height, fluxlimitpath, SM, oversizing, delta_r2, delta_r3, hst_w, hst_h, mirror_reflectivity, slope_error, sunshape='buie', sunshape_param=0.02, num_rays=1000000, latitude=34.85):
		"""
		casedir (str): case directory
		tower_h (float): tower height (m)
		Q_rec (float): receiver nominal output, e.g. P_pb/eta_pb*SM (W)
		T_in (float): receiver inlet temperature (K)
		T_out (float): receiver outlet temperature (K)
		HTF (str): heat transfer fluid, 'sodium' or 'salt'
		rec_material (str): receiver tube material, 'Haynes230' or 'Inconel740H'
		r_diameter (float): receiver diameter (m)
		r_height (float): receier height (m)
		fluxlimitpath (str): the directory of the flux limit file
		SM (float): solar multiple
		oversizing (float): the heliostat field oversizing factor
		delta_r2 (float): the radial distance in the second zone 
		delta_r3 (float): the radial distance in the third zone
		hst_w (float): heliostat width (m)
		hst_h (float): heliostat height (m)
		mirror_reflectivity (float): mirror reflectivity
		slope_error (float): slope error (rad)
		sunshape (str): 'pillbox' or 'buie'
		sunshape_param (float): the csr value if sunshape is Buie, the angular width (in degree) if sunshape is pillbox
		num_rays (int): number of rays in the ray tracing
		latitude (float): latitude of the plant location (degree)

		"""

	
		self.casedir=casedir
		self.fluxlimitpath=fluxlimitpath
		if not os.path.exists(casedir):
			os.makedirs(casedir)

		#shutil.copy('%s/SOLSTICE.py' % casedir, casedir) # for ray-tracing
		
		# power calculation
		#Q_rec = 111.e6/0.51*SM # receiver nominal output
		Q_rec_inc = Q_rec/0.88*oversizing # field nominal output
		num_hst=int(Q_rec_inc/0.6/980./(hst_w*hst_h)*1.3) # estimated number of hst for the large field
		print('Receiver nominal output', Q_rec)
		print('Field nominal output', Q_rec_inc)
		print('Estimated num hsts', num_hst)

		
		self.latitude=latitude # latitude of the field location
		
		# for the heliostat field
		self.hst_w=hst_w # heliostat width
		self.hst_h=hst_h # heliostat height
		self.delta_r2=delta_r2 # expanding for zone2
		self.delta_r3=delta_r3 # expanding for zone3
		dsep=0. # separation distance
		self.num_hst=num_hst # number of heliostat for the large field
		self.DM=np.sqrt(self.hst_w**2+self.hst_h**2)+dsep 
		self.csv='%s/pos_and_aiming.csv'%(self.casedir) # the large field
		self.csv_trimmed='%s/pos_and_aiming_trimmed.csv'%(self.casedir) # the trimmed field
		self.csv_aiming='%s/pos_and_aiming_new.csv'%(self.casedir) # the aiming field
		self.Q_rec_inc=Q_rec_inc
		self.oversizing=oversizing
		self.mirror_reflectivity=mirror_reflectivity
		self.slope_error=slope_error
		# for the receiver
		self.tower_h=tower_h # tower height
		self.r_diameter=r_diameter # receiver diameter
		self.r_height=r_height # receiver height
		self.bins=50 # vertical binning of receiver surface
		self.num_rays=num_rays # number of rays
		self.num_bundle=16 # A primary number of tube banks
		self.Q_rec=Q_rec
		self.T_in=T_in	
		self.T_out=T_out

		self.HTF=HTF
		if HTF=='sodium':
			self.HC=Na()
		elif HTF=='salt':
			self.HC=Solar_salt()
		else:
			print('ERROR: receiver HTF not found')			

		self.rec_material=rec_material	
		if rec_material=='Haynes230':
			self.rec_material_model=Haynes230()
		elif rec_material=='Inconel740H':
			self.rec_material_model=Inconel740H()
		elif rec_material=='Incoloy800H':
			self.rec_material_model=Incoloy800H()
		else:
			print('ERROR: receiver material not found')

		self.sunshape=sunshape
		self.sunshape_param=sunshape_param
		
		self.SM=SM # solar multiple
		
	def big_field_generation(self): # generate a large field 
		# generate and store 'pos_and_aiming.csv'
		pos_and_aim=radial_stagger(latitude=self.latitude, num_hst=self.num_hst, width=self.hst_w, height=self.hst_h, hst_z=0.7*self.hst_h, 
		  towerheight=self.tower_h+0.5*self.r_height, R1=150., delta_r_g=[0.866,self.delta_r2,self.delta_r3], 
		  dsep=0., field='surround', savedir=self.casedir, plot=False) # towerheight: the optical tower height
		
		# equatorial aiming 
		aiming_cylinder(self.r_height,self.r_diameter, pos_and_aim, self.casedir, c_aiming=0.) 
		return self.casedir
		
	def attenuation(self,csv): # calculate attenuation coefficient for SOLSTICE
		hst_info=np.loadtxt(csv,delimiter=',', skiprows=2)
		foc=hst_info[:,3]
		# to get the attenuation factor
		def func(x, b):
			return np.exp(-b * x)
		def fun_two(x):
			return 0.99321-0.0001176*x+1.97e-8*x**2
		xdata = np.linspace(0, float(np.max(foc)), int(np.max(foc)*100))
		y = fun_two(xdata)
		ydata = y
		popt, pcov = curve_fit(func, xdata, ydata)
		y2 = [func(i, popt[0]) for i in xdata]
		att_factor =popt[0]
		return att_factor
		
	def equinox(self,csv_equinox): # optical simulation at design point
		# run ray tracing simulation
		att_factor=self.attenuation(csv_equinox)		
		self.run_SOLSTICE(dni=980.,phi=0.,elevation=55.08,att_factor=att_factor,num_rays=self.num_rays,csv=csv_equinox)
		eta,q_results,eta_exc_intec=proces_raw_results('%s/simul'% self.casedir,'%s/'% self.casedir)
		print(eta,eta/eta_exc_intec)
		# calculate the efficiency of each heliostat
		hst_info=np.loadtxt(csv_equinox,delimiter=',', skiprows=2)
		num_hst=int(len(hst_info))
		Hst_data=np.arange(num_hst*10,dtype=float).reshape(num_hst,10)
		Hst_data[:,1:8]=hst_info
		q_hst_in,q_r_in=get_heliostat_to_receiver_data(simul='%s/simul'%self.casedir, DNI=980., receiver_name='cylinder')
		Hst_data[:,0]=np.arange(num_hst)
		Hst_data[:,8]=q_r_in
		Hst_data[:,9]=q_r_in[:]/q_hst_in[:] # heliostat optical efficiency
		return Hst_data
		
	def annual_big_field(self): # optical simulations of the large field at different sun positions
		att_factor=self.attenuation(self.csv)
		hst_info=np.loadtxt(self.csv,delimiter=',', skiprows=2)
		self.num_hst=int(len(hst_info))
		N=10  # for lammda, ecliptic angle
		M=24  # for omega, hour angle
		F=np.zeros((N,M,self.num_hst)) # the optical efficiency lookup table for each hst
		Lammda=np.linspace(-np.pi,np.pi,N+1)
		Omega=np.linspace(-np.pi,np.pi,M+1)
		sun=SunPosition()
		# to find the index of symmetric hst
		Sym_index=np.ones(len(hst_info))
		for i in range(len(hst_info)):
			Coor=np.array([-hst_info[i,0],hst_info[i,1]])
			Diff=Coor-hst_info[:,:2]
			Sym_index[i]=np.where(np.sqrt(Diff[:,0]**2+Diff[:,1]**2)==min(np.sqrt(Diff[:,0]**2+Diff[:,1]**2)))[0][0]
		
		for n in range(3,8):
			for m in range(int(0.5*M)+1):  # symmetry in East-West direction
				delta = 23.4556*np.sin(Lammda[n])
				theta=sun.zenith(self.latitude, delta, Omega[m]/np.pi*180.)
				phi=sun.azimuth(self.latitude, theta, delta, Omega[m]/np.pi*180.)
				elevation=90.-theta
				if elevation<=15.:
					continue
				dni=self.get_I(elevation)
				self.run_SOLSTICE(dni=dni,phi=phi,elevation=elevation,att_factor=att_factor,num_rays=self.num_rays,csv=self.csv)
				q_hst_in,q_r_in=get_heliostat_to_receiver_data(simul='%s/simul'%self.casedir, DNI=dni, receiver_name='cylinder')
				# calculate the efficiency of each heliostat
				F[n,m,:]=q_r_in[:]/q_hst_in[:]
			
			for m in range(int(0.5*M)+1,M): # in the afternoon
				for i in range(len(hst_info)):
					F[n,m,i]=F[n,M-m,int(Sym_index[i])]
					
		# symmetric due to season		
		for n in range(3):
			F[n,:,:]=F[5-n,:,:]
		for n in range(8,10):
			F[n,:,:]=F[15-n,:,:]
		
		In_energy=np.zeros(self.num_hst) # the input energy to the single hst, Wh 
		Out_field_energy=np.zeros(self.num_hst) # the annual energy at receiver,Wh
		A_hst=self.hst_w*self.hst_h # m2 for single hst
		opt_eff=0.
		for n in range(N):
			for m in range(M):
				delta = 23.4556*np.sin(Lammda[n])
				theta=sun.zenith(self.latitude, delta, Omega[m]/np.pi*180.)
				phi=sun.azimuth(self.latitude, theta, delta, Omega[m]/np.pi*180.)
				elevation=90.-theta
				if elevation<=0:
					continue
				dni=self.get_I(elevation)
				for h in range(self.num_hst):
					In_energy[h]+=dni*A_hst
					Out_field_energy[h]+=F[n,m,h]*dni*A_hst
		
		# output the annual efficiency of each heliostat
		Hst_data=np.arange(self.num_hst*11,dtype=float).reshape(self.num_hst,11)
		Hst_data[:,1:8]=np.loadtxt(self.csv,delimiter=',', skiprows=2)
		Hst_data[:,0]=np.arange(self.num_hst)
		Hst_data[:,8]=Out_field_energy
		Hst_data[:,9]=Out_field_energy[:]/In_energy[:] # heliostat annual optical efficiency
		
		# to get the output of single heliostats at design point
		Hst_data_equinox=self.equinox(csv_equinox=self.csv)
		Hst_data[:,10]=Hst_data_equinox[:,8]
		np.savetxt('%s/Hst_annual_eff.csv'%self.casedir, Hst_data, fmt='%s', delimiter=',')
		
	def determine_field(self): # field trimming
		Hst_data=np.loadtxt('%s/Hst_annual_eff.csv'%self.casedir,delimiter=',') 		
		self.plot_hst_eff(self.DM,self.csv,Hst_data[:,9])
		Hst_data_ranked = Hst_data[np.argsort(Hst_data[:,9])[::-1]] # ranked by the annual efficiency
		cumul_p = np.add.accumulate(Hst_data_ranked[:,10]) # summary of abs energy at design point
		Hst_data_partial_ranked = Hst_data_ranked[~(cumul_p>self.Q_rec_inc)]
		print(np.add.accumulate(Hst_data_partial_ranked[:,10])[-1])
		title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
		pos_and_aiming_new=np.append(title, Hst_data_partial_ranked[:,1:8])
		pos_and_aiming_new=pos_and_aiming_new.reshape(int(len(pos_and_aiming_new)/7), 7)
		np.savetxt('%s/pos_and_aiming_trimmed.csv'%(self.casedir), pos_and_aiming_new, fmt='%s', delimiter=',')
		#self.plot_hst_eff(self.DM,self.csv_trimmed,Hst_data_partial_ranked[:,-2])
		return Hst_data_partial_ranked[:,0]
	
	def flow_path_sodium(self): # flow path and pipe diameter selections for sodium receivers
		Q_demand=self.Q_rec
		velocity_limit=2.44
		oversizing=self.oversizing
		hst_info=np.loadtxt(self.csv_trimmed,delimiter=',', skiprows=2) 
		D0_group=np.array([33.40,42.16,48.26,60.33,73.03]) # nominal pipe diameters
		Candidate=np.array([])
		for n_p in range(1,2): # consider single, double and trible passes
			num_bundle=16
			num_fp=int(num_bundle/n_p)
			if n_p==3:
				num_bundle=12
				num_fp=int(num_bundle/n_p)

			rec = Cyl_receiver(
				radius=0.5*self.r_diameter, 
				height=self.r_height, 
				n_banks=num_bundle, 
				n_elems=self.bins, 
				D_tubes_o=60.33/1000., 
				D_tubes_i=60.33/1000.-2.*1.2e-3, 
    			abs_t=0.98, 
				ems_t=0.91, 
				k_coating=1.2, 
				D_coating_o=60.33/1000.+45e-6)

			pattern = 'cmvNit'
			Strt=np.array([8, 0, 7, 15, 9, 1, 6, 14, 10, 2, 5, 13, 11, 3, 4, 12])
			if n_p==3:
				pattern = 'cm3Nit'
				Strt=np.array([6,2,0,5,9,11,7,3,1,4,8,10])
			Azimuth_boundary=np.zeros(num_bundle+1) # from -90 to 270
			for i in range(num_bundle+1):
				Azimuth_boundary[i]=-90.+360./num_bundle*i	

			# number of heliostat in each field sector
			N_hst_sector=np.zeros(num_bundle,dtype=int)
			for i in range(num_bundle):
				for j in range(int(len(hst_info))):
					x=hst_info[j,0]
					y=hst_info[j,1]
					if abs(x)<1e-6:
						azimuth=90.
					else:	
						azimuth=np.arctan(y/x)/np.pi*180.
					if x<0:
						azimuth+=180.
					if azimuth>=Azimuth_boundary[i] and azimuth<Azimuth_boundary[i+1]:
						N_hst_sector[i]+=1

			# number of heliostats in each flow path
			N_hst_fp=np.zeros(num_fp,dtype=int) 
			for f in range(num_fp):
				for p in range(n_p):
					N_hst_fp[f]+=N_hst_sector[Strt[f*n_p+p]] 
			e_net_fp=Q_demand*max(N_hst_fp)/len(hst_info)

			vmax, n_t=eval_v_max(
				e_net_fp=e_net_fp, 
				HC=self.HC, 
				T_in=self.T_in, 
				T_out=self.T_out, 
				W_abs=np.pi*self.r_diameter, 
				n_b=self.num_bundle, 
				D_tube_o=D0_group/1000., 
				D_tube_in=D0_group/1000.-2.*1.2e-3, 
				pipe_spacing=1e-3)

			ID_1=np.asarray(np.where((vmax<velocity_limit) & (vmax>=velocity_limit/oversizing)))[0]
			if ID_1.size != 0:
				for i in range(int(len(ID_1))):
					Candidate=np.append(Candidate,[num_bundle,num_fp,D0_group[ID_1[i]],pattern])

			ID_2=np.asarray(np.where(vmax<velocity_limit/oversizing))[0]
			if ID_2.size != 0:
				for i in range(ID_2[0],min(ID_2[0]+1,len(D0_group))):
					Candidate=np.append(Candidate,[num_bundle,num_fp,D0_group[i],pattern])

		Candidate=Candidate.reshape(int(len(Candidate)/4),4)
		print("FLow path candidate", Candidate)
		
		Results=np.array([])
		for i in range(int(len(Candidate))):
			self.num_bundle=int(Candidate[i,0])
			self.num_fp=int(Candidate[i,1])
			self.num_pass=int(self.num_bundle/self.num_fp)
			self.D0=float(Candidate[i,2])
			self.pattern=Candidate[i,3]
			Defocus,MDBA_results=self.MDBA_aiming_new(dni=980,phi=0.,elevation=55.08)
			Results=np.append(Results,MDBA_results[0]) #MDBA_results[0] is the overall optic efficiency

		
		idx=np.where(Results==max(Results))[0][0]
		print('  candidate index with the max eta_opt', idx)

		self.num_bundle=int(Candidate[idx,0])
		self.num_fp=int(Candidate[idx,1])
		self.num_pass=int(self.num_bundle/self.num_fp)
		self.D0=float(Candidate[idx,2])
		self.pattern=Candidate[idx,3]
		np.savetxt('%s/flowpath.csv'%(self.casedir), np.array([self.num_bundle,self.num_fp,self.D0,self.pattern]), fmt='%s', delimiter=',')
	
	def flow_path_salt(self,num_bundle,num_fp,D0,pattern): 
		np.savetxt('%s/flowpath.csv'%(self.casedir), np.array([num_bundle,num_fp,D0,pattern]), fmt='%s', delimiter=',')
		self.num_bundle=num_bundle
		self.num_fp=num_fp
		self.D0=D0
		self.pattern=pattern
		self.num_pass=num_bundle/2
		
	def MDBA_aiming_new(self,dni,phi,elevation): # MDBA
		self.num_hst=int(len(np.loadtxt(self.csv_trimmed,delimiter=',', skiprows=2))) # the num hst for the large field
		att_factor=self.attenuation(self.csv_trimmed)
		GA_tot=dni*self.num_hst*self.hst_w*self.hst_h
		#print GA_tot
		print('')
		print('Start MDBA')

		# initialization for aiming
		C_aiming=np.zeros(self.num_bundle)
		C_aiming[:]=0.4
		Exp=np.zeros(self.num_bundle)
		Exp[:]=3.0
		A_f=np.zeros(self.num_bundle)
		if self.pattern =='cmvNit': # flow path for sodium, specific for num_bundle=16
			if self.num_bundle/self.num_fp == 1:
				A_f[:]=0.75
			elif self.num_bundle/self.num_fp == 2:
				A_f[:int(0.25*self.num_bundle)]=A_f[int(0.75*self.num_bundle):]=0.33
				A_f[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=0.67
		elif self.pattern=='NES-NWS': # this pattern is top injection from two northern banks
			for t in range(self.num_bundle):
				if t<0.5*self.num_bundle and t%2==0:
					A_f[t]=0.67
				elif t<0.5*self.num_bundle and t%2!=0:
					A_f[t]=0.33
				elif t>=0.5*self.num_bundle and t%2==0:
					A_f[t]=0.33
				else:
					A_f[t]=0.67
		self.C_aiming=C_aiming

		Hst_info,Hst_stand=aiming(
			self.casedir,
			self.r_height,
			self.r_diameter,
			C_aiming,
			self.csv_trimmed,
			self.tower_h,
			self.num_bundle,
			Exp,
			A_f,
			stand_by=False)

		usage=float(self.num_hst-sum(Hst_stand))/self.num_hst
		self.run_SOLSTICE(dni=dni,phi=phi,elevation=elevation,att_factor=att_factor,num_rays=self.num_rays,csv=self.csv_aiming)
		eta,q_results,eta_exc_intec=proces_raw_results('%s/simul'% self.casedir,'%s/'% self.casedir)
		cos_loss=q_results[1]/q_results[0]*usage # ratio
		ref_loss=q_results[3]/q_results[0]*usage
		sb_loss=(q_results[2]+q_results[4])/q_results[0]*usage
		att_loss=q_results[5]/q_results[0]*usage
		spi_loss=q_results[6]/q_results[0]*usage
		print('	Interception efficiency: ' + str(eta/eta_exc_intec))
		read_data(self.casedir,self.r_height,self.r_diameter,self.num_bundle,self.bins,flux_file=True,flux_map=False)
		results,aiming_results,vel_max,Strt=self.HT_model(35.,0.)
		Vel_bool=vel_max<2.44
		print('	Q_abs', q_results[-1])
		print('	Aiming extent:     [' + ','.join('%.2f'%x for x in C_aiming) + ']')
		print('	Shape exponent:    [' + ','.join('%.2f'%x for x in Exp) + ']')
		print('	Asymmetry factor:  [' + ','.join('%.2f'%x for x in A_f) + ']')
		print('	aiming_results[1]: %s/%s'%(np.sum(aiming_results[1]), len(aiming_results[1])))  # for each tube bank
		print('	Vel_bool: %s/%s'%(np.sum(Vel_bool), len(Vel_bool))) # for each flow path
		print('	vel_max:', np.max(vel_max))

		gap=0.01
		Defocus=np.all(aiming_results[1])==False
		title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])

		# search algorithm
		ite1=0
		pos_and_aiming_stand_by=np.array([])
		while ((np.all(aiming_results[1])==False or np.all(Vel_bool)==False) and ite1<100): #and np.all(C_aiming<1.):
			print('		Iteration', ite1)	

			if np.all(C_aiming<1.)==False:
				# instead of extent E, defocus high-foc heliostats
				Tude_index=np.full(self.num_bundle, True, dtype=bool)
				for i in range(self.num_bundle):
					if aiming_results[1][i]==False or Vel_bool[i]==False:
						Tude_index[int(Strt[i]-1)]=Tude_index[int(Strt[i])]=False
						if Strt[i]==self.num_bundle-1:
							Tude_index[0]=False
						else:
							Tude_index[int(Strt[i]+1)]=False
				pos_and_aiming_new=np.array([])
				for i in range(self.num_bundle):
					if Tude_index[i]==False:
						foc=Hst_info[i][:,3]
						pos_and_aiming_new=np.append(pos_and_aiming_new, Hst_info[i][np.argsort(foc)[::1]][:int(0.95*len(Hst_info[i]))])
						pos_and_aiming_stand_by=np.append(pos_and_aiming_stand_by, Hst_info[i][np.argsort(foc)[::1]][int(0.95*len(Hst_info[i])):])
					else:
						pos_and_aiming_new=np.append(pos_and_aiming_new, Hst_info[i])
					#print (i,len(Hst_info[i]),len(pos_and_aiming_new)/7,len(pos_and_aiming_stand_by)/7)
				pos_and_aiming_new=np.append(title,pos_and_aiming_new)
				pos_and_aiming_new=pos_and_aiming_new.reshape(int(len(pos_and_aiming_new)/7), 7)
				np.savetxt('%s/pos_and_aiming_new.csv' % self.casedir, pos_and_aiming_new, fmt='%s', delimiter=',')	
				print('		length: %s'%str(len(pos_and_aiming_new)-2))
			else:
				C_aiming_old=np.ones(self.num_bundle)
				C_aiming_old[:]=C_aiming[:]
				for i in range(self.num_bundle):
					if aiming_results[1][i]==False or Vel_bool[i]==False:
						if C_aiming[int(Strt[i])]>0.8:
							gap=0.01
						C_aiming[int(Strt[i])]+=gap
						if np.all(C_aiming<1.0):
							if Strt[i]==self.num_bundle-1:
								C_aiming[0]+=gap
							else:
								C_aiming[int(Strt[i]+1)]+=gap
							C_aiming[int(Strt[i]-1)]+=gap
			C_aiming[C_aiming-C_aiming_old>gap]=C_aiming_old[C_aiming-C_aiming_old>gap]+gap
			if aiming_results[5].size!=0:
				for i in range(self.num_bundle):
					if aiming_results[1][i]==False:
						# for A
						if A_f[int(Strt[i])]>0.5:
							if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
								A_f[int(Strt[i])]+=0.01
							elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
								A_f[int(Strt[i])]-=0.01
						else:
							if (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])<-0.1:
								A_f[int(Strt[i])]-=0.01
							elif (aiming_results[3][i]-aiming_results[4][i])/abs(aiming_results[4][i])>0.1:
								A_f[int(Strt[i])]+=0.01
						# for S
						if aiming_results[5][i]>0.55:
							Exp[int(Strt[i])]-=0.05
						elif aiming_results[5][i]<0.45:
							Exp[int(Strt[i])]+=0.05
			C_aiming[C_aiming>1.]=1.0
			if self.pattern=='NES-NWS':
				A_f[2:self.num_bundle-2]=0.5 # no tilted aiming
			Hst_info,Hst_stand=aiming(
				self.casedir,
				self.r_height,
				self.r_diameter,
				C_aiming,
				self.csv_aiming,
				self.tower_h,
				self.num_bundle,
				Exp,
				A_f,
				stand_by=False)
			usage=float(len(np.loadtxt(self.csv_aiming,delimiter=',', skiprows=2)))/self.num_hst
			self.run_SOLSTICE(dni=dni,phi=phi,elevation=elevation,att_factor=att_factor,num_rays=self.num_rays,csv=self.csv_aiming)
			eta,q_results,eta_exc_intec=proces_raw_results('%s/simul'% self.casedir,'%s/'% self.casedir)
			print('	 	Interception efficiency (iteration %s): %s'%(ite1, eta/eta_exc_intec))
			print('	 	Q_abs (iteration %s): %s'%(ite1, q_results[-1]))
			cos_loss=q_results[1]/q_results[0]*usage # ratio
			ref_loss=q_results[3]/q_results[0]*usage
			sb_loss=(q_results[2]+q_results[4])/q_results[0]*usage
			att_loss=q_results[5]/q_results[0]*usage
			spi_loss=q_results[6]/q_results[0]*usage
			read_data(self.casedir,self.r_height,self.r_diameter,self.num_bundle,self.bins,flux_file=True,flux_map=False)
			results,aiming_results,vel_max,Strt=self.HT_model(35.,0.,overflux=not np.all(aiming_results[1]))
			Vel_bool=vel_max<2.44
			print('		Aiming extent:     [' + ','.join('%.2f'%x for x in C_aiming) + ']')
			print('		Shape exponent:    [' + ','.join('%.2f'%x for x in Exp) + ']')
			print('		Asymmetry factor:  [' + ','.join('%.2f'%x for x in A_f) + ']')
			print('		aiming_results[1]: %s/%s'%(np.sum(aiming_results[1]), len(aiming_results[1])))
			print('		Vel_bool: %s/%s'%(np.sum(Vel_bool), len(Vel_bool)))
			print('		vel_max:', np.max(vel_max))
			ite1+=1

			
		# save the stand-by hst
		pos_and_aiming_stand_by=np.append(title,pos_and_aiming_stand_by)
		pos_and_aiming_stand_by=pos_and_aiming_stand_by.reshape(int(len(pos_and_aiming_stand_by)/7), 7)
		np.savetxt('%s/pos_and_aiming_stand_by.csv' % self.casedir, pos_and_aiming_stand_by, fmt='%s', delimiter=',')
		print('		length stand by: %s'%str(len(pos_and_aiming_stand_by)-2))

		MDBA_results=[q_results[-1]/GA_tot,1-usage,cos_loss,ref_loss,sb_loss,att_loss,spi_loss]
		
		return Defocus,MDBA_results

	
	def annual_trimmed_field(self): # OELT and RELT generations

		self.num_hst=int(len(np.loadtxt(self.csv_trimmed,delimiter=',', skiprows=2))) # the num hst for the large field
		flowpath=np.loadtxt('%s/flowpath.csv'%(self.casedir), dtype=str, delimiter=',')
		self.num_bundle=flowpath[0].astype(int)
		self.num_fp=flowpath[1].astype(int)
		self.num_pass=int(self.num_bundle/self.num_fp)
		self.D0=flowpath[2].astype(float)
		self.pattern=flowpath[3]

		months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
		
		N=10  # for lammda, ecliptic longitude
		M=24  # for omega
		F=np.arange((N+1)*(M+1)*7,dtype=float).reshape(N+1,M+1,7) # the OELT, 0-field eff, 1-unavail,2-cosine, 3-reflective,4-sb,5-att,6-spi
		
		Lammda=np.linspace(-np.pi,np.pi,N+1)
		Omega=np.linspace(-np.pi,np.pi,M+1)
		att_factor=self.attenuation(self.csv_trimmed)
		T_amb_g=np.linspace(-5.6,45,5)
		V_wind_g=np.linspace(0.,18.,5)
		results_table=np.array([]) # RELT
		#DNI_ratio=[1.24,1.00,0.76,0.52]
		DNI_ratio=[1.00]
		sun=SunPosition()
		F[:,:,:]=0.
		Defocus=F[:,:,0]==0
		d=0
		for d in range(int(len(DNI_ratio))):
			for n in range(3,8):

				for m in range(int(0.5*M)+1):
					if Defocus[n,m]==False:
						continue
					print('	DNI ratio',  d, 'n', n, 'm', m)
					delta = 23.4556*np.sin(Lammda[n])
					theta=sun.zenith(self.latitude, delta, Omega[m]/np.pi*180.)
					phi=sun.azimuth(self.latitude, theta, delta, Omega[m]/np.pi*180.)
					elevation=90.-theta
					if elevation<=8.:
						continue
					dni=self.get_I(elevation)
					Defocus[n,m],F[n,m,:]=self.MDBA_aiming_new(dni=dni*DNI_ratio[d],phi=phi,elevation=elevation)
					#print('	OPTICAL EFFICIENCY',F[n,m,:])


					# generate the RELT
					if d==0:
						if n==3 or n==5 or n==7:
							max_flux=read_data(self.casedir,self.r_height,self.r_diameter,self.num_bundle,self.bins,flux_file=True)
							for i in range(int(len(T_amb_g))):
								print('T_abm_g, 3', T_amb_g[i],3.)
								results,aiming_results,vel_max,Strt=self.HT_model(T_amb_g[i],3.)
								#print(results)
								if np.isnan(results[3]):
									continue 
								results_table=np.append(results_table,[results[0],T_amb_g[i],3.,dni,max_flux])
								results_table=np.append(results_table,results[3:])

							for i in range(int(len(V_wind_g))):
								print('20 V_wind_g', 20.,V_wind_g[i])
								results,aiming_results,vel_max,Strt=self.HT_model(20.,V_wind_g[i])
								#print(results)
								if np.isnan(results[3]):
									continue 
								results_table=np.append(results_table,[results[0],20.,V_wind_g[i],dni,max_flux])
								results_table=np.append(results_table,results[3:])
					#print('results_table', results_table)					



				for m in range(int(0.5*M)+1,M+1):
					F[n,m,:]=F[n,M-m,:]
					Defocus[n,m]=Defocus[n,M-m]
			
			# symmetric due to season		
			for n in range(3):
				F[n,:,:]=F[5-n,:,:]
				Defocus[n,:]=Defocus[5-n,:]
			
			for n in range(8,11):
				F[n,:,:]=F[15-n,:,:]
				Defocus[n,:]=Defocus[15-n,:]
			
			for n in range(N+1):
				for m in range(M+1):
					if F[n,m,0]==0:
						Defocus[n,m]=False
			np.savetxt('%s/Defocus_%s.csv'%(self.casedir,DNI_ratio[d]), Defocus, fmt='%s', delimiter=',')
			# to output F
			F_output=np.arange((N+2)*(M+2),dtype=float).reshape(N+2,M+2)
			F_output[0,1:]=Omega/np.pi*180.
			F_output[1:,0]=Lammda/np.pi*180.
			F_output[1:,1:]=F[:,:,0]
			np.savetxt('%s/F_optic_%s.csv'%(self.casedir,DNI_ratio[d]), F_output, fmt='%s', delimiter=',')
			F_output[1:,1:]=F[:,:,1]
			np.savetxt('%s/F_unavail_%s.csv'%(self.casedir,DNI_ratio[d]), F_output, fmt='%s', delimiter=',')
			#d+=1
		#print('****RELT****')
		#print(results_table)
		#print('')
		# to output RELT
		title=np.array(['Qin', 'T_amb', 'Wind_speed', 'DNI','Peak_flux','T_ext_mean','T_ext_mean','h_ext','q_refl','q_emi','q_conv','R_eff'])
		results_table=np.append(title, results_table)
		results_table=results_table.reshape(int(len(results_table)/12), 12)
		np.savetxt('%s/RELT.csv'%self.casedir, results_table, fmt='%s', delimiter=',')
		
		# to output design point data
		self.MDBA_aiming_new(dni=980,phi=0.,elevation=55.08)
		eta,q_results,eta_exc_intec=proces_raw_results('%s/simul'% self.casedir,'%s/'% self.casedir)
		self.num_hst=int(len(np.loadtxt('%s/pos_and_aiming_trimmed.csv'%self.casedir,delimiter=',', skiprows=2))) # for the big field
		eta=q_results[-1]/(self.num_hst*980.*self.hst_w*self.hst_h)
		print(eta)
		read_data(self.casedir,self.r_height,self.r_diameter,self.num_bundle,self.bins,flux_file=True)
		results,aiming_results,vel_max,Strt=self.HT_model(25.,3.)
		print(results[-1])
		Equinox=[eta,results[-1]]
		np.savetxt('%s/Equinox.csv'%self.casedir, Equinox, fmt='%s', delimiter=',')

		Equinox=np.loadtxt('%s/Equinox.csv'%self.casedir,delimiter=',')
		self.csv_trimmed='%s/pos_and_aiming_trimmed.csv'%self.casedir
		self.num_hst=int(len(np.loadtxt(self.csv_trimmed,delimiter=',', skiprows=2)))
		print(self.num_hst)
		table=np.loadtxt('%s/F_optic_1.0.csv'%self.casedir,delimiter=',')
		coefs_T,coefs,eff_abs,eff_emi,A_rec=receiver_correlation(self.r_diameter,self.r_height,folder=self.casedir)

		
		output_matadata_motab(table, field_type='surrounding', aiming='MDBA', n_helios=self.num_hst, A_helio=self.hst_w*self.hst_h, eff_design=Equinox[0], 
		  d_receiver=self.r_diameter, h_receiver=self.r_height, H_tower=self.tower_h, eff_rec_design=Equinox[1], coefs_T=coefs_T, coefs=coefs, eff_abs=eff_abs, eff_emi= eff_emi,SM=self.SM, savedir='%s/OELT_Solstice.motab'%self.casedir)
				
	def HT_model(self,T_amb,V_wind,overflux=True): # receiver heat balance model

		#print('D0', self.D0)
		#print('Flow pattern',self.pattern+str(self.num_fp))
		rec = Cyl_receiver(
			radius=0.5*self.r_diameter, 
			height=self.r_height, 
			n_banks=self.num_bundle, 
			n_elems=self.bins, 
			D_tubes_o=self.D0/1000., 
			D_tubes_i=self.D0/1000.-2.*1.5e-3, 
		    abs_t=0.94, 
			ems_t=0.88, 
			k_coating=1.2, 
			D_coating_o=self.D0/1000.+45e-6)
		if self.HTF=='sodium':
			Strt=rec.flow_path(option=self.pattern+str(self.num_fp),fluxmap_file=self.casedir+'/flux-table.csv') # 16 flow paths
		elif self.HTF=='salt':
			Strt=rec.flow_path(option=self.pattern,fluxmap_file=self.casedir+'/flux-table.csv')
		rec.balance(
			HC=self.HC, 
			material=self.rec_material_model, 
			T_in=self.T_in, 
			T_out=self.T_out, 
			T_amb=T_amb+273.15, 
			h_conv_ext='SK', 
			filesave=self.casedir+'/flux-table',air_velocity=V_wind)

		if self.rec_material=='Haynes230':
			material_name='N06230'
		elif self.rec_material=='Inconel740H':
			material_name='Incoloy800H'
		elif self.rec_material=='Incoloy800H':
			material_name='N08811'
	
		flux_limits_file='%s/%s_OD%.2f_WT1.20_peakFlux.csv'%(self.fluxlimitpath,material_name, self.D0)
	
		results,aiming_results,vel_max=tower_receiver_plots(
			files=self.casedir+'/flux-table', 
			efficiency=False, 
			maps_3D=False, 
			flux_map=False, 
			flow_paths=True,
			saveloc=None, 
			billboard=False, 
			flux_limits_file=flux_limits_file,
			C_aiming=self.C_aiming,overflux=overflux)
		
		vel_max_2=np.ones(self.num_bundle) # no considered velocity limit for salt receiver
		if self.HTF=='sodium':
			for i in range(self.num_fp):
				vel_max_2[self.num_pass*i]=vel_max[i]
				if self.num_pass>1:
					vel_max_2[self.num_pass*i+1]=vel_max[i]

			# plot the velocity at different banks
			plot_vel=False
			if plot_vel==True:
				X=np.arange(self.num_bundle+1)
				X_1=X[1:]-0.5
				X_2=X[1:]+0.5
				vel_max_new=np.zeros(self.num_bundle)
				vel_max_new[int(Strt[:])]=vel_max_2[:]
				fig = plt.figure(figsize=(8,6))
				ax=plt.subplot(111)
				ax.hlines(vel_max_new, X_1, X_2, colors='black')
				ax.hlines(y=2.44, xmin=0, xmax=20, colors='red')
				for i in range(self.num_bundle):
					plt.axvline(x=X_1[i],c='black',linestyle='--',linewidth=0.3)
				plt.axvline(X_2[-1], c='black',linestyle='--',linewidth=0.3)
				plt.xlabel('Tube bank index',fontsize=16)
				plt.ylabel('${\it{V}_\mathrm{HTF}}$ (m/s)',fontsize=16)
				ax.tick_params(axis='both', which='major', labelsize=16,direction='in')
				plt.ylim(min(vel_max_new)-0.2,3.0)
				plt.xlim(0,self.num_bundle+1)
				plt.savefig(open('%s/Velocity.png'%self.casedir,'w'), dpi=400)
				plt.close('all')
			
		return results,aiming_results,vel_max_2,Strt
	
	def get_I(self,elevation): # clear-sky DNI 
		I0=1363.
		zenith=90.-elevation
		AM=1./np.cos(zenith/180.*np.pi)
		I=I0*0.7**(AM**0.678)
		return I
	
	def run_SOLSTICE(self,dni,phi,elevation,att_factor,num_rays,csv): # the input is not a solstice style
		# transfer into SOLSTICE convention
		phi=270.-phi
		if phi > 360.:
			phi=phi-360.
		print('	Sun position (%.2f, %.2f), DNI %.1f'%(phi, elevation,dni))

		scene=SolsticeScene(mainfolder=self.casedir,num_rays=num_rays,dni=dni,azimuth=phi,zenith=elevation,att_factor=att_factor,csv=csv,tower_h=self.tower_h,r_cyl=self.r_diameter/2.,h_cyl=self.r_height,num_bundle=self.num_bundle, hst_w=self.hst_w, hst_h=self.hst_h, mirror_reflectivity=self.mirror_reflectivity, slope_error=self.slope_error, sunshape=self.sunshape, sunshape_param=self.sunshape_param)
		scene.gen_YAML()
		scene.runSOLSTICE(savefile=self.casedir, view=True)
	
	def plot_hst_eff(self,DM,csv,Eff_hst,plot_stand_by=False):
		hst_info=np.loadtxt(csv,delimiter=',', skiprows=2)
		N_hst=int(len(Eff_hst))
		fig, ax = plt.subplots(figsize=(16,9))
		cmap = cm.rainbow
		norm = mpl.colors.Normalize(vmin = 0.3, vmax = 0.9)
		smap = cm.ScalarMappable(norm = norm, cmap = cmap)
		smap.set_array([])
		color_bar = fig.colorbar(mappable = smap, ax = ax, orientation = 'vertical',aspect=30,pad=0.005)
		color_bar.set_label('$\it{\eta}_{\mathrm{hst}}$',fontsize=20)
		color_bar.ax.tick_params(labelsize=20) 
		
		Azimuth=np.arange(0.,1.+1./int(self.num_bundle),1./int(self.num_bundle))*2.*np.pi
		Lenght=np.ones(self.num_bundle)*1800.
		Lenght[int(0.25*self.num_bundle):int(0.75*self.num_bundle)]=1900
		Lenght[3]=Lenght[12]=1850
		for i in range(int(0.5*self.num_bundle)):
			plt.plot([2500*np.cos(Azimuth[i]),2500*np.cos(Azimuth[i]+np.pi)],[2500*np.sin(Azimuth[i]),2500*np.sin(Azimuth[i]+np.pi)],'--',linewidth=0.5,c='grey')
		for i in range(self.num_bundle):
			plt.text(x=Lenght[i]*np.cos((Azimuth[i]+Azimuth[i+1])/2.-0.5*np.pi), y=Lenght[i]*np.sin((Azimuth[i]+Azimuth[i+1])/2.-0.5*np.pi), s='%s'%(i+1), ha='center',va='center',fontsize=20)
			cell = Circle(xy = (Lenght[i]*np.cos((Azimuth[i]+Azimuth[i+1])/2.-0.5*np.pi),Lenght[i]*np.sin((Azimuth[i]+Azimuth[i+1])/2.-0.5*np.pi)), radius=100., edgecolor = 'black', 
			 linewidth=0.5,facecolor = 'white')
			ax.add_patch(cell)
		
		for i in range(N_hst):
			#print hst_info[i,0], hst_info[i,1]
			cell = Circle(xy = (hst_info[i,0], hst_info[i,1]), radius=DM/2., edgecolor = 'black', 
			 linewidth=0.1,facecolor = smap.to_rgba(Eff_hst[i]))
			ax.add_patch(cell)
		if plot_stand_by!=False:
			hst_info_stand_by=np.loadtxt(plot_stand_by,delimiter=',', skiprows=2)
			for i in range(int(len(hst_info_stand_by))):
				cell = Circle(xy = (hst_info_stand_by[i,0], hst_info_stand_by[i,1]), radius=DM/2., edgecolor = 'black', 
				 linewidth=0.5,facecolor = 'white')
				ax.add_patch(cell)
		#plt.grid(linestyle='--')
		
		#plt.xlim(-1900.,1900.)
		#plt.ylim(-1800.,2000.)
		plt.xlim(-2500.,2500.)
		plt.ylim(-2500.,2500.)
		ax.tick_params(axis='both', which='major', labelsize=20)
		ax.set_aspect(1)
		plt.savefig(self.casedir+'/field.png', bbox_inches='tight',dpi=100)
		plt.close('all')

		print('done')


