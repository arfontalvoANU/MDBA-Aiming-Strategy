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
import time

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

class gemasolar:
	def __init__(self,
		num_bundle=18, r_height=10.5, r_diameter=8.5, bins=50,
		tower_h=114.75, azimuth=0.0, elevation=52.38, DNI=950.0,
		num_fp=2, D0=22.40, num_rays=int(5e6), latitude=37.56,
		ndec=5, nhra=25, abs_t=0.94, ems_t=0.88,
		testcase='case-gemasolar',new_algorithm=False):

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
		self.ndec=ndec
		self.nhra=nhra
		self.abs_t=abs_t
		self.ems_t=ems_t
		# Declination angle
		self.Dec=np.linspace(-23.45, 23.45, num=self.ndec)
		# Solar hour angle
		self.Hra=np.linspace(-180.0, 180.0, num=self.nhra)
		self.N = int(self.num_bundle/self.num_fp*self.bins)
		self.DNI_ratio=[1.39,1.00,0.87,0.56]
		self.N_ratios = len(self.DNI_ratio)
		self.A_helio = 12.305*9.752

	def design_point(self,ratio,C_start,E_start,A_start):

		# define a unique case folder for the user
		self.basefolder = os.path.join(currentdir,'%s/DNI_ratio_%s'%(self.testcase,ratio))

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

		results,aiming_results,Strt=Model.HT_model(20.,0.)
		hst_info = np.loadtxt(Model.csv_trimmed, delimiter=',', skiprows=2)
		self.nhelios = int(len(hst_info[:,0]))
		self.eff_rec_design = results[9]

	def simulate_dni_ratio(self,C_start,E_start,A_start):

		for ratio in range(self.N_ratios):
			tinit = time.time()
			# define a unique case folder for the user
			self.basefolder = os.path.join(currentdir,'%s/DNI_ratio_%s'%(self.testcase,ratio))

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
			E = np.zeros((self.ndec, self.nhra))
			sun=SunPosition()
			irow = 0
			icol = 0

			self.sunpos = []
			self.irow = []
			self.icol = []
			# Estimating the azimuth and elevation vectors
			res=0
			f = open('%s/OELT_verification.csv'%(self.basefolder),'a')
			f.write('res,dec,hra,azi,ele,eta,sucess\n')
			for irow,dec in enumerate(self.Dec):
				for icol,hra in enumerate(self.Hra):
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
						DNI = self.DNI_ratio[ratio]*Model.get_I_Meinel(ele)
						if res >= 0:
							self.sunpos.append(res)
							self.irow.append(irow)
							self.icol.append(icol)
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
								self.D0, lat=self.latitude)
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
			np.savetxt('%s/OELT.txt'%(self.basefolder),E,fmt='%s', delimiter=',')
			del E

			# Print elapsed time
			seconds = time.time() - tinit
			m, s = divmod(seconds, 60)
			h, m = divmod(m, 60)
			print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))

			for fpath in range(1,3):
				# Writting the flux lookup table for this DNI ratio
				f = open('%s/flux_a230_salt_FP%s_DNIr%s.motab'%(os.path.join(currentdir,self.testcase),fpath,ratio),'w+')
				f.write('#1\n')
				E = np.zeros((self.ndec, self.nhra))
				F = []
				for k in range(self.N):
					f.write('double flux_%s(%d,%d)\n'%(k+1, self.ndec+1, self.nhra+1))
					f.write('0.0,')
					for res,irow,icol in zip(self.sunpos,self.irow,self.icol):
						filename = '%s/sunpos_%s/flux-table'%(self.basefolder, res)
						fileo = open(filename,'r')
						data = pickle.load(fileo)
						fileo.close()
						areas = data['areas']
						q_net = data['q_net']
						fp = data['fp']
						flux = q_net[fp[fpath-1]]/areas[fp[fpath-1]]/1e3
						E[irow,icol] = flux[k]
					x = np.linspace(-180.0, 180.0, self.nhra)
					y = np.linspace(-23.45, 23.45, self.ndec)
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

	def get_m_flow_lookup_tables(self):
		for fpath in range(1,3):
			f = open('%s/mflow_a230_salt_FP%s.motab'%(os.path.join(currentdir,self.testcase),fpath),'w+')
			f.write('#1\n')
			for ratio in range(self.N_ratios):
				E = np.zeros((self.ndec, self.nhra))
				f.write('double mflow_%s(%d,%d)\n'%(ratio, self.ndec+1, self.nhra+1))
				f.write('0.0,')
				for res,irow,icol in zip(self.sunpos, self.irow, self.icol):
					filename = '%s/DNI_ratio_%s/sunpos_%s/flux-table'%(os.path.join(currentdir,self.testcase), ratio, res)
					fileo = open(filename,'r')
					data = pickle.load(fileo)
					fileo.close()
					mflow = data['m']
					n_tubes = data['n_tubes']
					E[irow,icol] = mflow[int(fpath)-1]/n_tubes[0]
				x = np.linspace(-180.0, 180.0, self.nhra)
				y = np.linspace(-23.45, 23.45, self.ndec)
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

	def RLLT_ouput(self,C_start,E_start,A_start):
		tinit = time.time()
		# define a unique case folder for the user
		self.basefolder = os.path.join(currentdir,'%s/RLLT'%self.testcase)

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

		# Creating RLLT text file
		if not os.path.isfile('%s/RLLT.csv'%self.basefolder):
			f = open('%s/RLLT.csv'%self.basefolder,'a')
			header = \
				'row, h_conv_ext (W/m2/K), T_avg (K), (T**4_avg)**0.25 (K),'+\
				'conv_loss (W), rad_loss (W), receiver_loss (W),'+\
				'receiver_output (W), eta_th,C_aiming\n'
			f.write(header)
			f.close()

		# Importing RLLT inputs
		data = np.genfromtxt('%s/RLLT_input_gemasolar.csv'%currentdir, delimiter=',', skip_header=1)
		nrows = data.shape[0]
		for i in range(nrows):
			azimuth = data[i,4]
			elevation = data[i,5]
			dni = data[i,6]
			Tamb = data[i,7]
			Wspd = data[i,8]
			print(yellow('row: %s of %s\tdni: %s [W/m2]\tTamb: %s [deg C]\tWspd: %s [m/s]'%(i+1,nrows,dni,Tamb,Wspd)))
			casefolder = '%s/row_%s'%(self.basefolder,i)
			if not os.path.exists(casefolder):
				os.makedirs(casefolder)

			# Re-aiming
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
				dni,
				self.D0, lat=self.latitude)
			# Running aiming for design point
			Model.sweeping_algorithm(C_start,E_start,A_start)

			# Optical postprocessing
			eta,q_results,eta_exc_intec=proces_raw_results(
					'%s/vtk/simul'%casefolder,
					casefolder)

			results,aiming_results,Strt=Model.HT_model(Tamb,Wspd)
			f = open('%s/RLLT.csv'%self.basefolder,'a')
			f.write('%s,'%i)                                #row
			f.write('%s,'%results[5])                       #h_conv_ext
			f.write('%s,'%results[3])                       #T_avg
			f.write('%s,'%results[4])                       #(T**4_avg)**0.25
			f.write('%s,'%(results[8]*1e6))                 #conv_loss
			f.write('%s,'%(results[7]*1e6))                 #rad_loss
			f.write('%s,'%((results[7]+results[8])*1e6))    #receiver_loss
			f.write('%s,'%(results[0]*results[9]*1e6))      #receiver_output
			f.write('%s,'%results[9])                       #receiver_output
			f.write('%s\n'%Model.C_aiming[0])               #aiming extent
			f.close()

		# Print elapsed time
		seconds = time.time() - tinit
		m, s = divmod(seconds, 60)
		h, m = divmod(m, 60)
		print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))

	def coefficients(self,ratio,C_start,E_start,A_start):
		from pandas import DataFrame
		from sklearn import linear_model

		fpath = os.path.join(currentdir,'%s/RLLT'%self.testcase)

		# read from input file
		eff_abs=self.abs_t/(2./np.pi*(1.-self.abs_t)+self.abs_t) # effective absorptivity of the pipe
		eff_emi=self.ems_t/(2./np.pi*(1.-self.ems_t)+self.ems_t) # effective emissivity of the pipe

		R_input=np.loadtxt('RLLT_input.csv',delimiter=',', skiprows=1)
		R_output_old=np.loadtxt('%s/RLLT.csv' % fpath   ,delimiter=',', skiprows=1)
		R_output=R_output_old[~np.isnan(R_output_old).any(axis=1)] # to remove the rows with nan
		R_input=R_input[~np.isnan(R_output_old).any(axis=1)]

		R_input=R_input[~np.isnan(R_output).any(axis=1)]
		Q_in=R_output[:,-2]/R_output[:,-1]/eff_abs/1e6  # MW
		T_amb=R_input[:,-2] # C
		V_wind=R_input[:,-1] # m/s
		R_eff=R_output[:,-1]*eff_abs

		A_rec=R_output[0,3]/(R_output[0,0]*(R_output[0,1]-(R_input[0,-2]+273.15))) # receiver surface area

		Q_ref = (1-eff_abs)*Q_in           # the reflecive loss

		# linear regression for T_ext
		T_ext={'X2': Q_in, 'X3': (T_amb+273.15),'X4': V_wind,'T_ext': R_output[:,1]}
		df = DataFrame(T_ext,columns=['X2','X3','X4','T_ext'])
		X = df[['X2','X3','X4']] 
		Y = df['T_ext']
		regr = linear_model.LinearRegression()
		regr.fit(X, Y)
		C0=regr.intercept_
		C=regr.coef_
		T_ext_linear=C0+C[0]*Q_in+C[1]*(T_amb+273.15)+C[2]*V_wind

		# linear regression for T_ext_4_mean
		T_ext_4={'X2': Q_in, 'X3': (T_amb+273.15),'X4': V_wind,'T_ext_4': R_output[:,2]}
		df = DataFrame(T_ext_4,columns=['X2','X3','X4','T_ext_4'])
		X = df[['X2','X3','X4']] 
		Y = df['T_ext_4']
		regr = linear_model.LinearRegression()
		regr.fit(X, Y)
		C1=regr.intercept_
		C1_1=regr.coef_
		T_ext_4_linear=C1+C1_1[0]*Q_in+C1_1[1]*(T_amb+273.15)+C1_1[2]*V_wind

		coefs_T=[C0,C,C1,C1_1]

		# h with polynominal fitting
		coefs=np.polyfit(R_input[:,-1],R_output[:,0],4)
		h_conv=coefs[4]+coefs[3]*R_input[:,-1]+coefs[2]*(R_input[:,-1])**2+coefs[1]*(R_input[:,-1])**3+coefs[0]*(R_input[:,-1])**4

		Q_emi=eff_emi*5.67e-8*A_rec*(T_ext_4_linear**4-(T_amb+273.15)**4)/1e6

		Q_conv=h_conv*A_rec*(T_ext_linear-(T_amb+273.15))/1e6

		Qnet=eff_abs*Q_in-Q_conv-Q_emi
		Qnet_real= R_output[:,-2]/1e6
		Diff=(Qnet-Qnet_real)/Qnet_real


		data = np.loadtxt('%s/DNI_ratio_%s/OELT.txt'%(self.testcase,ratio), delimiter=',')
		self.design_point(ratio,C_start,E_start,A_start)
		# Extracting declination and hour angles
		nhra = data.shape[1]
		ndec = data.shape[0]
		dec = np.linspace(-23.45,23.45,ndec)
		hra = np.linspace(-180.0,180.0,nhra)
		eff_design = data[ndec-1,(nhra-1)/2+1]
		# Creating motab text file
		table = open('%s/gemasolar_H230_salt_MDBA_r%s.motab'%(os.path.join(currentdir,self.testcase),ratio),'w+')
		# Writting headers
		table.write('#1\n')
		table.write('#Comments\n')
		table.write("#METALABELS,n_helios,A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,CoT1,CoT2,CoT3,CoT4,CoT'1,CoT'2,CoT'3,CoT'4,Coh1,Coh2,Coh3,Coh4,Coh5,eff_abs,eff_emi,SM\n")
		table.write('#METAUNITS,integer,m2,real,m,m,m,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real\n')
		table.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(
			self.nhelios, self.A_helio, eff_design, self.r_diameter, self.r_height, self.tower_h, self.eff_rec_design,
			coefs_T[0], coefs_T[1][0], coefs_T[1][1], coefs_T[1][2], coefs_T[2], coefs_T[3][0], coefs_T[3][1], coefs_T[3][2],
			coefs[0], coefs[1], coefs[2], coefs[3], coefs[4],
			eff_abs, eff_emi,self.nhelios*self.A_helio*eff_design*self.DNI*self.eff_rec_design/(19.9e6/0.3774)
			)
			)
		table.write('double optics(%s,%s)\n'%(ndec+1,nhra+1))
		table.write('0.0,')
		for j in range(nhra):
			if j==nhra-1:
				table.write('%s\n'%hra[j])
			else:
				table.write('%s,'% hra[j])

		# writting values
		for i in range(ndec):
			table.write('%s,'% dec[i])
			for j in range(nhra):
				if j==nhra-1:
					table.write('%s\n'%data[i,j])
				else:
					table.write('%s,'%data[i,j])
		table.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Run representative sun positions of annual simulation for a specific DNI ratio')
	parser.add_argument('--ratio', type=int, default=1, help='DNI ratio to be simulated. Default=1')
	parser.add_argument('--fpath', type=int, default=1, help='Flowpath data to be compiled. Default=1')
	parser.add_argument('--C_start', type=float, default=0.5, help='The starting value of the aiming extent')
	parser.add_argument('--E_start', type=float, default=2.0, help='The starting value of the aiming exponent')
	parser.add_argument('--A_start', type=float, default=0.5, help='The starting value of the aiming asymetry factor')
	parser.add_argument('--design', type=bool, default=False, help='Run only the design point')
	parser.add_argument('--new_algorithm', type=bool, default=False, help='Run the new search algorithm')
	args = parser.parse_args()

	cyl_receiver = gemasolar(new_algorithm=args.new_algorithm)
	if args.design:
		cyl_receiver.design_point(args.ratio, args.C_start, args.E_start, args.A_start)
	else:
		cyl_receiver.simulate_dni_ratio(args.C_start, args.E_start, args.A_start)
		cyl_receiver.get_m_flow_lookup_tables()
		cyl_receiver.RLLT_ouput(args.C_start, args.E_start, args.A_start)
		cyl_receiver.coefficients(args.ratio, args.C_start, args.E_start, args.A_start)
