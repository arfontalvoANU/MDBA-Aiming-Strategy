import numpy as np
import os,shutil
import re
import solsticepy
from python_postprocessing import *
from solsticepy.master import Master
import colorama
colorama.init()

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def simul2data(simulfile,fluxfile,r_height,r_diameter,rec_slices,rec_stacks):
	rec_slices = int(rec_slices)
	rec_stacks = int(rec_stacks)
	elem_area=r_diameter*np.sin(np.pi/rec_slices)*r_height/rec_stacks
	simul=simulfile + '/simul'

	with open(simul) as f:
		lines = f.readlines()

	receiver_name='Front_faces_Absorbed_flux'
	for line in lines:
		if receiver_name in line:
			k=lines.index(line)

	for line in lines:
		if 'Back_faces_Incoming_flux' in line:
			kk=lines.index(line)

	lines_loc=lines[k:kk]

	for line in lines_loc:
		if 'LOOKUP_TABLE' in line:
			nn=lines_loc.index(line)+2
			break

	FluxTriangle=np.array([])
	for j in range(nn,nn+2*rec_stacks*rec_slices):
		FluxTriangle=np.append(FluxTriangle,(float(lines_loc[j].rstrip().split(' ', 1 )[0]))/1000.)

	Flux = np.zeros([rec_stacks,rec_slices])
	f = open(fluxfile,'w+')
	for x in range(rec_stacks):
		for y in range(rec_slices):
			Flux[x,y]=0.5*(FluxTriangle[2*rec_stacks*y+2*x]+FluxTriangle[2*rec_stacks*y+2*x+1])
			if y == (rec_slices - 1):
				f.write(str(Flux[x,y]) + '\n')
			else:
				f.write(str(Flux[x,y]) + ',')
	f.close()
	F_ave = np.sum(Flux)/np.size(Flux)
	return F_ave,Flux

class optical:
	def __init__(self,
					Dec,
					Hra,
					azi,
					ele,
					rec_slices=16.,
					rec_stacks=50.,
					rec_diameter=16.,
					rec_height=24.,
					DNI=980.,
					sunshape='pillbox',
					half_angle_deg = 0.2664,
					csr = 0.01,
					num_rays=5e6,
					layoutfile='demo_layout_and_aiming.csv',
					rec_pos_x=0.,
					rec_pos_y=0.,
					rec_pos_z=187.,
					rec_abs=0.98,
					receiver='cylinder',
					hemisphere='North',
					hst_w=12.2,
					hst_h=12.2,
					rho_refl=0.90,
					slope_error=1.5e-3,
					tower_h=187.,
					tower_r=0.01
					):

		self.rec_slices = rec_slices                 # Number of elements in the circumferetial direction
		self.rec_stacks = rec_stacks                 # Number of elements in the vertical direction
		self.hemisphere = hemisphere                 # North or South hemisphere

		# The sun
		# =========
		self.sun = solsticepy.Sun(dni=DNI, sunshape=sunshape, half_angle_deg=half_angle_deg, csr=csr)

		# S1. number of rays for the ray-tracing simulation
		self.num_rays=num_rays

		self.Azi = azi                               # Vector of solar azimuth angles from discretisation
		self.Ele = ele                               # Vector of solar elevation angles from discretisation
		self.Dec = Dec
		self.Hra = Hra
		self.efficiency = []

		# The field
		# ==========
		# F1. Heliostat
		self.hst_w=hst_w                             # Heliostat width [m]
		self.hst_h=hst_h                             # Heliostat height [m]
		self.rho_refl=rho_refl                       # mirror reflectivity
		self.slope_error=slope_error                 # Slope error [radians]

		# F2. Tower
		self.tower_h=tower_h                         # tower height
		self.tower_r=tower_r                         # tower radius

		# The receiver
		# ============
		# R1. shape
		self.receiver=receiver                       # Receiver type [str]: 'flat' or 'cylinder'
		# R2. Size
		self.rec_r=rec_diameter                      # Receiver diameter [m]
		self.rec_h=rec_height                        # Receiver height [m]
		# R3. position
		self.loc_x=rec_pos_x                         # Receiver position along the x axis [m]
		self.loc_y=rec_pos_y                         # Receiver position along the y axis [m]
		self.loc_z=rec_pos_z                         # Receiver position along the z axis [m]
		# R4. Abosrptivity
		self.rec_abs=rec_abs

		# Extracting the heliostat positions from the loaded CSV file
		# ===========
		layout=np.loadtxt(layoutfile, delimiter=',', skiprows=2)
		self.hst_pos=layout[:,:3]
		self.hst_foc=layout[:,3] 
		self.hst_aims=layout[:,4:]
		self.one_heliostat=False
		self.rec_param=np.array([
								self.rec_r, 
								self.rec_h, 
								self.rec_slices,
								self.rec_stacks, 
								self.loc_x, 
								self.loc_y, 
								self.loc_z])

	def sim_group(self,
			userdefinedfolder=False,
			folder='case-Thu-13-48'
			):

		if userdefinedfolder:
			basefolder=folder
		else:
			# define a unique case folder for the user
			snum = 0
			suffix = ""
			while 1:
				import datetime
				dt = datetime.datetime.now()
				ds = dt.strftime("%a-%H-%M")
				npos = np.size(self.Azi)
				basefolder = os.path.join(os.getcwd(),'case-%s%s-npos%s'%(ds,suffix,npos))
				if os.path.exists(basefolder):
					snum+=1
					suffix = "-%d"%(snum,)
					if snum > 200:
						raise RuntimeError("Some problem with creating basefolder")
				else:
					# good, we have a new case dir
					break
		res = 0
		# We simulate only morning to solar noon 
		# and mirror the optical efficiencies for the afternoon
		E = np.zeros((len(self.Dec), 2*len(self.Hra)-1))
		irow = 0
		icol = 0
		if not os.path.exists(basefolder):
			for azimuth,elevation in zip(self.Azi,self.Ele):

				# We don't simulate elevations below the threshold
				if elevation > 0.0:
					casefolder = os.path.join(basefolder,'./design_pos_%s'%(res))
					if os.path.exists(casefolder):
						shutil.rmtree(casefolder)

					master=Master(casedir=casefolder)
					outfile_yaml = master.in_case(folder=casefolder, fn='input.yaml')
					outfile_recv = master.in_case(folder=casefolder, fn='input-rcv.yaml')

					# generate the YAML file from the input parameters specified above
					solsticepy.gen_yaml(
								self.sun, 
								self.hst_pos, 
								self.hst_foc, 
								self.hst_aims,
								self.hst_w, 
								self.hst_h,
								self.rho_refl, 
								self.slope_error, 
								self.receiver, 
								self.rec_param, 
								self.rec_abs, 
								outfile_yaml=outfile_yaml, 
								outfile_recv=outfile_recv, 
								hemisphere=self.hemisphere, 
								tower_h=self.tower_h, 
								tower_r=self.tower_r, 
								spectral=False, 
								medium=0, 
								one_heliostat=self.one_heliostat)

					# run Solstice using the generate inputs, and run all required post-processing
					master.run(
								azimuth, 
								elevation, 
								int(self.num_rays), 
								self.rho_refl,
								self.sun.dni, 
								folder=casefolder, 
								gen_vtk=False, 
								verbose=True
								)

					efficiency_total, performance_hst = solsticepy.process_raw_results(
								'%s/simul'%casefolder,
								'%s/simul'%casefolder,
								self.rho_refl,
								self.sun.dni)

					print('optical efficiency: %s'%(yellow(str(efficiency_total))))
					self.efficiency.append(efficiency_total.nominal_value)
					E[irow,icol] = efficiency_total.nominal_value

				res += 1
				if icol+1==len(self.Hra):
					icol  =0
					irow +=1
				else:
					icol +=1

		# Writting outputs to OELT file
		E[:,len(self.Hra):2*len(self.Hra)-1]=np.fliplr(E)[:,len(self.Hra):2*len(self.Hra)-1]
		np.savetxt('%s/OELT.csv'%basefolder,E,fmt='%s', delimiter=',')
		self.basefolder = basefolder

	def sim_single(self,
				azimuth=0.0,
				elevation=45.0,
				userdefinedfolder=False,
				folder='case-Thu-13-48'
				):

		if userdefinedfolder:
			basefolder=folder
		else:
			# define a unique case folder for the user
			snum = 0
			suffix = ""
			while 1:
				import datetime
				dt = datetime.datetime.now()
				ds = dt.strftime("%a-%H-%M")
				npos = np.size(self.Azi)
				basefolder = os.path.join(os.getcwd(),'case-%s%s-npos%s'%(ds,suffix,npos))
				if os.path.exists(basefolder):
					snum+=1
					suffix = "-%d"%(snum,)
					if snum > 200:
						raise RuntimeError("Some problem with creating basefolder")
				else:
					# good, we have a new case dir
					break
		if not os.path.exists(basefolder):
			# We don't simulate elevations below the threshold
			if elevation > 0.0:
				casefolder = os.path.join(basefolder,'./design_pos')
				if os.path.exists(casefolder):
					shutil.rmtree(casefolder)

				master=Master(casedir=casefolder)
				outfile_yaml = master.in_case(folder=casefolder, fn='input.yaml')
				outfile_recv = master.in_case(folder=casefolder, fn='input-rcv.yaml')

				# generate the YAML file from the input parameters specified above
				solsticepy.gen_yaml(
							self.sun, 
							self.hst_pos, 
							self.hst_foc, 
							self.hst_aims,
							self.hst_w, 
							self.hst_h,
							self.rho_refl, 
							self.slope_error, 
							self.receiver, 
							self.rec_param, 
							self.rec_abs, 
							outfile_yaml=outfile_yaml, 
							outfile_recv=outfile_recv, 
							hemisphere=self.hemisphere, 
							tower_h=self.tower_h, 
							tower_r=self.tower_r, 
							spectral=False, 
							medium=0, 
							one_heliostat=self.one_heliostat)

				# run Solstice using the generate inputs, and run all required post-processing
				master.run(
							azimuth, 
							elevation, 
							int(self.num_rays), 
							self.rho_refl,
							self.sun.dni, 
							folder=casefolder, 
							gen_vtk=False, 
							verbose=True
							)

				simul2data(
							casefolder,
							'%s/flux-table.csv'%casefolder,
							self.rec_h,
							self.rec_r,
							self.rec_slices,
							self.rec_stacks
							)

		self.basefolder = basefolder


