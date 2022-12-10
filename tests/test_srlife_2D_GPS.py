import time
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.family'] = 'Times'
import matplotlib.pyplot as plt
import os, sys, math, scipy.io, argparse
import gc
from tqdm import tqdm

sys.path.append('../..')

from srlife import receiver, solverparams, library, thermal, structural, system, damage, managers
from mdbapy.nitrateSaltHeuristics import *

def setup_problem(Ro, th, H_rec, Nr, Nt, Nz, times, fluid_temp, h_flux, pressure, T_base, folder=None, days=1):
	# Setup the base receiver
	period = 24.0                                                         # Loading cycle period, hours
	days = days                                                           # Number of cycles represented in the problem 
	panel_stiffness = "disconnect"                                        # Panels are disconnected from one another
	model = receiver.Receiver(period, days, panel_stiffness)              # Instatiating a receiver model

	# Setup each of the two panels
	tube_stiffness = "rigid"
	panel_0 = receiver.Panel(tube_stiffness)

	# Basic receiver geometry (Updated to Gemasolar)
	r_outer = Ro*1000                                                     # Panel tube outer radius (mm)
	thickness = th*1000                                                   # Panel tube thickness (mm)
	height = H_rec*1000                                                   # Panel tube height (mm)

	# Tube discretization
	nr = Nr                                                               # Number of radial elements in the panel tube cross-section
	nt = Nt                                                               # Number of circumferential elements in the panel tube cross-section
	nz = Nz                                                               # Number of axial elements in the panel tube

	# Setup Tube 0 in turn and assign it to the correct panel
	tube_0 = receiver.Tube(r_outer, thickness, height, nr, nt, nz, T0 = T_base)
	tube_0.set_times(times)
	tube_0.set_bc(receiver.ConvectiveBC(r_outer-thickness, height, nz, times, fluid_temp), "inner")
	tube_0.set_bc(receiver.HeatFluxBC(r_outer, height, nt, nz, times, h_flux), "outer")
	tube_0.set_pressure_bc(receiver.PressureBC(times, pressure))

	# Tube 1
	tube_1 = receiver.Tube(r_outer, thickness, height, nr, nt, nz, T0 = T_base)
	tube_1.set_times(times)
	tube_1.set_bc(receiver.ConvectiveBC(r_outer-thickness, height, nz, times, fluid_temp), "inner")
	tube_1.set_bc(receiver.HeatFluxBC(r_outer, height, nt, nz, times, h_flux), "outer")
	tube_1.set_pressure_bc(receiver.PressureBC(times, pressure))

	# Assign to panel 0
	panel_0.add_tube(tube_0, "tube0")
	panel_0.add_tube(tube_1, "tube1")

	# Assign the panels to the receiver
	model.add_panel(panel_0, "panel0")

	# Save the receiver to an HDF5 file
	if folder==None:
		fileName = 'model.hdf5'
	else:
		fileName = '%s/model.hdf5'%folder
	model.save('model.hdf5')

def run_problem(zpos,nz,progress_bar=True,folder=None,nthreads=4,load_state0=False,savestate=False,savefolder='.',loadfolder='.',debug=False):
	# Load the receiver we previously saved
	model = receiver.Receiver.load('model.hdf5')

	# Choose the material models
	fluid_mat = library.load_fluid("nitratesalt", "base")
	mat =     "A230"
	thermat = "base"                               # base
	defomat = "const_base"                         # base | elastic_creep | elastic_model | const_elastic_creep | const_base
	damat =   "base"                               # base
	thermal_mat, deformation_mat, damage_mat = library.load_material(mat, thermat, defomat, damat)

	# Cut down on run time for now by making the tube analyses 1D
	# This is not recommended for actual design evaluation
	for panel in model.panels.values():
		for tube in panel.tubes.values():
			tube.make_2D(tube.h/nz*zpos)
			tube.savefolder=savefolder
			tube.loadfolder=loadfolder
			tube.load_state0=load_state0
			tube.savestate=savestate

	# Setup some solver parameters
	params = solverparams.ParameterSet()
	params['progress_bars'] = progress_bar         # Print a progress bar to the screen as we solve
	params['nthreads'] = nthreads                  # Solve will run in multithreaded mode, set to number of available cores

	params["thermal"]["steady"] = False            # Ignore thermal mass and use conduction only
	params["thermal"]["rtol"] = 1.0e-6             # Iteration relative tolerance
	params["thermal"]["atol"] = 1.0e-4             # Iteration absolute tolerance
	params["thermal"]["miter"] = 20                # Maximum iterations
	params["thermal"]["substep"] = 1               # Divide user-provided time increments into smaller values

	params["structural"]["rtol"] = 1.0e-6          # Relative tolerance for NR iterations
	params["structural"]["atol"] = 1.0e-6          # Absolute tolerance for NR iterations
	params["structural"]["miter"] = 50             # Maximum newton-raphson iterations
	params["structural"]["verbose"] = False        # Verbose solve

	params["system"]["rtol"] = 1.0e-6              # Relative tolerance
	params["system"]["atol"] = 1.0e-6              # Absolute tolerance
	params["system"]["miter"] = 20                 # Number of permissible nonlinear iterations
	params["system"]["verbose"] = debug            # Print a lot of debug info

	# Choose the solvers, i.e. how we are going to solve the thermal,
	# single tube, structural system, and damage calculation problems.
	# Right now there is only one option for each
	# Define the thermal solver to use in solving the heat transfer problem
	thermal_solver = thermal.FiniteDifferenceImplicitThermalSolver(params["thermal"])
	# Define the structural solver to use in solving the individual tube problems
	structural_solver = structural.PythonTubeSolver(params["structural"])
	# Define the system solver to use in solving the coupled structural system
	system_solver = system.SpringSystemSolver(params["system"])
	# Damage model to use in calculating life
	damage_model = damage.TimeFractionInteractionDamage(params['damage'])

	# The solution manager
	solver = managers.SolutionManager(model, thermal_solver, thermal_mat, fluid_mat,
		structural_solver, deformation_mat, damage_mat,
		system_solver, damage_model, pset = params)

	# Actually solve for life
	solver.solve_heat_transfer()#
	solver.solve_structural()
#	mydict = {}
#	for tube in solver.tubes:
#		mydict["temperature"] = tube.results["temperature"]
#	scipy.io.savemat(os.path.join(savefolder,"quadrature_results.mat"),mydict)
	result = 1
	return result


class fmflow:
	def __init__(self, T_in, Tamb, CG, h_ext, model):
		self.Tamb = Tamb
		self.CG = CG
		self.h_ext = h_ext
		self.model = model
		self.Tf = T_in*np.ones(model.nz+1)
		self.qnet = np.zeros((2*model.nt-1,model.nz))

	def fun(self, m_flow):

		for k in range(self.model.nz):
			Qnet = self.model.Temperature(m_flow, self.Tf[k], self.Tamb, self.CG[k], self.h_ext)
			C = self.model.specificHeatCapacityCp(self.Tf[k])*m_flow
			if C>0:
				self.Tf[k+1] = self.Tf[k] + Qnet/C
			else:
				self.Tf[k+1] = self.Tf[k]
			self.qnet[:,k] = self.model.qnet/1e6
		res = abs(self.Tf[-1] - self.model.T_out)
		return res


def secant(fun, x0, tol=1e-2, maxiter=100, debug=False):
	x1 = (1.-1e-5)*x0
	f0 = fun(x0)
	f1 = fun(x1)
	i = 0
	xs = [x0,x1]
	fs = [f0,f1]
	alpha = 1
	while abs(f1) > tol and i < maxiter:
		x2 = float(f1 - f0)/(x1 - x0)
		x = x1 - alpha*float(f1)/x2

		if x<0:
			alpha=alpha/2.
		elif np.isnan(fun(x)):
			alpha=alpha/2.
		else:
			x0 = x1
			x1 = x
			f0 = f1
			f1 = fun(x1)
			alpha = 1
			xs.append(x1)
			fs.append(f1)
#			print(x1,f1)

	fmin = min(fs)
	imin = fs.index(fmin)
	xmin = xs[imin]
	return xmin

def run_gemasolar(panel,position,days,nthreads,clearSky,load_state0,savestate,step,debug,case="tmy",rfile = 'GemasolarSystemOperationCS_res.mat'):

	print(yellow('	Verification of inputs:'))
	print('	panel %s, pos %s, days %s-%s, nthreads=%s, clearSky=%s, load_state0=%s, savestate=%s, step:%s'%(panel,position,days[0],days[1],nthreads,clearSky,load_state0,savestate,step))

	# Material
	materials = {
		'A230':{'kp':17.97,'alpha':15.61e-06,'Young':184e9,'poisson':0.31},
		'800H':{'kp':18.30,'alpha':18.28e-06,'Young':171e9,'poisson':0.31}
	}

	DO=22.4;WT=1.2;mat="A230"

	# Instantiating model with the gemasolar geometry
	model = receiver_cyl(
		 Ri = DO/2000.0 - WT/1000.0
		,Ro = DO/2000.0
		,R_fouling=8.808e-5
		,alpha = materials[mat]['alpha']
		,Young = materials[mat]['Young']
		,poisson = materials[mat]['poisson']
		,ab = 0.93
		,em = 0.87
		,kp = materials[mat]['kp']
		,Dittus=False
		,mat=mat)
	# Number of radial nodes
	nr = 9

	model.import_mat(rfile)                                                                                 # Importing SolarTherm output
	times = model.data[:,0]                                                                                 # Get times
	CG = model.data[:,model._vars['CG[1]'][2]:model._vars['CG[450]'][2]+1]                                  # Get fluxes
	m_flow_tb = model.data[:,model._vars['m_flow_tb'][2]]                                                   # Get mass flow rates
	Tamb = model.data[:,model._vars['data.Tdry'][2]]                                                        # Get ambient temperatures
	h_ext = model.data[:,model._vars['h_conv'][2]]                                                          # Get heat transfer coefficient due to external convection
	con_state = model.data[:,model._vars['con_state'][2]]                                                   # Concentrator states

	# Filtering times
	index = []
	_times = []
	time_lb = days[0]*86400
	time_ub = days[1]*86400
	print(yellow('	Sorting times'))
	step_ramp=min(60.,step)
	for i in tqdm(range(len(times))):
		if con_state[i]==1:
			if times[i]%7200.==0 and times[i] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i)
				_times.append(times[i])
		elif con_state[i]==2 and con_state[i-1]==1:
			if times[i-1] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i-1)
				_times.append(times[i-1])
		elif con_state[i]==2 or con_state[i]==5:
			if times[i]%step_ramp==0 and times[i] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i)
				_times.append(times[i])
		elif con_state[i]==5 and con_state[i+1]==1:
			if times[i+1] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i+1)
				_times.append(times[i+1])
		else:
			if times[i]%step==0 and times[i] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i)
				_times.append(times[i])
	del model.data
	del _times
	gc.collect()

	times = times[index].flatten()/3600.                                                                  # Filtered times, in hours
	CG = CG[index,:]                                                                                      # Filtered fluxes
	m_flow_tb = m_flow_tb[index]                                                                          # Filtered mass flow rates
	Tamb = Tamb[index]                                                                                    # Filtered ambient temperature
	h_ext = h_ext[index]                                                                                  # Filtered heat transfer coefficient due to external convection
	con_state = con_state[index]                                                                          # Filtered concentrator states
	ramp_up = np.where(con_state==2)[0]
	full_ld = np.where(con_state==3)[0]
	part_ld = np.where(con_state==4)[0]
	ramp_dw = np.where(con_state==5)[0]

	# Instantiating variables
	Tf = 293.15*np.ones((times.shape[0],model.nz+1))
	Tf[full_ld,:] = model.T_in*np.ones_like(Tf[full_ld,:])
	Tf[part_ld,:] = model.T_in*np.ones_like(Tf[part_ld,:])
	Tf[ramp_dw,:] = 533.15*np.ones_like(Tf[ramp_dw,:])
	Tf[ramp_up,:] = 533.15*np.ones_like(Tf[ramp_up,:])
	qnet = np.zeros((times.shape[0],2*model.nt-2,model.nz))
	pressure = 0.1*np.ones_like(times)
	quarter=np.where(model.thetas>np.pi/2)[0]

	# Running thermal model
	print(yellow('	Running thermal model'))

	for k in tqdm(range(model.nz)):
		Qnet = model.Temperature(m_flow_tb, Tf[:,k], Tamb, CG[:,k], h_ext)
		C = model.specificHeatCapacityCp(Tf[:,k])*m_flow_tb
		Tf[:,k+1] = Tf[:,k] + np.divide(Qnet, C, out=np.zeros_like(C), where=C!=0)
		model.qnet[:,quarter]=0.0
		qnet[:,:,k] = np.concatenate((model.qnet[:,1:-1],model.qnet[:,::-1]),axis=1)/1e6

	# Getting internal pressure
	pressure = np.where(m_flow_tb>0, 0.1, m_flow_tb)
	pressure = pressure.flatten()
	lb = model.nbins*(panel-1)
	ub = lb + model.nbins

	ndays = (days[1]-days[0])

	loadfolder = os.path.join(os.getcwd(),'results_%s_d%s'%(case,days[0]))
	savefolder = os.path.join(os.getcwd(),'results_%s_d%s'%(case,days[1]))
	if not os.path.isdir(savefolder):
		os.mkdir(savefolder)

	Ro = model.Ro
	thickness = model.thickness
	H_rec = model.H_rec
	nt = 2*model.nt-2
	nbins = model.nbins

	# Creating the hdf5 model
	setup_problem(Ro,
	              thickness,
	              H_rec,
	              nr,
	              nt,
	              nbins,
	              times,
	              Tf[:,lb:ub],
	              qnet[:,:,lb:ub],
	              pressure,
	              T_base = 293.15,
	              days=ndays)

	# Running srlife
	if load_state0==1:
		load_state0=True
	life = run_problem(
		          position,
		          nbins,
		          nthreads=nthreads,
		          load_state0=load_state0,
		          savestate=True,
		          loadfolder=loadfolder,
		          savefolder=savefolder,
		          debug=debug)

	scipy.io.savemat('%s/inputs.mat'%savefolder,{
	              'times':times,
	              'qnet':qnet,
	              'Tf':Tf,
	              'pressure':pressure})

	# Plotting thermal results
	fig, axes = plt.subplots(2,3, figsize=(18,8))

	axes[0,0].plot(times, Tf[:,lb:ub])
	axes[0,0].set_ylabel(r'$T_\mathrm{f}$ [K]')
	axes[0,0].set_xlabel(r'$t$ [h]')

	axes[0,1].plot(times, qnet[:,0,lb:ub])
	axes[0,1].set_ylabel(r'$q^{\prime\prime}_\mathrm{net}$ [MW/m$^2$]')
	axes[0,1].set_xlabel(r'$t$ [h]')

	axes[0,2].plot(times, pressure)
	axes[0,2].set_ylabel(r'$P$ [MPa]')
	axes[0,2].set_xlabel(r'$t$ [h]')

	quadrature_results = scipy.io.loadmat('%s/quadrature_results.mat'%(savefolder))

	vm = np.sqrt((
	              (quadrature_results['stress_xx'] - quadrature_results['stress_yy'])**2.0 + 
	              (quadrature_results['stress_yy'] - quadrature_results['stress_zz'])**2.0 + 
	              (quadrature_results['stress_zz'] - quadrature_results['stress_xx'])**2.0 + 
	              6.0 * (quadrature_results['stress_xy']**2.0 + 
	              quadrature_results['stress_yz']**2.0 + 
	              quadrature_results['stress_xz']**2.0))/2.0)

	em = np.sqrt((
	              (quadrature_results['mechanical_strain_xx'] - quadrature_results['mechanical_strain_yy'])**2.0 + 
	              (quadrature_results['mechanical_strain_yy'] - quadrature_results['mechanical_strain_zz'])**2.0 + 
	              (quadrature_results['mechanical_strain_zz'] - quadrature_results['mechanical_strain_xx'])**2.0 + 
	              6.0 * (quadrature_results['mechanical_strain_xy']**2.0 + 
	              quadrature_results['mechanical_strain_yz']**2.0 + 
	              quadrature_results['mechanical_strain_xz']**2.0))/2.0)

	axes[1,0].plot(times,vm[:, 0,0],label='Inner')
	axes[1,0].plot(times,vm[:,-1,0],label='Outer')
	axes[1,0].set_xlabel(r'$t$ [h]')
	axes[1,0].set_ylabel(r'$\sigma_\mathrm{crown,eq}$ [MPa]')
	axes[1,0].legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

	axes[1,1].plot(times,quadrature_results['temperature'][:, 0,0]-273.15,label='Inner')
	axes[1,1].plot(times,quadrature_results['temperature'][:,-1,0]-273.15,label='Outer')
	axes[1,1].set_xlabel(r'$t$ [h]')
	axes[1,1].set_ylabel(r'$T_\mathrm{crown}$ [\textdegree C]')
	axes[1,1].set_ylim([-0.05,700])
	axes[1,1].legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

	axes[1,2].plot(times,em[:, 0,0],label='Inner')
	axes[1,2].plot(times,em[:,-1,0],label='Outer')
	axes[1,2].set_xlabel(r'$t$ [h]')
	axes[1,2].set_ylabel(r'$\epsilon_\mathrm{crown,eq}$ [mm/mm]')
	axes[1,2].legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

	plt.tight_layout()
	plt.savefig('%s/results_%s_d%s'%(savefolder,case,days[1]))
	plt.close(fig)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Estimates average damage of a representative tube in a receiver panel')
	parser.add_argument('--panel', type=int, default=7, help='Panel to be simulated. Default=1')
	parser.add_argument('--position', type=float, default=26, help='Panel position to be simulated. Default=1')
	parser.add_argument('--days', nargs=2, type=int, default=[0,1], help='domain of days to simulate')
	parser.add_argument('--nthreads', type=int, default=4, help='Number of processors. Default=4')
	parser.add_argument('--clearSky', type=bool, default=False, help='Run clear sky DNI (requires to have the solartherm results)')
	parser.add_argument('--load_state0', type=int, default=0, help='Load state from a previous simulation')
	parser.add_argument('--savestate', type=bool, default=True, help='Save the last state of the last simulated day')
	parser.add_argument('--step', type=float, default=900, help='Simulation step. Default=900')
	parser.add_argument('--debug', type=bool, default=False, help='Debug option. Default=False')
	parser.add_argument('--rfile', type=str, default='/home/arfontalvo/solartherm/examples/SimpleSystemOperation_res_1m_new.mat')
	args = parser.parse_args()

	tinit = time.time()
	lb=args.days[0]
	ub=args.days[1]
	for i in range(lb,ub):
		args.days=[i,i+1]
		try:
			run_gemasolar(args.panel,args.position,args.days,args.nthreads,args.clearSky,args.load_state0,args.savestate,args.step,args.debug,rfile = args.rfile)
		except RuntimeError:
			try:
				time_step = 600.0
				run_gemasolar(args.panel,args.position,args.days,args.nthreads,args.clearSky,args.load_state0,args.savestate,time_step,args.debug,rfile = args.rfile)
			except RuntimeError:
				time_step = 60.0
				run_gemasolar(args.panel,args.position,args.days,args.nthreads,args.clearSky,args.load_state0,args.savestate,time_step,args.debug,rfile = args.rfile)
	seconds = time.time() - tinit
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
