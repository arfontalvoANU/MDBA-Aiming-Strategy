#!/usr/bin/env python3

import sys, time, os
from math import exp, log, sqrt, pi, ceil, floor, asin

import numpy as np
import scipy.optimize as opt

import nashTubeStress as nts
from nashTubeStress import coolant
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm

LIBRARY_DIR = os.path.abspath(
	os.path.join(os.path.dirname(nts.__file__), 'mats'))

################################### PLOTTING ###################################

import matplotlib as mpl
import matplotlib.pyplot as plt
## uncomment following if not in ~/.config/matplotlib/matplotlibrc already
plt.rc('text', usetex = True)
plt.rcParams['text.latex.preamble'] = r"\usepackage{newtxtext,newtxmath,siunitx}"
mpl.rc('figure.subplot', bottom=0.13, top=0.95)
mpl.rc('figure.subplot', left=0.15, right=0.95)
from matplotlib import colors, ticker, cm
from matplotlib.transforms import Affine2D
from matplotlib.lines import Line2D
from matplotlib.projections import PolarAxes
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
from mpl_toolkits.axisartist.grid_finder import \
	(FixedLocator, MaxNLocator, DictFormatter)

################################### FUNCTIONS ##################################
def findFlux(flux, s, f, i, point):
	"""
	Helper for finding optimum flux for certain stress condition
	"""
	s.CG = flux

	ret = s.solve(eps=1e-6)
	s.postProcessing()

	if point=='max':
		# T_max|sigmaEqMax:
		sigmaEqMax = np.interp(np.max(s.T), f[:,0], f[:,i])
		return sigmaEqMax - np.max(s.sigmaEq)
	elif point=='inside':
		# T_i:
		sigmaEqMax = np.interp(s.T[0,0], f[:,0], f[:,i])
		return sigmaEqMax - s.sigmaEq[0,0]
	elif point=='outside':
		# T_o
		sigmaEqMax = np.interp(s.T[0,-1], f[:,0], f[:,i])
		return sigmaEqMax - s.sigmaEq[0,-1]
	elif point=='membrane':
		# Assuming shakedown has occured (membrane stress remains):
		sigmaEqMax = np.interp(np.average(s.T[0,:]), f[:,0], f[:,i])
		return sigmaEqMax - np.average(s.sigmaEq[0,:])
	elif point=='corrosion':
		# T_i:
		return 630 + 273.15 - s.T[0,0]
	else: sys.exit('Variable point {} not recognised'.format(point))

##################################### MAIN #####################################

def fluxLim(OD,WT,fname,mat):
	if not os.path.exists('{0}_OD{1:.2f}_WT{2:.2f}_peakFlux.csv'.format(fname, OD, WT)):
		h_ext=10. # convective loss due to wind W/(m2.K)
		salt = coolant.nitrateSalt(False)
		salt.update(723.15)
		iterator='inline'
		nr=6; nt=61

		T_int = np.linspace(290, 565, 12)
		T_int = np.append(T_int,600) + 273.15

		""" Instantiating figure and subplots"""
		fig = plt.figure(figsize=(3.5, 3.5))
		ax = fig.add_subplot(111)
		csv = np.c_[T_int,]
		header  = '0'

		b = OD/2e3		 # outside tube radius [mm->m]
		a = (b-WT*1e-3)	 # inside tube radius [mm->m]
		g = nts.Grid(nr=nr, nt=nt, rMin=a, rMax=b) # nr, nt -> resolution
		vfs = np.array([0, 1., 2., 3., 4.])

		""" Looping velocities"""
		for vf in vfs[1:]:
			props = np.genfromtxt(os.path.join(LIBRARY_DIR, 'props/{0}'.format(mat)), delimiter=';')
			props[:,0] += 273.15 # degC to K
			props[:,1] *= 1e9  # convert GPa -> Pa
			props[:,3] *= 1e-6 # convert 1e-6 mm/mm -> mm/mm
			E     = np.interp(723.15,props[:,0],props[:,1])
			k     = np.interp(723.15,props[:,0],props[:,2])
			alpha = np.interp(723.15,props[:,0],props[:,3])
			nu = 0.31
			s = nts.Solver(g, debug=False, CG=0.85e6, k=k, T_int=723.15, R_f=8.808e-5,
							A=0.968, epsilon=0.87, T_ext=293.15, h_ext=h_ext,
							P_i=0e5, alpha=alpha, E=E, nu=nu, n=1,
							bend=False)
			s.extBC = s.extTubeHalfCosFluxRadConv
			s.intBC = s.intTubeConv
			fv = np.genfromtxt(os.path.join(LIBRARY_DIR, mat), delimiter=';')
			fv[:,0] += 273.15 # degC to K
			fv[:,2] *= 3e6 # apply 3f criteria to Sm and convert MPa->Pa
			TSod_met = np.zeros(len(T_int))
			fluxSalt = np.zeros(len(T_int))
			print('vf: {0:.2f} m/s'.format(vf))
			for i in tqdm(range(len(T_int))):
				s.E =     np.interp(T_int[i],props[:,0],props[:,1])
				s.k =     np.interp(T_int[i],props[:,0],props[:,2])
				s.alpha = np.interp(T_int[i],props[:,0],props[:,3])
				s.T_int = T_int[i]
				salt.update(T_int[i])
				s.h_int, dP = coolant.HTC(False, salt, a, b, s.k, 'Gnielinski', 'velocity', vf)
				fluxSalt[i] = opt.newton(
					findFlux, 1e5,
					args=(s, fv, 2, 'outside'),
					maxiter=1000, tol=1e-2
				)
				TSod_met[i] = np.max(s.T)

			ax.plot(T_int-273.15,fluxSalt*1e-6, label=r'U = {0} m/s'.format(vf))
			ax.set_xlabel(r'\textsc{fluid temperature}, '+\
						  '$T_\mathrm{f}$ (\si{\celsius})')
			ax.set_ylabel(
				r'\textsc{incident flux}, $\vec{\phi_\mathrm{q}}$ '+\
				'(\si{\mega\watt\per\meter\squared})'
			)
			ax.legend(loc='best')
			fig.tight_layout()
			fig.savefig('{0}_OD{1:.2f}_WT{2:.2f}_peakFlux.pdf'.format(fname, OD, WT),
						transparent=True)
			plt.close(fig)
			## Dump peak flux results to CSV file:
			csv  = np.c_[csv, fluxSalt,]

		csv = np.concatenate((vfs.reshape(1,len(vfs)),csv))
		np.savetxt('{0}_OD{1:.2f}_WT{2:.2f}_peakFlux.csv'.format(fname, OD, WT), csv, fmt='%s', delimiter=',')
