#!/usr/bin/env python3

import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.family'] = 'Times'
import os, sys, math, scipy.io, argparse
from scipy.interpolate import interp1d, RegularGridInterpolator
import scipy.optimize as opt
import time, ctypes
from numpy.ctypeslib import ndpointer
from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import warnings
import pickle

import colorama
colorama.init()
def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

sys.path.append('../..')

from srlife import receiver, solverparams, library, thermal, structural, system, damage, managers
from neml import uniaxial

strMatNormal = lambda a: [''.join(s).rstrip() for s in a]
strMatTrans  = lambda a: [''.join(s).rstrip() for s in zip(*a)]
sign = lambda x: math.copysign(1.0, x)

SODIR = os.path.join(os.getenv('HOME'),'.local','lib')

def print_table_latex(model, data):
	print('')
	print('Cumulative creep and fatigue damage and projected receiver life for a nitrate-salt Gemasolar-like receiver.')
	for i in range(int(model.nz/model.nbins)):
		lb = model.nbins*i
		ub = lb + model.nbins - 1
		j = np.argmin(data['max_cycles'][lb:ub])
		print('%d & %4.2f & %.1e & %.1e & %4.3f \\\\'%(
			i+1,
			model.H_rec*(j+1)/model.nbins,
			np.cumsum(data['Dc'],axis=0)[-1,lb+j],
			np.cumsum(data['Df'],axis=0)[-1,lb+j],
			data['max_cycles'][lb+j]/365,
		))

def id_cycles(times, period, days):
	"""
		Helper to separate out individual cycles by index

		Parameters:
			times       Tube times
			period      Period of a single cycle
			days        Number of days in simulation
	"""
	tm = np.mod(times, period)
	inds = list(np.where(tm == 0)[0])
	if len(inds) != (days + 1):
		raise ValueError("Tube times not compatible with the receiver number of days and cycle period!")
	return inds

def destring_array(string):
	"""
	Make an array from a space separated string
	"""
	return np.array(list(map(float, string.split(" "))))

def cycles_to_fail(material, pname, temp, erange):
	"""
		Returns fatigue cycles to failure at a given temperature and strain range

		Parameters:
		  pname:       property name ("nominalFatigue")
		  erange:      strain range in mm/mm
		  temp:        temperature in K
	"""
	pdata = material.data[pname]
	T, a, n, cutoff = [],[],[],[]

	for i in pdata:
		T.append(destring_array(pdata[i]["T"]))
		a.append(destring_array(pdata[i]["a"]))
		n.append(destring_array(pdata[i]["n"]))
		cutoff.append(destring_array(pdata[i]["cutoff"]))

	if np.array(a).shape != np.array(n).shape:
		raise ValueError("\tThe lists of a and n must have equal lengths!")

	inds=np.array(T).argsort(axis=0)
	T = np.array(T)[inds]
	a = np.array(a)[inds]
	n = np.array(n)[inds]
	cutoff = np.array(cutoff)[inds]

	if temp > max(T):
		print('	T: '+','.join('%s'%pdata[i]["T"] for i in pdata))
		print('	temp:  %s'%temp)
		raise ValueError("\ttemperature is out of range for cycle to failure determination")

	for i in range(np.size(T, axis=0)):
		if temp<=T[i]:
			polysum = 0.0
			if erange<=cutoff[i]:
				erange = cutoff[i][0][0]
			for (b,m) in zip(a[i][0],n[i][0]):
				polysum+=b*np.log10(erange)**m
			break

	return 10**polysum


def cycle_fatigue(strains, temperatures, material, nu = 0.5):
	"""
		Calculate fatigue damage for a single cycle

		Parameters:
			strains         single cycle strains
			temperatures    single cycle temperatures
			material        damage model

		Additional parameters:
			nu              effective Poisson's ratio to use
	"""
	pt_temps = np.max(temperatures, axis = 0)

	pt_eranges = np.zeros(pt_temps.shape)

	nt = strains.shape[1]
	for i in range(nt):
		for j in range(nt):
			de = strains[:,j] - strains[:,i]
			eq = np.sqrt(2) / (2*(1+nu)) * np.sqrt(
					(de[0] - de[1])**2 + (de[1]-de[2])**2 + (de[2]-de[0])**2.0
					+ 3.0/2.0 * (de[3]**2.0 + de[4]**2.0 + de[5]**2.0)
					)
			pt_eranges = np.maximum(pt_eranges, eq)

	dmg = np.zeros(pt_eranges.shape)
	for ind in np.ndindex(*dmg.shape):
		dmg[ind] = 1.0 / cycles_to_fail(material,"nominalFatigue", pt_temps[ind], pt_eranges[ind])

	return dmg


def calculate_max_cycles(Dc, Df, material, rep_min = 1, rep_max = 1e6):
	"""
		Actually calculate the maximum number of repetitions for a single point

		Parameters:
			Dc          creep damage per simulated cycle
			Df          fatigue damage per simulated cycle
			material    damaged material properties
	"""
	if not material.inside_envelope("cfinteraction", Df(rep_min), Dc(rep_min)):
		return 0

	if material.inside_envelope("cfinteraction", Df(rep_max), Dc(rep_max)):
		return np.inf

	return opt.brentq(lambda N: material.inside_envelope("cfinteraction", Df(N), Dc(N)) - 0.5,
			rep_min, rep_max)


def make_extrapolate(D, extrapolate="lump",order=1):
	"""
		Return a damage extrapolation function based on extrapolate
		giving the damage for the nth cycle

		Parameters:
			D:      raw, per cycle damage
	"""
	if extrapolate == "lump":
		return lambda N, D = D: N * np.sum(D) / len(D)
	elif extrapolate == "last":
		def Dfn(N, D = D):
			N = int(N)
			if N < len(D)-1:
				return np.sum(D[:N])
			else:
				return np.sum(D[:-1]) + D[-1] * N

		return Dfn
	elif extrapolate == "poly":
		p = np.polyfit(np.array(list(range(len(D))))+1, D, order)
		return lambda N, p=p: np.polyval(p, N)
	else:
		raise ValueError("Unknown damage extrapolation approach %s!" % extrapolate)

class Factors:
	def __init__(self,nodes=2,Ro=45./2000.):
		self.nt=nodes
		self.Ro=Ro
		self.dt=np.pi/(nodes-1)
		self.theta=np.linspace(0.0,np.pi,nodes)
	def alphas(self,theta):
		alpha = np.pi/2.
		alpha += -np.arcsin(np.cos(theta)/np.sqrt(5.-4*np.sin(theta)))
		alpha += -np.arcsin(2*np.sqrt(1-np.sin(theta))/np.sqrt(5.-4*np.sin(theta)))
		return alpha
	def string(self,thetai,thetaj):
		ai=self.alphas(thetaj)
		aj=self.alphas(thetai)
		ti=thetai
		tj=thetaj
		R=self.Ro
		if ai<=ti:
			#T1 is not used
			if aj<=tj:
				#T2 is not used
				L=R*np.sqrt(pow(2-np.sin(tj)-np.sin(ti),2) + pow(np.cos(tj)-np.cos(ti),2))
			else:
				#T2 is used
				L =R*np.sqrt(pow(2-np.sin(aj)-np.sin(ti),2) + pow(np.cos(aj)-np.cos(ti),2))
				L+=R*(aj-tj)
		else:
			#T1 is used
			if aj<=tj:
				#T2 is not used
				L =R*(ai-ti)
				L+=R*np.sqrt(pow(2-np.sin(tj)-np.sin(ai),2) + pow(np.cos(tj)-np.cos(ai),2))
			else:
				#T2 is used
				L =R*(ai-ti)
				L+=R*np.sqrt(pow(2-np.sin(aj)-np.sin(ai),2) + pow(np.cos(aj)-np.cos(ai),2))
				L+=R*(aj-tj)
		return L

	def factor(self,i,j):
		if self.theta[j]>=np.pi/2. and self.theta[i]<np.pi/2.:
			ti=max(self.theta[i]-self.dt/2.,0.)
			ti1=min(self.theta[i]+self.dt/2.,np.pi/2.)
			tj=max(self.theta[j]-np.pi/2.-self.dt/2.,0.)
			tj1=min(self.theta[j]-np.pi/2.+self.dt/2.,np.pi/2.)
			s1=self.string(ti,tj)
			d1=self.string(ti,tj1)
			d2=self.string(ti1,tj)
			s2=self.string(ti1,tj1)
			num=d1+d2-s1-s2
			den=2*self.Ro*(ti1-ti)
			F=num/den
		elif self.theta[j]<np.pi/2. and self.theta[i]>=np.pi/2.:
			ti=max(self.theta[i]-np.pi/2.-self.dt/2.,0.)
			ti1=min(self.theta[i]-np.pi/2.+self.dt/2.,np.pi/2.)
			tj=max(self.theta[j]-self.dt/2.,0.)
			tj1=min(self.theta[j]+self.dt/2.,np.pi/2.)
			s1=self.string(ti,tj)
			d1=self.string(ti,tj1)
			d2=self.string(ti1,tj)
			s2=self.string(ti1,tj1)
			num=d1+d2-s1-s2
			den=2*self.Ro*(ti1-ti)
			F=num/den
		else:
			F=0.0
		return F

	def aperture(self,i,inverse=False):
		ti=max(self.theta[i]-self.dt/2.,0.)
		ti1=min(self.theta[i]+self.dt/2.,np.pi/2.)
		if self.theta[i]>=np.pi/2:
			ti=max(self.theta[i]-np.pi/2.-self.dt/2.,0.)
			ti1=min(self.theta[i]-np.pi/2.+self.dt/2.,np.pi/2.)
		ai=self.alphas(ti)
		ai1=self.alphas(ti1)
		s1=self.Ro*ti
		d1=self.Ro*np.sqrt(pow(2-np.sin(ai)-np.sin(ti),2)+pow(np.cos(ai)-np.cos(ti),2)) + self.Ro*ai
		d2=self.Ro*ti1
		s2=self.Ro*np.sqrt(pow(2-np.sin(ai1)-np.sin(ti1),2)+pow(np.cos(ai1)-np.cos(ti1),2)) + self.Ro*ai1
		num=d1+d2-s1-s2
		den=2*self.Ro*(ti1-ti)
		if inverse:
			F=(ti1-ti)/2*num/den
		else:
			F=num/den
		return F

	def factors(self,debug=False):
		F=np.zeros((self.nt+1,self.nt+1))
		for i in range(self.nt+1):
			if i>0:
				F[i,0]=self.aperture(i-1)
		for j in range(self.nt+1):
			if j>0:
				F[0,j]=self.aperture(j-1,inverse=True)
		for i in range(1,self.nt+1):
			for j in range(1,self.nt+1):
				if not i==j:
					F[i,j]=self.factor(i-1,j-1)
		if debug:
			csv = np.sum(F,axis=1)
			F=np.c_[F,csv]
			np.savetxt("Factors.csv",F,delimiter=",")
		return F

class receiver_cyl:
	def __init__(self,coolant = 'salt', Ri = 57.93/2000, Ro = 60.33/2000, T_in = 290, T_out = 565,
                      nz = 450, nt = 46, nr = 18, R_fouling = 0.0, ab = 0.94, em = 0.88, kp = 16.57, H_rec = 10.5, D_rec = 8.5,
                      nbins = 50, alpha = 15.6e-6, Young = 186e9, poisson = 0.31, T0=0.,
                      thermat = "base",defomat = "const_base",damat = "base",mat='800H',
                      debugfolder = os.path.expanduser('~'), debug = False, verification = False, Dittus=True,maxiter=100):
		self.coolant = coolant
		self.Ri = Ri
		self.Ro = Ro
		self.thickness = Ro - Ri
		self.T_in = T_in + 273.15
		self.T_out = T_out + 273.15
		self.nz = nz
		self.nt = nt
		self.nr = nr
		self.R_fouling = R_fouling
		self.ab = ab
		self.em = em
		self.kp = kp
		self.H_rec = H_rec
		self.D_rec = D_rec
		self.dz = H_rec/nbins
		self.nbins=nbins
		self.debugfolder = debugfolder
		self.debug = debug
		self.verification = verification
		self.sigma = 5.670374419e-8
		# Discretisation parameters
		self.dt = np.pi/(nt-1)
		# Tube section diameter and area
		self.d = 2.*Ri                               # Tube inner diameter (m)
		self.area = 0.25 * np.pi * pow(self.d,2.)    # Tube flow area (m2)
		self.ln = np.log(Ro/Ri)                      # Log of Ro/Ri simplification
		#Auxiliary variables
		self.thetas=np.linspace(0.0,np.pi,self.nt)
		# View factors
		model=Factors(nodes=self.nt,Ro=self.Ro)
		self.F=model.factors()
		self.maxiter=maxiter

		# Vector Ho
		self.Ho=np.zeros(self.nt+1)
		for i in range(1,self.nt+1):
			if self.thetas[i-1]<np.pi/2.:
				ti=max(self.thetas[i-1]-self.dt/2.,0)
				tf=min(self.thetas[i-1]+self.dt/2.,np.pi/2)
			else:
				ti=max(self.thetas[i-1]-np.pi/2.-self.dt/2.,0)
				tf=min(self.thetas[i-1]-np.pi/2.+self.dt/2.,np.pi/2)
			self.Ho[i]=1/(tf-ti)*(np.sin(tf)-np.sin(ti))
		self.C_ab=self.matrixC(ab)
		self.C_em=self.matrixC(em)
		self.A_eb=self.matrixA()
		self.theta = np.linspace(-np.pi, np.pi,self.nt*2-2)
		self.n = 3
		self.l = alpha
		self.E = Young
		self.nu = poisson
		self.T0 = T0
		l = self.E*self.nu/((1+self.nu)*(1-2*self.nu));
		m = 0.5*self.E/(1+self.nu);
		props = l*np.ones((3,3)) + 2*m*np.identity(3)
		self.invprops = np.linalg.inv(props)
		# Loading dynamic library
		so_file = "%s/stress.so"%SODIR
		stress = ctypes.CDLL(so_file)
		self.fun = stress.curve_fit
		self.fun.argtypes = [ctypes.c_int,
						ctypes.c_double,
						ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
						ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
		self.Dittus = Dittus
		# Choose the material models
		self.thermal_mat, self.deformation_mat, self.damage_mat = library.load_material(mat, thermat, defomat, damat)
		pdata = self.damage_mat.data["nominalFatigue"]
		# Getting maximum temperature for fatigue damage regression
		T = [];
		for i in pdata:
			T.append(float(pdata[i]['T']))
		self.T_max_fatigue = max(T)

	def deltak(self,i,j):
		if i==j:
			return 1
		else:
			return 0

	def matrixC(self,ab):
		# Matrix C
		C=np.zeros((self.nt+1,self.nt+1))
		for i in range(self.nt+1):
			for j in range(self.nt+1):
				if j==0:
					C[i,j]=self.deltak(i,j)
				else:
					C[i,j]=self.deltak(i,j)/ab-(1/ab-1)*self.F[i,j]
		# Inverse
		Cinv=np.linalg.inv(C)
		return Cinv

	def matrixA(self):
		# Matrix A
		A=np.zeros((self.nt+1,self.nt+1))
		for i in range(self.nt+1):
			for j in range(self.nt+1):
				A[i,j]=self.deltak(i,j)-self.F[i,j]
		return A

	def import_mat(self,fileName):
		mat = scipy.io.loadmat(fileName, chars_as_strings=False)
		names = strMatTrans(mat['name']) # names
		descr = strMatTrans(mat['description']) # descriptions
		self.data = np.transpose(mat['data_2'])

		self._vars = {}
		self._blocks = []
		for i in range(len(names)):
			d = mat['dataInfo'][0][i] # data block
			x = mat['dataInfo'][1][i]
			c = abs(x)-1  # column
			s = sign(x)   # sign
			if c:
				self._vars[names[i]] = (descr[i], d, c, s)
				if not d in self._blocks:
					self._blocks.append(d)
				else:
					absc = (names[i], descr[i])
		del mat

	def density(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			d = 2090.0 - 0.636 * (T - 273.15)
		else:
			d = 219.0 + 275.32 * (1.0 - T / 2503.7) + 511.58 * np.sqrt(1.0 - T / 2503.7)
		return d

	def dynamicViscosity(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			eta = 0.001 * (22.714 - 0.120 * (T - 273.15) + 2.281e-4 * pow((T - 273.15),2) - 1.474e-7 * pow((T - 273.15),3))
			emax = 0.001 * (22.714 - 0.120 * 600. + 2.281e-4 * pow(600.,2) - 1.474e-7 * pow(600.,3))
			eta = np.maximum(eta,emax)
		else:
			eta = np.exp(-6.4406 - 0.3958 * np.log(T) + 556.835/T)
		return eta

	def thermalConductivity(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			k = 0.443 + 1.9e-4 * (T - 273.15)
		else:
			k = 124.67 - 0.11381 * T + 5.5226e-5 * pow(T,2) - 1.1842e-8 * pow(T,3);
		return k;

	def T_h(self,h):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			delta = np.abs(pow(1396.0182,2.) + 4*0.086*h)
			T = (-1396.0182 + np.sqrt(delta)) / (2 * 0.086)
		else:
			T = 1000 * (1.6582 - 8.4790e-4 * h + 4.4541e-7 * pow(h,2) - 2992.6 * pow(h,-2))
		return T

	def enthalpy(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			h = 1396.0182*T + 0.086 * pow(T,2.)
		else:
			h = 1000 * (1.6582 - 8.4790e-4 * T + 4.4541e-7 * pow(T,2) - 2992.6 * pow(T,-2))
		return h

	def specificHeatCapacityCp(self,T):
		"""
		   Thermal and transport properties of nitrate salt and sodium
		   Nitrate salt:  Zavoico, A. B. Solar power tower design basis document, revision 0; Topical. Sandia National Labs., 2001.
		   Liquid sodium: Fink, J. K.; Leibowitz, L. Thermodynamic and transport properties of sodium liquid and vapor. Argonne National Lab., 1995.
		"""
		if self.coolant == 'salt':
			C = 1396.0182 + 0.172 * T
		else:
			C = 1000 * (1.6582 - 8.4790e-4 * T + 4.4541e-7 * pow(T,2) - 2992.6 * pow(T,-2))
		return C

	def phis(self,x):
		return np.exp(-1./x)/(np.exp(-1./x) + np.exp(-1./(1-x)))

	def Temperature(self, m_flow, Tf, Tamb, CG, h_ext):
		"""
		    Flow and thermal variables:
		    hf: Heat transfer coefficient due to internal forced-convection
		    mu: HTF dynamic viscosity (Pa-s)
		    kf: HTF thermal conductivity (W/m-K)
		    C:  HTF specific heat capacity (J/kg-K)
		    Re: HTF Reynolds number
		    Pr: HTF Prandtl number
		    Nu: Nusselt number due to internal forced convection
		"""

		Tf,temp = np.meshgrid(np.ones(self.nt),Tf)
		Tf = Tf*temp

		# HTF thermo-physical properties
		mu = self.dynamicViscosity(Tf)                 # HTF dynamic viscosity (Pa-s)
		kf = self.thermalConductivity(Tf)              # HTF thermal conductivity (W/m-K)
		C = self.specificHeatCapacityCp(Tf)            # HTF specific heat capacity (J/kg-K)

		m_flow,temp = np.meshgrid(np.ones(self.nt), m_flow)
		m_flow = m_flow*temp

		Tamb,temp = np.meshgrid(np.ones(self.nt), Tamb)
		Tamb = Tamb*temp

		h_ext,temp = np.meshgrid(np.ones(self.nt), h_ext)
		h_ext = h_ext*temp

		# HTF internal flow variables
		Re = m_flow * self.d / (self.area * mu)    # HTF Reynolds number
		Pr = mu * C / kf                           # HTF Prandtl number
		Re_zer = np.where(Re<=0)[0]
		Re_pos = np.where(Re>0)[0]
		Re_lam = np.where((Re>0)     & (Re<=2e3))[0]
		Re_tra = np.where((Re>2e3) & (Re<4e3))[0]
		Re_tur = np.where(Re>=4e3)[0]
		f = np.zeros_like(Re)
		Nu= 4.36*np.ones_like(Re)
		Pr_neg = np.where(Pr<=0)[0]

		if self.coolant == 'salt':
			if self.Dittus:
				Nu = 0.023 * pow(Re, 0.8) * pow(Pr, 0.4)
			else:
				f[Re_lam] = 64/Re[Re_lam]
				f[Re_tur] = pow(1.82*np.log10(Re[Re_tur]) - 1.64, -2)
				f[Re_tra] = 64/Re[Re_tra] + self.phis((Re[Re_tra]-2e3)/2e3)*(pow(1.82*np.log10(Re[Re_tra]) - 1.64, -2.) - 64/Re[Re_tra])

				Nu[Re_lam] = 4.36
				Nu[Re_tur] = (f[Re_tur]/8)*(Re[Re_tur] - 1000)*Pr[Re_tur]/(1 + 12.7*pow(f[Re_tur]/8, 0.5)*(pow(Pr[Re_tur],0.66)-1))
				Nu[Re_tra] = 4.36 + self.phis((Re[Re_tra]-2e3)/2e3)*((f[Re_tra]/8.)*(Re[Re_tra] - 1000.)*Pr[Re_tra]/(1. + 12.7*pow(f[Re_tra]/8., 0.5)*(pow(Pr[Re_tra], 0.66) -1.)) - 4.36)

		else:
			Nu = 5.6 + 0.0165 * pow(Re*Pr, 0.85) * pow(Pr, 0.01)

		# HTF internal heat transfer coefficient
		hf = Nu * kf / self.d
		if self.R_fouling>0:
			hf[Re_pos] = 1./(1./hf[Re_pos] + self.R_fouling)

		# Calculating heat flux at circumferential nodes
		self.cosines = np.dot(self.C_ab,self.Ho)[1:]
		cosinesm,fluxes = np.meshgrid(self.cosines,CG)
		qabs = fluxes*cosinesm
		a = -((self.em*(self.kp + hf*self.ln*self.Ri)*self.Ro*self.sigma)/((self.kp + hf*self.ln*self.Ri)*self.Ro*(qabs + self.em*self.sigma*pow(Tamb,4)) + hf*self.kp*self.Ri*Tf + (self.kp + hf*self.ln*self.Ri)*self.Ro*Tamb*(h_ext)))
		b = -((hf*self.kp*self.Ri + (self.kp + hf*self.ln*self.Ri)*self.Ro*(h_ext))/((self.kp + hf*self.ln*self.Ri)*self.Ro*(qabs + self.em*self.sigma*pow(Tamb,4)) + hf*self.kp*self.Ri*Tf + (self.kp + hf*self.ln*self.Ri)*self.Ro*Tamb*(h_ext)))
		c1 = 9.*a*pow(b,2.) + np.sqrt(3.)*np.sqrt(-256.*pow(a,3.) + 27.*pow(a,2)*pow(b,4))
		c2 = (4.*pow(2./3.,1./3.))/pow(c1,1./3.) + pow(c1,1./3.)/(pow(2.,1./3.)*pow(3.,2./3.)*a)
		To = -0.5*np.sqrt(c2) + 0.5*np.sqrt((2.*b)/(a*np.sqrt(c2)) - c2)
		Ti = (To + hf*self.Ri*self.ln/self.kp*Tf)/(1 + hf*self.Ri*self.ln/self.kp)
		qnet = hf*(Ti - Tf)
		self.nt_half=int(self.nt/2)
		_qnet = qnet
		Qnet = _qnet.sum(axis=1)*self.Ri*self.dt*self.dz
		net_zero = np.where(Qnet<0)[0]
		Qnet[net_zero] = 0.0
		_qnet[net_zero,:] = 0.0
		self.qnet = _qnet
		self.qnet[:,self.nt_half:self.nt]=0.0
		self.hf = hf
		self.Q_abs = qabs.sum(axis=1)*self.Ro*self.dt*self.dz
		convergence = np.ones_like(To)
		it=0
		while (convergence>1e-4).any():
			# Radiation model
			Eb = self.sigma*pow(np.c_[Tamb[:,0],To],4)
			q_rad=np.dot(self.C_em,np.dot(self.A_eb,np.transpose(Eb)))[1:]
			q_rad=np.transpose(q_rad)
			self.Q_rad = q_rad.sum(axis=1)*self.Ro*self.dt*self.dz
			# Convection losses
			q_cnv = h_ext*(To - Tamb)
			self.Q_cnv = q_cnv.sum(axis=1)*self.Ro*self.dt*self.dz
			# Temperature update
			Ti = Tf + (qabs - q_rad - q_cnv)*self.Ro/self.Ri/hf
			To_new = Ti*(1 + hf*self.Ri*self.ln/self.kp) - hf*self.Ri*self.ln/self.kp*Tf
			convergence=np.abs(To-To_new)/To
			To = np.copy(To_new)
			it+=1
			if it>self.maxiter:
				raise(RuntimeError("The external temperature did not converge"))

		Qnet = self.Q_abs - self.Q_rad - self.Q_cnv

		for t in range(Ti.shape[0]):
			BDp = self.Fourier(Ti[t,:])

		# Fourier coefficients
		Ti[:,self.nt_half:self.nt]=Tf[:,self.nt_half:self.nt]
		To[:,self.nt_half:self.nt]=Tf[:,self.nt_half:self.nt]
#		self.stress, self.epsilon = self.crown_stress(Ti,To)
		self.Ti = Ti
		self.To = To
		return Qnet

	def Fourier(self,T):
		coefs = np.empty(21)
		self.fun(self.nt, self.dt, T, coefs)
		return coefs

	def crown_stress(self, Ti, To):
		stress = np.zeros((Ti.shape[0],self.nt,self.nr,6))
		strain = np.zeros((Ti.shape[0],self.nt,self.nr,6))
		ntimes = Ti.shape[0]
		for t in range(ntimes):
			thetas = np.linspace(0.0, np.pi, self.nt)
			BDp = self.Fourier(Ti[t,:])
			BDpp = self.Fourier(To[t,:])
			r = np.linspace(self.Ri,self.Ro,self.nr)
			for i in range(self.nt):
				T = np.linspace(Ti[t,i],To[t,i],self.nr)
				for j in range(self.nr):
					stress[t,i,j,:], strain[t,i,j,:] = self.Thermoelastic(T[j], r[j], thetas[i], BDp, BDpp)
		return stress,strain

	def Thermoelastic(self, T, r, theta, BDp, BDpp):
		Tbar_i = BDp[0]; BP = BDp[1]; DP = BDp[2];
		Tbar_o = BDpp[0]; BPP = BDpp[1]; DPP = BDpp[2];
		a = self.Ri; b = self.Ro; a2 = a*a; b2 = b*b; r2 = r*r; r4 = pow(r,4);

		C = self.l*self.E/(2.*(1. - self.nu));
		D = 1./(2.*(1. + self.nu));
		kappa = (Tbar_i - Tbar_o)/np.log(b/a);
		kappa_theta = r*a*b/(b2 - a2)*((BP*b - BPP*a)/(b2 + a2)*np.cos(theta) + (DP*b - DPP*a)/(b2 + a2)*np.sin(theta));
		kappa_tau   = r*a*b/(b2 - a2)*((BP*b - BPP*a)/(b2 + a2)*np.sin(theta) - (DP*b - DPP*a)/(b2 + a2)*np.cos(theta));

		T_theta = T - ((Tbar_i - Tbar_o) * np.log(b/r)/np.log(b/a)) - Tbar_o;

		Qr = kappa*C*(0 -np.log(b/r) -a2/(b2 - a2)*(1 -b2/r2)*np.log(b/a) ) \
					+ kappa_theta*C*(1 - a2/r2)*(1 - b2/r2);
		Qtheta = kappa*C*(1 -np.log(b/r) -a2/(b2 - a2)*(1 +b2/r2)*np.log(b/a) ) \
					+ kappa_theta*C*(3 -(a2 +b2)/r2 -a2*b2/r4);
		Qz = kappa*C*(1 -2*np.log(b/r) -2*a2/(b2 - a2)*np.log(b/a) ) \
					+ kappa_theta*2*self.nu*C*(2 -(a2 + b2)/r2) -self.l*self.E*T_theta;
		Qrtheta = kappa_tau*C*(1 -a2/r2)*(1 -b2/r2);

		Q_Eq = np.sqrt(0.5*(pow(Qr -Qtheta,2) + pow(Qr -Qz,2) + pow(Qz -Qtheta,2)) + 6*pow(Qrtheta,2));
		Q = np.zeros((6,))
		Q[0] = Qr; Q[1] = Qtheta; Q[2] = Qz; Q[3] = Qrtheta
		e = np.zeros((6,))
		e[0] = 1/self.E*(Qr - self.nu*(Qtheta + Qz))
		e[1] = 1/self.E*(Qtheta - self.nu*(Qr + Qz))
		e[2] = -self.l*T_theta

		if self.verification:
			print("=============== NPS Sch. 5S 1\" S31609 at 450degC ===============")
			print("Biharmonic coefficients:")
			print("Tbar_i [K]:      749.6892       %4.4f"%Tbar_i)
			print("  B'_1 [K]:      45.1191        %4.4f"%BP)
			print("  D'_1 [K]:      -0.0000        %4.4f"%DP)
			print("Tbar_o [K]:      769.7119       %4.4f"%Tbar_o)
			print(" B''_1 [K]:      79.4518        %4.4f"%BPP)
			print(" D''_1 [K]:      0.0000         %4.4f\n"%DPP)
			print("Stress at outside tube crown:")
			print("Q_r [MPa]:       0.0000         %4.4f"%(Qr/1e6))
			print("Q_rTheta [MPa]:  0.0000         %4.4f"%(Qrtheta/1e6))
			print("Q_theta [MPa]:  -101.0056       %4.4f"%(Qtheta/1e6))
			print("Q_z [MPa]:      -389.5197       %4.4f"%(Qz/1e6))
			print("Q_Eq [MPa]:      350.1201       %4.4f"%(Q_Eq/1e6))

		return Q,e

	def run_heuristics(self,days,step,rfile,case,casedir='.'):
		'''
		Run salt receiver heuristics and find each section damage
		   Inputs:
		      -days: Temporal bounds of the simulation (list of int [a,b])
		      -step: Simulation time step (float)
		      -rfile: SolarTherm output .mat file (str)
		   Output:
		      -
		'''

		self.import_mat(rfile)                                                                                # Importing SolarTherm output
		times = self.data[:,0]                                                                                # Get times
		CG = self.data[:,self._vars['heliostatField.CG[1]'][2]:self._vars['heliostatField.CG[450]'][2]+1]    # Get fluxes
		m_flow_tb = self.data[:,self._vars['heliostatField.m_flow_tb'][2]]                                    # Get mass flow rates
		Tamb = self.data[:,self._vars['receiver.Tamb'][2]]                                                    # Get ambient temperatures
		h_ext = self.data[:,self._vars['receiver.h_conv'][2]]                                                 # Get heat transfer coefficient due to external convection

		# Filtering times
		index = []
		_times = []
		time_lb = days[0]*86400
		time_ub = days[1]*86400
		print(yellow('	Sorting times'))
		for i in tqdm(range(len(times))):
			if times[i]%step==0 and times[i] not in _times and time_lb<=times[i] and times[i]<=time_ub:
				index.append(i)
				_times.append(times[i])

		times = times[index].flatten()/3600.                                                                  # Filter times
		CG = CG[index,:]                                                                                      # Filter fluxes
		m_flow_tb = m_flow_tb[index]                                                                          # Filter mass flow rates
		CG[np.where(m_flow_tb==0)[0],:] = 0.0                                                                 # Matching zero flow with zero flux
		Tamb = Tamb[index]                                                                                    # Filter ambient temperature
		h_ext = h_ext[index]                                                                                  # Filter heat transfer coefficient due to external convection

		# Instantiating variables
		field_off = [0]
		start = []
		stop = []
		for i in range(1,times.shape[0]-1):
			if m_flow_tb[i]==0 and m_flow_tb[i+1]==0 and m_flow_tb[i-1]==0:
				field_off.append(i)
			if m_flow_tb[i]==0 and m_flow_tb[i+1]>0 and m_flow_tb[i-1]==0:
				start.append(i)
			if m_flow_tb[i]==0 and m_flow_tb[i+1]==0 and m_flow_tb[i-1]>0:
				stop.append(i)

		field_off.append(times.shape[0]-1)
		sigma_o = np.zeros((times.shape[0],self.nz,6))
		sigma_i = np.zeros((times.shape[0],self.nz,6))
		epsilon_o = np.zeros((times.shape[0],self.nz,6))
		epsilon_i = np.zeros((times.shape[0],self.nz,6))
		T_o = np.zeros((times.shape[0],self.nz))
		T_i = np.zeros((times.shape[0],self.nz))
		Tf = self.T_in*np.ones((times.shape[0],self.nz+1))

		for i in field_off:
			Tf[i,:] = 293.15*np.ones((self.nz+1,))
		for i in start:
			Tf[i,:] = 533.15*np.ones((self.nz+1,))
		for i in stop:
			Tf[i,:] = 533.15*np.ones((self.nz+1,))
		qnet = np.zeros((times.shape[0],2*self.nt-1,self.nz))

		# Running thermal model
		print(yellow('	Running thermal model'))
		for k in tqdm(range(self.nz)):
			Qnet = self.Temperature(m_flow_tb, Tf[:,k], Tamb, CG[:,k], h_ext)
			C = self.specificHeatCapacityCp(Tf[:,k])*m_flow_tb
			Tf[:,k+1] = Tf[:,k] + np.divide(Qnet, C, out=np.zeros_like(C), where=C!=0)
			sigma_o[:,k,:] = self.stress[:,0,:]/1e6
			sigma_i[:,k,:] = self.stress[:,1,:]/1e6
			epsilon_o[:,k,:] = self.epsilon[:,0,:]
			epsilon_i[:,k,:] = self.epsilon[:,1,:]
			qnet[:,:,k] = self.qnet/1e6
			T_o[:,k] = self.To
			T_i[:,k] = self.Ti

		sigmaEq_o = np.sqrt((
			(sigma_o[:,:,0] - sigma_o[:,:,1])**2.0 + 
			(sigma_o[:,:,1] - sigma_o[:,:,2])**2.0 + 
			(sigma_o[:,:,2] - sigma_o[:,:,0])**2.0 + 
			6.0 * (sigma_o[:,:,3]**2.0 + sigma_o[:,:,4]**2.0 + sigma_o[:,:,5]**2.0))/2.0)

		sigmaEq_i = np.sqrt((
			(sigma_i[:,:,0] - sigma_i[:,:,1])**2.0 + 
			(sigma_i[:,:,1] - sigma_i[:,:,2])**2.0 + 
			(sigma_i[:,:,2] - sigma_i[:,:,0])**2.0 + 
			6.0 * (sigma_i[:,:,3]**2.0 + sigma_i[:,:,4]**2.0 + sigma_i[:,:,5]**2.0))/2.0)

		inv = np.where(T_o>self.T_max_fatigue)[0]
		temp,T_o[inv,:] = np.meshgrid(np.ones((1,450)),Tf[inv,0].flatten())
		temp,T_i[inv,:] = np.meshgrid(np.ones((1,450)),Tf[inv,0].flatten())
		temp,Tf[inv,1:] = np.meshgrid(np.ones((1,450)),Tf[inv,0].flatten())
		sigmaEq_o[inv,:] = 0.0
		sigmaEq_i[inv,:] = 0.0
		epsilon_o[inv,:,:] = 0.0
		epsilon_i[inv,:,:] = 0.0

		# Time to rupture
		period = 24
		tR = self.damage_mat.time_to_rupture("averageRupture", T_o, sigmaEq_o)
		dts = np.diff(times)
		time_dmg = dts[:,np.newaxis]/tR[1:]

		# Break out to cycle damage
		inds = id_cycles(times, period, days[1]-days[0])

		# Saving
		data = {}
		data['sigma_o'] = sigma_o
		data['sigma_i'] = sigma_i
		data['epsilon_o'] = epsilon_o
		data['epsilon_i'] = epsilon_i
		data['T_o'] = T_o
		data['T_i'] = T_i
		data['times'] = times
		data['Tf'] = Tf
		data['m_flow_tb'] = m_flow_tb
		data['CG'] = CG
		data['nbins'] = self.nbins
		scipy.io.savemat(os.path.join(casedir,'%s.mat'%case),data)

		# Cycle damage
		Dc = np.array([np.sum(time_dmg[inds[i]:inds[i+1]], axis = 0) for i in range(days[1]-days[0])])

		### Fatigue cycles ###
		temperatures = T_o
		strains = epsilon_o

		# Run through each cycle and ID max strain range and fatigue damage
		strain_names = [0, 1, 2, 3, 4, 5] #[ex,ey,ez,exy,exz,eyz]
		strain_factors = [1.0,1.0,1.0,2.0, 2.0, 2.0]

		try:
			Df =  np.array([cycle_fatigue(np.array([ef*strains[inds[i]:inds[i+1],:,en] for en,ef in zip(strain_names, strain_factors)]),
			                          temperatures[inds[i]:inds[i+1]], 
			                          self.damage_mat) for i in range(days[1]-days[0])])
		except:
			Df = 1e-6*np.ones_like(Dc)
			array = np.c_[times,m_flow_tb,T_o]
			header = 'times,m_flow' + ''.join(',z[%d]'%(i+1) for i in range(self.nz))
			np.savetxt(os.path.join(casedir,'debug%s.csv'%case),array,fmt='%s',delimiter=',',header=header)

		### Calculating the number of cycles

		# Defining the number of columns as the number of days
		# This is used to create an array with nrows = nelements x nquad,
		# and ncols = number of days
		nc = days[1] - days[0]
		max_cycles = []

		for c,f in zip(Dc.reshape(nc,-1).T, Df.reshape(nc,-1).T):
			# The damage is extrapolated and the number of cycles is determined
			# There are three extrapolation approaches. Here we use the 'lump' one
			max_cycles.append(calculate_max_cycles(make_extrapolate(c), make_extrapolate(f), self.damage_mat))

		max_cycles = np.array(max_cycles)

		# Saving
		data['Dc'] = Dc
		data['Df'] = Df
		data['max_cycles'] = max_cycles
		data['min_cycle'] = np.argmin(max_cycles)
		data['max_creep'] = np.argmin(np.cumsum(Dc, axis=0))
		data['max_fatig'] = np.argmin(np.cumsum(Df, axis=0))
		print_table_latex(self,data)
		scipy.io.savemat(os.path.join(casedir,'%s.mat'%case),data)

		if (days[1] - days[0])==1:
			# Creating subplots
			fig, axes = plt.subplots(2,4, figsize=(18,8))
			# Tube front stress (inner)
			axes[0,0].plot(times,sigmaEq_i)
			axes[0,0].set_xlabel(r'$t$ [h]')
			axes[0,0].set_ylabel(r'$\sigma_\mathrm{max,i}$ [MPa]')
			# Tube front stress (outer)
			axes[0,1].plot(times,sigmaEq_o)
			axes[0,1].set_xlabel(r'$t$ [h]')
			axes[0,1].set_ylabel(r'$\sigma_\mathrm{max,o}$ [MPa]')
			# Tube front temperatures
			axes[0,2].plot(times,T_i-273.15)
			axes[0,2].set_xlabel(r'$t$ [h]')
			axes[0,2].set_ylabel(r'$T_\mathrm{i}$ [\textdegree C]')
			# Tube front temperatures
			axes[0,3].plot(times,T_o-273.15)
			axes[0,3].set_xlabel(r'$t$ [h]')
			axes[0,3].set_ylabel(r'$T_\mathrm{o}$ [\textdegree C]')
			# Tube front temperature vs z
			z = np.linspace(0,self.H_rec*self.nz/self.nbins,self.nz)
			axes[1,0].plot(z,np.transpose(T_i)-273.15)
			axes[1,0].set_xlabel(r'$z$ [m]')
			axes[1,0].set_ylabel(r'$T_\mathrm{front}$ [\textdegree C]')
			# Tube front temperature vs z
			axes[1,1].plot(z,np.transpose(T_o)-273.15)
			axes[1,1].set_xlabel(r'$z$ [m]')
			axes[1,1].set_ylabel(r'$T_\mathrm{front}$ [\textdegree C]')
			# Fluid temperature vs z
			axes[1,2].plot(z,np.transpose(Tf)[1:,:]-273.15)
			axes[1,2].set_xlabel(r'$z$ [m]')
			axes[1,2].set_ylabel(r'$T_\mathrm{fluid}$ [\textdegree C]')
			# Fluid temperature vs z
			axes[1,3].plot(times,Tf[:,-1]-273.15)
			axes[1,3].set_xlabel(r'$t$ [h]')
			axes[1,3].set_ylabel(r'$T_\mathrm{fluid}$ [\textdegree C]')
			plt.tight_layout()
			plt.savefig(os.path.join(casedir,'%s.png'%case),dpi=300)

	def thermal_verification(self,filename):

		# Getting the data

		fileo = open(filename,'rb')
		data = pickle.load(fileo)
		fileo.close()

		num_fp = len(data['fp'])
		z = np.linspace(0,1,self.nz)

		Ti = []
		To = []
		si = []
		so = []

		for fp in range(num_fp):
			# Instantiating receiver model
			h_ext = data['h_conv_ext']
			Tamb = data['T_amb']
			CG = flux_in = data['flux_in'][fp]
			m_flow_tb = data['m'][fp]/data['n_tubes'][fp]

			# Instantiating variables
			Tf = self.T_in*np.ones(self.nz+1)
			T_o = np.zeros(self.nz)
			T_i = np.zeros(self.nz)
			s_o = np.zeros(self.nz)
			s_i = np.zeros(self.nz)

			# Running thermal model
			for k in range(self.nz):
				Qnet = self.Temperature(m_flow_tb, Tf[k], Tamb, CG[k], h_ext)
				C = self.specificHeatCapacityCp(Tf[k])*m_flow_tb
				Tf[k+1] = Tf[k] + Qnet/C
				T_o[k] = self.To
				T_i[k] = self.Ti
				sigmaEq_o = np.sqrt((
				(self.stress[:,:,0] - self.stress[:,:,1])**2.0 + 
				(self.stress[:,:,1] - self.stress[:,:,2])**2.0 + 
				(self.stress[:,:,2] - self.stress[:,:,0])**2.0 + 
				6.0 * (self.stress[:,:,3]**2.0 + self.stress[:,:,4]**2.0 + self.stress[:,:,5]**2.0))/2.0)

				s_o[k] = sigmaEq_o[0,0]/1e6
				s_i[k] = sigmaEq_o[0,1]/1e6

			Ti.append(T_i)
			To.append(T_o)
			si.append(s_i)
			so.append(s_o)
		return Ti,To,si,so

class heuristics:
	def __init__(self, res='heuristics_res.mat', mat='800H', thermat = 'base', defomat='const_base', damat='base', folder='.'):
		self.folder = folder
		self.mat = mat
		self.defomat = defomat                                    # base | elastic_creep | elastic_model | const_elastic_creep | const_base
		self.thermal_mat, self.deformation_mat, self.damage_mat = library.load_material(mat, thermat, defomat, damat)
		self.res = res
		data = scipy.io.loadmat('%s/%s'%(self.folder,self.res))
		self.times = data['times'].flatten()
		self.nbins = data['nbins'].flatten()[0]
		self.epsilon_o = data['epsilon_o']
		self.sigma_o = data['sigma_o']
		self.T_o = data['T_o']
		self.strain = np.sqrt(2.0)/3.0 * np.sqrt(
			  (self.epsilon_o[:,:,0] - self.epsilon_o[:,:,1])**2.0
			+ (self.epsilon_o[:,:,1] - self.epsilon_o[:,:,2])**2.0
			+ (self.epsilon_o[:,:,2] - self.epsilon_o[:,:,0])**2.0
			+ 6.0 * (self.epsilon_o[:,:,3]**2.0
			+ self.epsilon_o[:,:,4]**2.0
			+ self.epsilon_o[:,:,5]**2.0)
			)
		self.stress = np.sqrt((
			(self.sigma_o[:,:,0] - self.sigma_o[:,:,1])**2.0 + 
			(self.sigma_o[:,:,1] - self.sigma_o[:,:,2])**2.0 + 
			(self.sigma_o[:,:,2] - self.sigma_o[:,:,0])**2.0 + 
			6.0 * (self.sigma_o[:,:,3]**2.0 + 
			self.sigma_o[:,:,4]**2.0 + 
			self.sigma_o[:,:,5]**2.0))/2.0)
		del data

	def uniaxial_neml(self,panel,position):
		lb = int(self.nbins*(panel-1) + position)
		stress_corr = self.unneml(lb)
		filename = '%s/%s'%(self.folder,self.res)
		# Calculating damage
		tR,time_dmg,Dc,Df,max_cycles = self.creepdmg(self.T_o[:,lb],stress_corr)
		print('%s,%d,%d,%.6f,%.6f,%.2f'%(self.res,panel,int(position),np.cumsum(Dc)[-1],np.cumsum(Df)[-1],max_cycles[0]))
		# Saving data
		csv = np.c_[self.times,self.T_o[:,lb],self.stress[:,lb],stress_corr,self.epsilon_o[:,lb,:]]
		headers = ['times','temp_o','stress_eq_o','stress_eq_corr','strain_xx','strain_yy','strain_zz','strain_xy','strain_xz','strain_yz']
		header = ''
		for s in headers:
			header += '%s_neml_%d,'%(s,panel)
		np.savetxt('neml_panel_%d_position_%s.csv'%(panel,int(position)),csv,delimiter=',',header=header)

	def creepdmg(self,T,stress_corr):
		period = 24
		tR = self.damage_mat.time_to_rupture("averageRupture", T, stress_corr)
		dts = np.diff(self.times)
		time_dmg = dts/tR[1:]

		# Break out to cycle damage
		nt = int((self.times[-1]-self.times[0])/24)
		inds = id_cycles(self.times, period, nt)

		# Cycle damage
		Dc = np.array([np.sum(time_dmg[inds[i]:inds[i+1]], axis = 0) for i in range(nt)])

		temperatures = self.T_o
		strains = self.epsilon_o

		# Run through each cycle and ID max strain range and fatigue damage
		strain_names = [0, 1, 2, 3, 4, 5] #[ex,ey,ez,exy,exz,eyz]
		strain_factors = [1.0,1.0,1.0,2.0, 2.0, 2.0]

		Df =  np.array([cycle_fatigue(np.array([ef*strains[inds[i]:inds[i+1],:,en] for en,ef in zip(strain_names, strain_factors)]),
			temperatures[inds[i]:inds[i+1]], 
			self.damage_mat) for i in range(nt)])

		### Calculating the number of cycles

		nc = nt
		max_cycles = []

		for c,f in zip(Dc.reshape(nc,-1).T, Df.reshape(nc,-1).T):
			max_cycles.append(calculate_max_cycles(make_extrapolate(c), make_extrapolate(f), self.damage_mat))

		max_cycles = np.array(max_cycles)

		return tR,time_dmg,Dc,Df,max_cycles

	def tablecsv(self,filename,stress_corr,tR,time_dmg,lb):
		e = self.epsilon_o[:,lb,:]
		f = open('%s.csv'%(filename),'+w')
		f.write('time_neml,Temp_neml,stress_nts,stress_nml,e_xx_neml,e_yy_neml,e_zz_neml,e_xy_neml,e_xz_neml,e_yz_neml,tR,time_dmg\n')
		for t,T,sn,scorr,e1,e2,e3,e4,e5,e6,tr,td in zip(self.times,self.T_o[:,lb],self.stress[:,lb],stress_corr,e[:,0],e[:,1],e[:,2],e[:,3],e[:,4],e[:,5],tR,time_dmg):
			f.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(t,T,sn,scorr,e1,e2,e3,e4,e5,e6,tr,td))
		f.close()

	def unneml(self, lb):
		# Choose the material models
		deformation_mat = library.load_deformation(self.mat, self.defomat)
		neml_model = deformation_mat.get_neml_model()
		umodel = uniaxial.UniaxialModel(neml_model, verbose = False)

		hn = umodel.init_store()
		en = self.strain[0,lb]
		sn = self.stress[0,lb]
		Tn = self.T_o[0,lb]
		tn = self.times[0]
		un = 0.0
		pn = 0.0;
		es = [en]
		ss = [sn]
		ep=0

		for enp1, Tnp1, tnp1, snt in tqdm(zip(self.strain[1:,lb],self.T_o[1:,lb],self.times[1:],self.stress[1:,lb]),total=len(self.times[1:])):
			snp1, hnp1, Anp1, unp1, pnp1 = umodel.update(enp1, en, Tnp1, Tn, tnp1, tn, sn, hn, un, pn)
			ss.append(abs(snp1))
			sn = snp1
			hn = np.copy(hnp1)
			en = enp1
			Tn = Tnp1
			tn = tnp1
			un = unp1
			pn = pnp1
		return np.array(ss)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Estimates average damage of a representative tube in a receiver panel')
	# Heuristics
	parser.add_argument('--days', nargs=2, type=int, default=[0,1])
	parser.add_argument('--step', type=float, default=1800)
	parser.add_argument('--heuristics', type=int, default=1)
	parser.add_argument('--file', type=str, default=os.path.join(os.getenv('HOME'),'solartherm/examples/GemasolarSystemOperation_res.mat'))
	parser.add_argument('--case', type=str, default='heuristics_res')
	parser.add_argument('--mat', type=str, default='800H')
	parser.add_argument('--DO', type=float, default=22.4)
	parser.add_argument('--WT', type=float, default=1.20)
	parser.add_argument('--casedir', type=str, default=os.path.join(os.getcwd(),'HEURISTICS-SALT'))
	# UniaxialModel
	parser.add_argument('--neml', type=int, default=1)
	parser.add_argument('--neml_res', type=str, default='heuristics_res.mat')
	parser.add_argument('--panel', type=float, default=4)
	parser.add_argument('--position', type=float, default=30)
	args = parser.parse_args()

	tinit = time.time()

	# Material
	materials = {
		'A230':{'kp':17.97,'alpha':15.61e-06,'Young':184e9,'poisson':0.31},
		'800H':{'kp':18.30,'alpha':18.28e-06,'Young':171e9,'poisson':0.31}
	}

	# Heuristics
	model = receiver_cyl(
		Ri = args.DO/2000.0 - args.WT/1000.0
		,Ro = args.DO/2000.0
		,R_fouling=8.808e-5
		,alpha = materials[args.mat]['alpha']
		,Young = materials[args.mat]['Young']
		,poisson = materials[args.mat]['poisson']
		,ab = 0.93
		,em = 0.87
		,kp = materials[args.mat]['kp']
		,Dittus=False
		,mat=args.mat)

	args.case=args.file.split('/')[-1].split('.mat')[0]
	print('	Ri: %.2f'%(model.Ri*1000),
		'	Ro: %.2f'%(model.Ro*1000),
		'	alpha: %g'%model.l,
		'	Young: %.1f'%(model.E/1e9),
		'	poisson: %.2f'%model.nu,
		'	kp: %s'%model.kp)
	print('	SolarTherm res: %s'%args.file)
	print('	Material: %s'%args.mat)
	print('	Case: %s'%args.case)
	print('	Casedir: %s'%args.casedir)
	model.run_heuristics(days=args.days,step=args.step,rfile=args.file,case=args.case,casedir=args.casedir)

	seconds = time.time() - tinit
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
