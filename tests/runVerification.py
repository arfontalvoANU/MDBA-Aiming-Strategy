#!/usr/bin/env python3

import os
from nitrateSaltHeuristics import *
import numpy as np
from mdbapy.cal_sun import SunPosition

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 14

if __name__=='__main__':

	mat = '800H'

	# Material
	materials = {
		'A230':{'kp':17.97,'alpha':15.61e-06,'Young':184e9,'poisson':0.31},
		'800H':{'kp':18.30,'alpha':18.28e-06,'Young':171e9,'poisson':0.31}
	}

	N=10  # for lammda, ecliptic longitude
	M=24  # for omega

	Lammda=np.linspace(-np.pi,np.pi,N+1)
	Omega=np.linspace(-np.pi,np.pi,M+1)
	sun=SunPosition()
	dni = [1.00,1.39,0.87,0.56]
	latitude = 37.56

	parent=os.path.join(os.getenv('HOME'),'ownCloud/phd_update/damage')
	for material in ['N08810']:
		for case,DO,WT in zip(['N08810_OD22.40_WT1.20_600','N08810_OD45.00_WT1.50_600'],[22.4,45.],[1.2,1.5]):
			casedir=os.path.join(parent,material,case)

			# Heuristics
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

			z = np.linspace(0,1,model.nz)

			for r,d in enumerate(dni):
				print('	case: %s 	d: %s'%(case,d))
				with PdfPages('%s/flux_table_d%s.pdf'%(casedir,d)) as pdf:
					for n in range(3,8):
						for m in range(int(0.5*M)+1):
							delta = 23.4556*np.sin(Lammda[n])
							theta=sun.zenith(latitude, delta, Omega[m]/np.pi*180.)
							phi=sun.azimuth(latitude, theta, delta, Omega[m]/np.pi*180.)
							elevation=90.-theta
							if elevation<=8.:
								continue
							try:
								DNI=d*1363.*pow(0.7,pow(1./np.cos(theta/180.*np.pi),0.678))
								Ti,To,si,so = model.thermal_verification('%s/flux_table_n%d_m%d_d%s'%(casedir,n,m,d))

								fig = plt.figure(figsize=(12, 8))
								fig.suptitle(r'$G_b=%.2f$ $\mathrm{W}\cdot\mathrm{m}^{-2}$, $\gamma=%.2f$\textdegree, $\alpha=%.2f$\textdegree'%(DNI,phi,elevation))

								axes1 = fig.add_subplot(2,2,1)
								axes1.plot(z, Ti[0]-273.15, label=r'$T_i$')
								axes1.plot(z, To[0]-273.15, label=r'$T_o$')
								axes1.set_xlabel(r'$z$ [m]')
								axes1.set_ylabel(r'$T$ [\textdegree C]')
								axes1.legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

								axes2 = fig.add_subplot(2,2,2)
								axes2.plot(z, si[0], label=r'$\sigma_{i}$')
								axes2.plot(z, so[0], label=r'$\sigma_{o}$')
								axes2.set_xlabel(r'$z$ [m]')
								axes2.set_ylabel(r'$\sigma_\mathrm{eq}$ [MPa]')
								axes2.legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

								axes3 = fig.add_subplot(2,2,3)
								axes3.plot(z, Ti[1]-273.15, label=r'$T_i$')
								axes3.plot(z, To[1]-273.15, label=r'$T_o$')
								axes3.set_xlabel(r'$z$ [m]')
								axes3.set_ylabel(r'$T$ [\textdegree C]')
								axes3.legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

								axes4 = fig.add_subplot(2,2,4)
								axes4.plot(z, si[1], label=r'$\sigma_{i}$')
								axes4.plot(z, so[1], label=r'$\sigma_{o}$')
								axes4.set_xlabel(r'$z$ [m]')
								axes4.set_ylabel(r'$\sigma_\mathrm{eq}$ [MPa]')
								axes4.legend(loc="best", borderaxespad=0, ncol=1, frameon=False)

								fig.tight_layout()
#								fig.subplots_adjust(top=0.9)
								pdf.savefig()  # saves the current figure into a pdf page
								plt.close('all')

							except FileNotFoundError:
								continue
