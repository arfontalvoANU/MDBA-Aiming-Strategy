#!/usr/bin/env python3

import os, pickle, colorama
import numpy as np
from tqdm import tqdm
from mdbapy.cal_sun import SunPosition
from scipy.interpolate import interp2d

colorama.init()
def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

class PostProcess:
	def __init__(self,dni,num_fp,num_el,casedir='.',latitude=37.56):
		self.latitude = latitude
		self.casedir = casedir
		self.dni=dni
		self.num_fp=num_fp
		self.num_el=num_el

	def get_fluxes(self): # OELT and RELT generation

		print(yellow('	Getting fluxes'))

		N=10  # for lammda, ecliptic longitude
		M=24  # for omega
		F=np.arange((N+1)*(M+1)*7,dtype=float).reshape(N+1,M+1,7) # the OELT, 0-field eff, 1-unavail,2-cosine, 3-reflective,4-sb,5-att,6-spi

		Lammda=np.linspace(-np.pi,np.pi,N+1)
		Omega=np.linspace(-np.pi,np.pi,M+1)
		sun=SunPosition()

		for r,d in enumerate(self.dni):
			for fp in range(self.num_fp):
				print('	d,fp: ({0:.2f},{1:g})'.format(d,fp))
				f = open('%s/FLUX_fp%s_d%s.motab'%(self.casedir,fp+1,d),'w+')
				f.write('#1\n')
				for k in tqdm(range(self.num_el)):
					F[:,:,:]=0.
					for n in range(3,8):
						for m in range(int(0.5*M)+1):
							delta = 23.4556*np.sin(Lammda[n])
							theta=sun.zenith(self.latitude, delta, Omega[m]/np.pi*180.)
							phi=sun.azimuth(self.latitude, theta, delta, Omega[m]/np.pi*180.)
							elevation=90.-theta
							if elevation<=8.:
								continue

							try:
								fileo = open('%s/flux_table_n%d_m%d_d%s'%(self.casedir,n,m,d),'rb')
								data = pickle.load(fileo)
								fileo.close()
								F[n,m,:] = data['flux_in'][fp][k]
							except:
								continue

						for m in range(int(0.5*M)+1,M+1):
							F[n,m,:]=F[n,M-m,:]

					for n in range(3):
						F[n,:,:]=F[5-n,:,:]

					for n in range(8,11):
						F[n,:,:]=F[15-n,:,:]

					F_output=np.arange((N+2)*(M+2),dtype=float).reshape(N+2,M+2)
					F_output[0,1:]=Omega/np.pi*180.
					F_output[1:,0]=Lammda/np.pi*180.
					F_output[1:,1:]=F[:,:,0]

					f.write('double flux_%s(%d,%d)\n'%(k+1, N+2, M+2))
					for i in range(F_output.shape[0]):
						for j in range(F_output.shape[1]):
							if i==0:
								f.write('%.1f '%F_output[i,j])
							elif j==0:
								f.write('%.1f '%F_output[i,j])
							else:
								f.write('%.1f '%F_output[i,j])
						f.write('\n')
					f.write('\n')
				f.close()

	def get_mflow(self): # OELT and RELT generations

		print(yellow('	Getting mass flow rates'))

		N=10  # for lammda, ecliptic longitude
		M=24  # for omega
		F=np.arange((N+1)*(M+1)*7,dtype=float).reshape(N+1,M+1,7) # the OELT, 0-field eff, 1-unavail,2-cosine, 3-reflective,4-sb,5-att,6-spi

		Lammda=np.linspace(-np.pi,np.pi,N+1)
		Omega=np.linspace(-np.pi,np.pi,M+1)
		sun=SunPosition()

		for fp in range(self.num_fp):
			f = open('%s/MFLOW_Solstice_fp%s.motab'%(self.casedir,fp+1),'w+')
			f.write('#1\n')
			for r,d in enumerate(self.dni):
				print('	d,fp: ({0:.2f},{1:g})'.format(d,fp))
				F[:,:,:]=0.
				for n in range(3,8):
					for m in range(int(0.5*M)+1):
						delta = 23.4556*np.sin(Lammda[n])
						theta=sun.zenith(self.latitude, delta, Omega[m]/np.pi*180.)
						phi=sun.azimuth(self.latitude, theta, delta, Omega[m]/np.pi*180.)
						elevation=90.-theta
						if elevation<=8.:
							continue
						try:
							fileo = open('%s/flux_table_n%d_m%d_d%s'%(self.casedir,n,m,d),'rb')
							data = pickle.load(fileo)
							fileo.close()
							F[n,m,:] =  data['m'][fp]/data['n_tubes'][fp]
						except:
							continue

					for m in range(int(0.5*M)+1,M+1):
						F[n,m,:]=F[n,M-m,:]

				for n in range(3):
					F[n,:,:]=F[5-n,:,:]

				for n in range(8,11):
					F[n,:,:]=F[15-n,:,:]

				F_output=np.arange((N+2)*(M+2),dtype=float).reshape(N+2,M+2)
				F_output[0,1:]=Omega/np.pi*180.
				F_output[1:,0]=Lammda/np.pi*180.
				F_output[1:,1:]=F[:,:,0]

				f.write('# DNI ratio: %s\n'%(d))
				f.write('double mflow_%s(%d,%d)\n'%(r+1, N+2, M+2))
				for i in range(F_output.shape[0]):
					for j in range(F_output.shape[1]):
						if i==0:
							f.write('%.1f '%F_output[i,j])
						elif j==0:
							f.write('%.1f '%F_output[i,j])
						else:
							f.write('%s '%F_output[i,j])
					f.write('\n')
				f.write('\n')
			f.close()

	def Tables(self): # OELT and RELT generations
		motab = open('%s/OELT_Solstice.motab'%self.casedir,'r')
		p = motab.readlines()
		motab.close()

		f = open('%s/OELTs_Solstice.motab'%(self.casedir),'w+')
		g = open('%s/HALTs_Solstice.motab'%(self.casedir),'w+')
		for k in range(5):
			f.write('%s'%p[k])
			g.write('%s'%p[k])

		for r,d in enumerate(self.dni):
			E = np.genfromtxt('%s/F_optic_%s.csv'%(self.casedir,d),delimiter=',')
			F = np.genfromtxt('%s/F_unavail_%s.csv'%(self.casedir,d),delimiter=',')

			f.write('double optics_%d(%d,%d)\n'%(r,E.shape[0],E.shape[1]))
			g.write('double optics_%d(%d,%d)\n'%(r,E.shape[0],E.shape[1]))
			for i in range(E.shape[0]):
				for j in range(E.shape[1]):
					if i==0:
						f.write('%.1f '%E[i,j])
						g.write('%.1f '%F[i,j])
					elif j==0:
						f.write('%.1f '%E[i,j])
						g.write('%.1f '%F[i,j])
					else:
						f.write('%s '%E[i,j])
						g.write('%s '%F[i,j])
				f.write('\n')
				g.write('\n')
			f.write('\n')
			g.write('\n')

		f.close()
		g.close()

