#!/usr/bin/env python3

import os
import numpy as np
import argparse
import pickle
from mdbapy.python_postprocessing import proces_raw_results
from mdbapy.one_key_co_optimisation import one_key_start
from mdbapy.Flux_reader import read_data
from mdbapy.Loss_analysis import receiver_correlation
from tqdm import tqdm

import colorama
colorama.init()
def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

class postprocess:
	def __init__(self, T=565, D0=22.4, WT=1.2, latitude=37.56, material='N06230'):
		self.dnir = [0.56,0.87,1.0,1.2,1.39]
		self.T = T
		self.D0 = D0
		self.WT = WT
		self.latitude=latitude
		self.csv_trimmed = os.path.join(os.getcwd(),'dnir0.56','job0','pos_and_aiming_trimmed.csv')
		self.num_hst=int(len(np.loadtxt(self.csv_trimmed,delimiter=',', skiprows=2))) # the num hst for the large field
		self.hst_w=10.5623289
		self.hst_h=10.5623289
		cases = np.genfromtxt('cases.csv',delimiter=',')
		self.dni = cases[:,5]
		self.jobs= cases[:,0].astype(int)
		self.n= cases[:,3].astype(int)
		self.m= cases[:,4].astype(int)
		self.r_diameter,self.r_height=8.5,10.5
		self.tower_h=114.75
		self.SM=2.4

	def get_fluxes(self): # OELT and RELT generations
		N=10
		M=24
		for fp in range(2):
			g = open('MFLOW_Solstice_fp{0}.motab'.format(fp+1),'w+')
			g.write('#1\n')
			for dr,d in enumerate(self.dnir):
				print('	d,fp: ({0:.2f},{1:g})'.format(d,fp))
				f = open('FLUX_fp{0}_d{1}.motab'.format(fp+1,d),'w+')
				f.write('#1\n')
				for k in tqdm(range(450)):
					F=np.zeros((N+1,M+1))
					G=np.zeros((N+1,M+1))
					for i,n,m in zip(self.jobs,self.n,self.m):
						casedir = os.path.join(os.getcwd(),'dnir{0}'.format(d),'job{0}'.format(i))
						fileo = open('%s/flux-table'%casedir,'rb')
						data = pickle.load(fileo)
						fileo.close()
						if abs(data['T_HC'][0][449]-data['T_out'])<1.0 and abs(data['T_HC'][1][449]-data['T_out'])<1.0:
							# Flux map
							F[n,m] = data['flux_in'][fp][k]
							# Mass flow rate
							G[n,m] = data['m'][fp]/data['n_tubes'][fp]

					for n in range(3):
						F[n,:]=F[5-n,:]
						G[n,:]=G[5-n,:]

					for n in range(8,11):
						F[n,:]=F[15-n,:]
						G[n,:]=G[15-n,:]

#					F[:,13:M+1]=F[:,0:12][:,::-1]
#					G[:,13:M+1]=G[:,0:12][:,::-1]

					F_output=np.zeros((N+2,M+2))
					F_output[0,1:]=np.linspace(-180,180,M+1)
					F_output[1:,0]=np.linspace(-180,180,N+1)
					F_output[1:,1:]=F[:,:]

					f.write('double flux_%s(%d,%d)\n'%(k+1, N+2, M+2))
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

				G_output=np.zeros((N+2,M+2))
				G_output[0,1:]=np.linspace(-180,180,M+1)
				G_output[1:,0]=np.linspace(-180,180,N+1)
				G_output[1:,1:]=G[:,:]

				g.write('# DNI ratio: %s\n'%(d))
				g.write('double mflow_%s(%d,%d)\n'%(dr+1, N+2, M+2))
				for i in range(G_output.shape[0]):
					for j in range(G_output.shape[1]):
						if i==0:
							g.write('%.1f '%G_output[i,j])
						elif j==0:
							g.write('%.1f '%G_output[i,j])
						else:
							g.write('%s '%G_output[i,j])
					g.write('\n')
				g.write('\n')
			g.close()

	def postprocess(self): # OELT and RELT generations
		N=10
		M=24
		T_amb_g=np.linspace(-5.6,45,5)
		V_wind_g=np.linspace(0.,18.,5)
		H=np.zeros((N+1,M+1))

		results_table=np.array(['Qin', 'T_amb', 'Wind_speed', 'DNI','Peak_flux','T_ext_mean','T_ext_mean','h_ext','q_refl','q_emi','q_conv','R_eff'])

		h = open('OELTs_Solstice.motab','w+')
		h.write('#1\n')
		h.write('#Comments\n')
		h.write("#METALABELS,n_helios,A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,CoT1,CoT2,CoT3,CoT4,CoT'1,CoT'2,CoT'3,CoT'4,Coh1,Coh2,Coh3,Coh4,Coh5,eff_abs,eff_emi,SM\n")
		h.write('#METAUNITS,integer,m2,real,m,m,m,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real\n')
		for i,n,m in zip(self.jobs,self.n,self.m):
			casedir = os.path.join(os.getcwd(),'dnir1.0','job{0}'.format(i))
			if n==3 or n==5 or n==7:
				results_table=np.vstack((results_table,np.genfromtxt('%s/RELT.csv'%casedir, delimiter=',')))

		Equinox=np.loadtxt('equinox/Equinox.csv',delimiter=',')
		np.savetxt('RELT.csv', results_table, fmt='%s', delimiter=',')
		coefs_T,coefs,eff_abs,eff_emi,A_rec=receiver_correlation(self.r_diameter,self.r_height,folder='.')
		h.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(
			int(self.num_hst), self.hst_h*self.hst_w, Equinox[0], self.r_diameter, self.r_height, self.tower_h, Equinox[1],
			coefs_T[0],coefs_T[1][0],coefs_T[1][1],coefs_T[1][2],coefs_T[2],coefs_T[3][0],coefs_T[3][1],coefs_T[3][2],
			coefs[0],coefs[1],coefs[2],coefs[3],coefs[4],eff_abs,eff_emi,self.SM))

		for dr,d in enumerate(self.dnir):
			print('	Optics d: ({0:.2f})'.format(d))
			for i,n,m in zip(self.jobs,self.n,self.m):
				casedir = os.path.join(os.getcwd(),'dnir{0}'.format(d),'job{0}'.format(i))
				fileo = open('%s/flux-table'%casedir,'rb')
				data = pickle.load(fileo)
				fileo.close()
				dni=self.dni[i]
				if abs(data['T_HC'][0][449]-data['T_out'])<1.0 and abs(data['T_HC'][1][449]-data['T_out'])<1.0:
					eta,q_results,eta_exc_intec=proces_raw_results('%s/simul'%casedir,'%s/'%casedir)
					GA_tot = d*dni*self.num_hst*self.hst_w*self.hst_h
					H[n,m] = q_results[-1]/GA_tot

			for n in range(3):
				H[n,:]=H[5-n,:]

			for n in range(8,11):
				H[n,:]=H[15-n,:]

#			H[:,13:M+1]=H[:,0:12][:,::-1]

			H_output=np.zeros((N+2,M+2))
			H_output[0,1:]=np.linspace(-180,180,M+1)
			H_output[1:,0]=np.linspace(-180,180,N+1)
			H_output[1:,1:]=H[:,:]

			h.write('# DNI ratio: %s\n'%(d))
			h.write('double optics_%s(%d,%d)\n'%(dr, N+2, M+2))
			for i in range(H_output.shape[0]):
				for j in range(H_output.shape[1]):
					if i==0:
						h.write('%.1f '%H_output[i,j])
					elif j==0:
						h.write('%.1f '%H_output[i,j])
					else:
						h.write('%s '%H_output[i,j])
				h.write('\n')
			h.write('\n')
		h.close()


if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Miscelaneous utilities.')
	parser.add_argument('--T', type=int, default=565, help='')
	parser.add_argument('--folder', type=str, default='./gemasolar_annual_N08810_565_elastic', help='')

	args = parser.parse_args()
	model=postprocess()
	model.postprocess()
	model.get_fluxes()
