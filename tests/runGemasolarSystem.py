#! /bin/env python2

from solartherm import simulation
from solartherm import postproc

import numpy as np
import scipy.io as sio
import DyMat
import os

class runModelica():
	def __init__(self):
		fn = os.path.join(os.path.expanduser('~'),'solartherm/examples/GemasolarSystemOperation.mo')
		self.sim = simulation.Simulator(fn)
		self.sim.compile_model(args=['-d=nonewInst'])
		self.sim.compile_sim(args=['-s'])
		self.sim.load_init()

	def simulate(self,casedir,T,savedir):
		h = 1396.0182*T + 0.086*pow(T,2)
		par_n = [
			'T_hot_set'
			,'T_hot_start'
			,'state_hot_set.h'
			,'controlCold.state_ref.h'
			,'opt_file'
			,'hav_file'
			,'file_dni1'
			,'file_dni2'
			,'file_dni3'
			,'file_dni4'
			,'file_mflow'
		]
		par_v = [
			str(T)
			,str(T)
			,str(h)
			,str(h)
			,os.path.join(casedir,'OELTs_Solstice.motab')
			,os.path.join(casedir,'HALTs_Solstice.motab')
			,os.path.join(casedir,'FLUX_fp1_d0.56.motab')
			,os.path.join(casedir,'FLUX_fp1_d0.87.motab')
			,os.path.join(casedir,'FLUX_fp1_d1.0.motab')
			,os.path.join(casedir,'FLUX_fp1_d1.39.motab')
			,os.path.join(casedir,'MFLOW_Solstice_fp1.motab')
		]
		self.sim.update_pars(par_n,par_v)
		self.sim.simulate(start=0, stop='1y', step='5m', solver='dassl', nls='homotopy')
		self.res = postproc.SimResultElec(self.sim.res_fn)
		self.perf = self.res.calc_perf()

		print('Simulation finished')
		print('Starting post-processing')
		header = ['epy (MWh/year)','lcoe ($/MWh)','capf (%)','srev ($)']
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(header)))
		print(' '.join('%s'%('%s'%x + ' '*(max(len(header[i]),len('%s'%self.perf[i])) - len('%s'%x))) for i,x in enumerate(self.perf)))
		os.system('cp GemasolarSystemOperation*.mat %s'%savedir)

	def clean(self):
		os.system('rm GemasolarSystemOperation*')

if __name__ == '__main__':
	model=runModelica()

	parent=os.path.join(os.getenv('HOME'),'ownCloud/phd_update/damage')
	for material in ['N08810']:
		for case,T in zip(['N08810_OD22.40_WT1.20_600','N08810_OD45.00_WT1.50_600'],[600,600]):
			casedir=os.path.join(parent,material,case)
			print(casedir)
			savedir='/mnt/fb7cc2c9-e328-4f3f-a6f8-918195722408/Gemasolar/%s.mat'%case
			model.simulate(casedir,T,savedir)
	model.clean()
