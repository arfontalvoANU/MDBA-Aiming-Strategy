#! /bin/env python3

from __future__ import division
import unittest

import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

class TestCooptimisationSalt(unittest.TestCase):
	def setUp(self):

		from mdbapy.one_key_co_optimisation import one_key_start
		from mdbapy.cal_sun import SunPosition
		from mdbapy.cal_layout_r import radial_stagger, aiming_cylinder
		from mdbapy.Deviation_aiming_new3 import aiming
		from mdbapy.Open_CSPERB import eval_v_max, Cyl_receiver
		from mdbapy.Open_CSPERB_plots import tower_receiver_plots
		from mdbapy.HC import Na
		from mdbapy.Tube_materials import Inconel740H
		from mdbapy.Flux_reader import read_data
		from mdbapy.Loss_analysis import receiver_correlation
		from mdbapy.output_motab import output_motab, output_matadata_motab
		from mdbapy.python_postprocessing import proces_raw_results, get_heliostat_to_receiver_data
		from mdbapy.SOLSTICE import SolsticeScene
		from mdbapy import ThermoElasticPeakFlux

		casedir=os.path.abspath(
			os.path.join(os.path.dirname(__file__), 'TEST-COOPTIMISATION-SALT'))

		if not os.path.exists(casedir):
			os.makedirs(casedir)

		Model=one_key_start(
			casedir=casedir, 
			tower_h=114.75, 
			Q_rec=111.e6/0.51*2.4,
			T_in=290+273.15,
			T_out=565+273.15,
			HTF='salt',
			rec_material='Incoloy800H',
			r_diameter=8.5,
			r_height=10.5,
			fluxlimitpath=casedir,
			SM=2.4,
			oversizing=1., 	
			delta_r2=0.9,
			delta_r3=1.9,
			hst_w=10.5623289,
			hst_h=10.5623289,
			mirror_reflectivity=0.88*0.95,
			slope_error=0.5*np.sqrt(4*pow(2.6e-3,2)+pow(2.1e-3,2)),
			sunshape='pillbox',
			sunshape_param=np.degrees(4.65e-3),
			num_rays=int(10e6),
			latitude=37.56,
			)
		print('	Equivalent slope error: %.2f [mrad]'%(0.5*1e3*np.sqrt(4*pow(2.6e-3,2)+pow(2.1e-3,2))))
		# input the number of tube bundles, number of flowpaths, pipe outer diameter and flow path pattern
		Model.flow_path_salt(num_bundle=18,num_fp=2,D0=45.0,WT=1.5,pattern='NES-NWS')
		Model.MDBA_aiming_new(dni=930.,phi=0.,elevation=75.89)

		fileo = open(os.path.join(casedir,'flux-table'),'rb')
		data = pickle.load(fileo)
		fileo.close()

		fp = data['fp']
		q_net = data['q_net']
		areas = data['areas']
		pipe_lengths = data['pipe_lengths']

		print('	m_flow fp 1: %.2f'%(data['m'][0]/data['n_tubes'][0]))
		print('	m_flow fp 2: %.2f'%(data['m'][1]/data['n_tubes'][0]))

		x = pipe_lengths[0]
		x = 0.5*(x[1:] + x[:-1])

		flux=q_net[fp[0]]/areas[fp[0]]/1e6

		sanchez_file = '/home/arfontalvo/ownCloud/phd_update/damage/gemasolar_verification/sanchez_flux.csv'
		sanchez=np.genfromtxt(sanchez_file,delimiter=',',skip_header=1)

		armando_file = '/home/arfontalvo/ownCloud/phd_update/damage/gemasolar_verification/mdba_flux_565.csv'
		armando=np.genfromtxt(armando_file,delimiter=',',skip_header=1)

		lw=1.0
		fig,axes = plt.subplots(1,1, figsize=(6,3))
		axes.plot(sanchez[:,0],sanchez[:,1],lw=lw,label='Sanchez')
		axes.plot(armando[:,0],armando[:,1],lw=lw,label='MDBA Armando')
		axes.plot(x,flux,lw=lw,label='MDBA Shuang')
		axes.set_ylim([0,1])
		axes.set_yticks(np.linspace(0,1,11))
		axes.set_xlim([0,94.5])
		axes.set_xticks(np.linspace(0,94.5,10))
		axes.legend(loc='best')
		plt.savefig(os.path.join(casedir,'verification.png'),dpi=300)

	def test_touching(self):
		successed=1
		self.assertEqual(successed,1)
		#os.system('rm *.vtk')


if __name__ == '__main__':
	unittest.main()

