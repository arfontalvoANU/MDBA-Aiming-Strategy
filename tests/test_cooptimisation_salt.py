#! /bin/env python3

from __future__ import division
import unittest

import os
import numpy as np


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

		casedir='TEST-COOPTIMISATION-SALT'
		self.tablefile=casedir+'/OELT_Solstice.motab'

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
			fluxlimitpath='../data/201015_N08811_thermoElasticPeakFlux_velocity_salt',
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
			num_rays=int(1e6),
			latitude=37.56,
			)
		print('	Equivalent slope error: %.2f [mrad]'%(0.5*1e3*np.sqrt(4*pow(2.6e-3,2)+pow(2.1e-3,2))))
		# input the number of tube bundles, number of flowpaths, pipe outer diameter and flow path pattern
		Model.flow_path_salt(num_bundle=18,num_fp=2,D0=45.,pattern='NES-NWS') 
		Model.MDBA_aiming_new(dni=930.,phi=0.,elevation=75.89)
		C_aiming=np.zeros(Model.num_bundle)
		C_aiming[:]=0.4
		material_name='N08811'
		flux_limits_file='%s/%s_OD%.2f_WT1.50_peakFlux.csv'%(Model.fluxlimitpath,material_name, Model.D0)
		results,aiming_results,vel_max=tower_receiver_plots(
			files=Model.casedir+'/flux-table', 
			efficiency=False, 
			maps_3D=False, 
			flux_map=False, 
			flow_paths=True,
			saveloc=None, 
			billboard=False, 
			flux_limits_file=flux_limits_file,
			C_aiming=C_aiming,overflux=False)

		import pickle
		from scipy.interpolate import interp2d
		fileo = open('%s/flux-table'%(casedir),'rb')
		data = pickle.load(fileo)
		fileo.close()

		fdata = np.loadtxt(flux_limits_file, delimiter=',')
		fit = interp2d(fdata[0,1:], fdata[1:,0], fdata[1:,1:])
		for j in range(2):
			fp = data['fp'][j]
			flux_lim = fit(data['V'][j], data['T_HC'][j])
			np.savetxt('%s/flux_lim_%s.csv'%(casedir,j), flux_lim[:,0]/1e3, delimiter=',')
			table=data['q_net'][fp]/data['areas'][fp]*1e-6
			np.savetxt('%s/flux_mdba_%s.csv'%(casedir,j),table,delimiter=',')

		print('	n_tubes:     %d  '%(data['n_tubes'][0]))
		print('	m_flow fp 1: %.2f'%(data['m'][0]/data['n_tubes'][0]))
		print('	m_flow fp 2: %.2f'%(data['m'][1]/data['n_tubes'][0]))
		print(' ')
		print('	T_out  fp 1: %.2f'%(data['T_HC'][0][-1]-273.15))
		print('	T_out  fp 2: %.2f'%(data['T_HC'][1][-1]-273.15))
		print('	h_conv_ext:  %s  '%data['h_conv_ext'])
		print('	T_ext_ave:   %.2f'%(np.average(data['T_ext'])))

	def test_touching(self):
		if os.path.exists(self.tablefile):
			successed=1
		self.assertEqual(successed,1)
		#os.system('rm *.vtk')


if __name__ == '__main__':
	unittest.main()

