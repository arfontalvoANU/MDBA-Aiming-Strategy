#! /bin/env python3

from __future__ import division
import unittest

import os
import numpy as np
import pickle
from scipy.interpolate import interp2d
import nitrateSaltPeakFlux as nspf
import PostProcessing as pp

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

		casedir=os.path.abspath(
			os.path.join(os.path.dirname(__file__), 'TEST-COOPTIMISATION-SALT'))
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
			num_rays=int(5e6),
			latitude=37.56,
			)
		print('	Equivalent slope error: %.2f [mrad]'%(0.5*1e3*np.sqrt(4*pow(2.6e-3,2)+pow(2.1e-3,2))))
		# input the number of tube bundles, number of flowpaths, pipe outer diameter and flow path pattern
		Model.flow_path_salt(num_bundle=18,num_fp=2,D0=22.4,WT=1.2,pattern='NES-NWS')
		nspf.fluxLim(Model.D0,Model.WT,os.path.join(casedir,Model.material_name),Model.mat)
		Model.annual_trimmed_field([1.00,1.39,0.87,0.56])
		res=pp.PostProcess([1.00,1.39,0.87,0.56],num_fp=2,num_el=450,casedir=casedir,latitude=Model.latitude)
		res.get_fluxes()
		res.get_mflow()
		res.Tables()

	def test_touching(self):
		if os.path.exists(self.tablefile):
			successed=1
		self.assertEqual(successed,1)
		#os.system('rm *.vtk')


if __name__ == '__main__':
	unittest.main()

