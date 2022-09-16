#! /bin/env python3

import os
import shutil
import numpy as np
import argparse
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

class SaltAnnualTrimmedField:
	def __init__(self,material='N08810',D0=45.0,WT=1.5,T=565):

		matdict = {'N06230':'Haynes230'
			,'N08810':'Incoloy800H'
			,'N07740':'Inconel740H'}

		case = '%s_OD%.2f_WT%.2f_%s'%(material,D0,WT,T)
		print('	',case)
		casedir=os.path.abspath(
			os.path.join(os.path.dirname(__file__), case))
		trimmed_field = os.path.join(
			os.path.abspath(os.path.dirname(__file__)),'pos_and_aiming_trimmed.csv')

		if not os.path.exists(casedir):
			os.makedirs(casedir)
			shutil.copy(trimmed_field,os.path.join(casedir,'pos_and_aiming_trimmed.csv'))

		Model=one_key_start(
			casedir=casedir, 
			tower_h=114.75, 
			Q_rec=111.e6/0.51*2.4,
			T_in=290+273.15,
			T_out=T+273.15,
			HTF='salt',
			rec_material=matdict[material],
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
		Model.flow_path_salt(num_bundle=18,num_fp=2,D0=D0,WT=WT,pattern='NES-NWS')
		vfs = np.array([0.1, 1., 2., 3., 4.])
		T_int = np.linspace(290.,565.,12)
		T_int = np.append(T_int,600.) + 273.15
		ThermoElasticPeakFlux.fluxLim(Model.D0,Model.WT,os.path.join(casedir,Model.material_name),Model.mat,vfs,T_int)
		Model.annual_trimmed_field()
		Model.get_fluxes()
		Model.get_mflow()
		Model.Tables()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run annual simulation')
	parser.add_argument('--material', type=str, default='N08810')
	parser.add_argument('--OD', type=float, default=45.0)
	parser.add_argument('--WT', type=float, default=1.5)
	parser.add_argument('--TF', type=float, default=565)
	args = parser.parse_args()

	SaltAnnualTrimmedField(material=args.material,D0=args.OD,WT=args.WT,T=args.TF)

