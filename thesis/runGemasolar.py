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

LIMITS_DIR=os.path.dirname(__file__)

class SaltAnnualTrimmedField:
	def __init__(self,material='N08810',D0=45.0,WT=1.5,T=565,index=0,dnir=1.0,sf_vector=[1.,1.,1.,0.97,0.85,0.85,0.81,0.81,0.78]):

		matdict = {'N06230':'Haynes230'
			,'N08810':'Incoloy800H'
			,'N07740':'Inconel740H'}

		case = '%s_OD%.2f_WT%.2f_%d'%(material,D0,WT,T)
		print('	Case is: ',case,'	Index:',index)
		casedir=os.getcwd()

		print(dnir)
		T_amb_g=np.linspace(-5.6,45,5)
		V_wind_g=np.linspace(0.,18.,5)
		results_table=np.array([]) # RELT

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
			fluxlimitpath=LIMITS_DIR,
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
			num_rays=int(1e7),
			latitude=37.56,
			sf_vector=sf_vector
			)
		print('	AFD SF: %s'%Model.sf_vector)

		print('	Equivalent slope error: %.2f [mrad]'%(0.5*1e3*np.sqrt(4*pow(2.6e-3,2)+pow(2.1e-3,2))))
		# input the number of tube bundles, number of flowpaths, pipe outer diameter and flow path pattern
		Model.flow_path_salt(num_bundle=18,num_fp=2,D0=D0,WT=WT,pattern='NES-NWS')
		vfs = np.array([0.1, 1., 2., 3., 4.])
		T_int = np.linspace(290.,565.,12)
		T_int = np.append(T_int,600.) + 273.15
		ThermoElasticPeakFlux.fluxLim(Model.D0,Model.WT,os.path.join(LIMITS_DIR,Model.material_name),Model.mat,vfs,T_int)
		data = np.genfromtxt(os.path.join(os.path.dirname(__file__),'cases.csv'),delimiter=',')
		Model.MDBA_aiming_new(dni=dnir*data[index,5],phi=data[index,6],elevation=data[index,7])
		# RLLT
		n = data[index,3]
		m = data[index,4]
		if dnir==1:
			if n==3 or n==5 or n==7:
				dni=dnir*data[index,5]
				max_flux=read_data(casedir,Model.r_height,Model.r_diameter,Model.num_bundle,Model.bins,flux_file=True)
				for i in range(int(len(T_amb_g))):
					print('T_abm_g, 3', T_amb_g[i],3.)
					results,aiming_results,vel_max,Strt=Model.HT_model(T_amb_g[i],3.)
					if np.isnan(results[3]):
						continue 
					results_table=np.append(results_table,[results[0],T_amb_g[i],3.,dni,max_flux])
					results_table=np.append(results_table,results[3:])

				for i in range(int(len(V_wind_g))):
					print('20 V_wind_g', 20.,V_wind_g[i])
					results,aiming_results,vel_max,Strt=Model.HT_model(20.,V_wind_g[i])
					if np.isnan(results[3]):
						continue 
					results_table=np.append(results_table,[results[0],20.,V_wind_g[i],dni,max_flux])
					results_table=np.append(results_table,results[3:])
				results_table=results_table.reshape(int(len(results_table)/12), 12)
				np.savetxt('%s/RELT.csv'%casedir, results_table, fmt='%s', delimiter=',')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Run annual simulation')
	parser.add_argument('--material', type=str, default='N06230') # CHANGE THIS MATERIAL NAME
	parser.add_argument('--D0', type=float, default=42.16)        # CHANGE THIS MATERIAL NAME
	parser.add_argument('--WT', type=float, default=1.2)          # CHANGE THIS MATERIAL NAME
	parser.add_argument('--T', type=int, default=600)             # CHANGE THIS TEMPERATURE
	parser.add_argument('--case', type=int, default=0)
	parser.add_argument('--dnir', type=float, default=1.0)
	parser.add_argument('--sf_vector', type=float, nargs=9, default=[1,1,1,1,1,1,1,1,1])
	args = parser.parse_args()

	SaltAnnualTrimmedField(material=args.material,D0=args.D0,WT=args.WT,T=args.T,index=args.case,dnir=args.dnir,sf_vector=args.sf_vector)

