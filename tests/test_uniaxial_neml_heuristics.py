#!/usr/bin/env python3
from time import time
from mdbapy.nitrateSaltHeuristics import *

tinit = time.time()

cases = {
	 'N06230_OD22.40_WT1.20_565':{'p':6,'z':5.67}
	,'N06230_OD22.40_WT1.20_600':{'p':7,'z':5.67}
	,'N06230_OD45.00_WT1.50_565':{'p':1,'z':5.67}
	,'N06230_OD45.00_WT1.50_600':{'p':1,'z':4.41}
	}

for case in cases:
	if 'N08810' in case:
		mat = '800H'
	else:
		mat = 'A230'
	casedir = '/home/arfontalvo/ownCloud/phd_update/damage/N06230/MDBA2ST_5DNIR_53SUNPOS/HEURISTICS-SALT/'
	filename = '{0}.mat'.format(case)
	p=cases[case]['p']
	z=int(cases[case]['z']/10.5*50)
	print('	Control')
	print('	material: {0}	panel: {1:g}	height:{2:g}'.format(mat,p,z))
	model = heuristics(res=filename,folder=casedir,mat=mat)
	model.uniaxial_neml(p,z)

seconds = time.time() - tinit
m, s = divmod(seconds, 60)
h, m = divmod(m, 60)
print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
