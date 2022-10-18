#!/usr/bin/env python3
from time import time
from mdbapy.nitrateSaltHeuristics import *

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Estimates average damage of a representative tube in a receiver panel')
	# Heuristics
	parser.add_argument('--days', nargs=2, type=int, default=[0,1])
	parser.add_argument('--step', type=float, default=1800)
	parser.add_argument('--heuristics', type=int, default=1)
	parser.add_argument('--file', type=str, default=os.path.join(os.getenv('HOME'),'solartherm/examples/GemasolarSystemOperation_res.mat'))
	parser.add_argument('--case', type=str, default='heuristics_res')
	parser.add_argument('--mat', type=str, default='800H')
	parser.add_argument('--DO', type=float, default=22.4)
	parser.add_argument('--WT', type=float, default=1.20)
	parser.add_argument('--casedir', type=str, default=os.path.join(os.getcwd(),'HEURISTICS-SALT'))
	args = parser.parse_args()
	tinit = time.time()

	# Material
	materials = {
		'A230':{'kp':17.97,'alpha':15.61e-06,'Young':184e9,'poisson':0.31},
		'800H':{'kp':18.30,'alpha':18.28e-06,'Young':171e9,'poisson':0.31}
	}

	# Heuristics
	model = receiver_cyl(
		Ri = args.DO/2000.0 - args.WT/1000.0
		,Ro = args.DO/2000.0
		,R_fouling=8.808e-5
		,alpha = materials[args.mat]['alpha']
		,Young = materials[args.mat]['Young']
		,poisson = materials[args.mat]['poisson']
		,ab = 0.93
		,em = 0.87
		,kp = materials[args.mat]['kp']
		,Dittus=False
		,mat=args.mat)

	args.case=args.file.split('/')[-1].split('.mat')[0]
	print('	Ri: %.2f'%(model.Ri*1000),
		'	Ro: %.2f'%(model.Ro*1000),
		'	alpha: %g'%model.l,
		'	Young: %.1f'%(model.E/1e9),
		'	poisson: %.2f'%model.nu,
		'	kp: %s'%model.kp)
	print('	SolarTherm res: %s'%args.file)
	print('	Material: %s'%args.mat)
	print('	Case: %s'%args.case)
	print('	Casedir: %s'%args.casedir)
	model.run_heuristics(days=args.days,step=args.step,rfile=args.file,case=args.case,casedir=args.casedir)

	seconds = time.time() - tinit
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	print('Simulation time: {:d}:{:02d}:{:02d}'.format(int(h), int(m), int(s)))
