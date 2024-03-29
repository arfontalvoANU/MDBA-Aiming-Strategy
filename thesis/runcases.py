#! /bin/env python3

import os
import shutil
import argparse
import numpy as np

def run_cases(args):
	FILEDIR=os.path.dirname(os.path.abspath(__file__))

	data=np.genfromtxt(os.path.join(FILEDIR,'cases.csv'), delimiter=',')

	mydict={'N06230':'A230','N08810':'A800'}
	s = ''
	for item in args.sf_vector:
		s+=f'{item} '
	args.sf_vector = s

	MATDIR=os.path.join('/scratch',args.project,'af5590','mdba_jobs',args.casename) # CHANGE THIS FOLDER NAME

	if not os.path.exists(MATDIR):
		os.mkdir(MATDIR)

	for dnir in [0.56,0.87,1.0,1.2,1.39]:
		DNIDIR=os.path.join(MATDIR,'dnir%s'%dnir)
		if not os.path.exists(DNIDIR):
			os.mkdir(DNIDIR)

		for case in range(data.shape[0]):

			CASEDIR=os.path.join(DNIDIR,'job%d'%case)

			s='#!/bin/sh\n'
			s+='#PBS -S /bin/sh\n'
			s+=f'#PBS -P {args.project}\n'
			s+='#PBS -q normal\n'
			s+='#PBS -l walltime=03:00:00,mem=8GB,ncpus=4\n'
			s+='#PBS -N runGemasolar\n'
			s+='#PBS -o runGemasolar.out\n'
			s+='#PBS -e runGemasolar.err\n'
			s+='\n'
			s+='echo ""\n'
			s+='echo "start job script"\n'
			s+='\n'
			s+='module load python3/3.8.5\n'
			s+='\n'
			s+=f'source /scratch/{args.project}/software-package/solstice-0.9.0/etc/solstice.profile\n'
			s+='\n'
			s+='cd %s\n'%CASEDIR
			s+=f'python3.8 {FILEDIR}/runGemasolar.py --material {args.material} --D0 {args.D0} --WT {args.WT} --T {args.T:d} --dnir {dnir} --case {case}  --sf_vector {args.sf_vector}'

			if not os.path.exists(CASEDIR):
				os.mkdir(CASEDIR)

			os.system('cp %s/pos_and_aiming_trimmed.csv %s/'%(FILEDIR, CASEDIR))

			f=open(os.path.join(CASEDIR,'jobscript'),'w+')
			f.write(s)
			f.close()

	# Equinox
	CASEDIR=os.path.join(MATDIR,'equinox')
	if not os.path.exists(CASEDIR):
		os.mkdir(CASEDIR)

	os.system('cp %s/pos_and_aiming_trimmed.csv %s/'%(FILEDIR, CASEDIR))

	s='#!/bin/sh\n'
	s+='#PBS -S /bin/sh\n'
	s+=f'#PBS -P {args.project}\n'
	s+='#PBS -q normal\n'
	s+='#PBS -l walltime=03:00:00,mem=8GB,ncpus=4\n'
	s+='#PBS -N runGemasolar\n'
	s+='#PBS -o runGemasolar.out\n'
	s+='#PBS -e runGemasolar.err\n'
	s+='\n'
	s+='echo ""\n'
	s+='echo "start job script"\n'
	s+='\n'
	s+='module load python3/3.8.5\n'
	s+='\n'
	s+=f'source /scratch/{args.project}/software-package/solstice-0.9.0/etc/solstice.profile\n'
	s+='\n'
	s+='cd %s\n'%CASEDIR
	s+=f'python3.8 {FILEDIR}/runEquinox.py --material {args.material} --D0 {args.D0} --WT {args.WT} --T {args.T:d} --dnir 1.0 --case 0  --sf_vector {args.sf_vector}'

	f=open(os.path.join(CASEDIR,'jobscript'),'w+')
	f.write(s)
	f.close()

	if args.submit:
		for dni,dnir in enumerate([0.56,0.87,1.0,1.2,1.39]):
			DNIDIR=os.path.join(MATDIR,'dnir%s'%dnir)
			for case in range(data.shape[0]):
				CASEDIR=os.path.join(DNIDIR,'job%d'%case)
				os.chdir(CASEDIR)
				os.system(f'qsub -N {args.casename[0]}{args.casename[-2::]}-d{dnir}-j{case} jobscript') # CHANGE THIS JOB NAMEs

		CASEDIR=os.path.join(MATDIR,'equinox')
		os.chdir(CASEDIR)
		os.system(f'qsub -N {args.casename[0]}{args.casename[-2::]}-eqnx jobscript') # CHANGE THIS JOB NAME

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to prepare folders and jobscripts for independent simulations')
	parser.add_argument('casename', type=str, default='fatigue6')
	parser.add_argument('--material', type=str, default='N06230')
	parser.add_argument('--project', type=str, default='xa1')
	parser.add_argument('--D0', type=float, default=22.4)
	parser.add_argument('--WT', type=float, default=1.2)
	parser.add_argument('--T', type=int, default=565)
	parser.add_argument('--sf_vector', type=float, nargs=9, default=[1,1,1,0.9725,0.8375,0.8475,0.8125,0.8125,0.7575])
	parser.add_argument('--submit',action='store_true')
	parser.add_argument('--no-submit', dest='submit', action='store_false')
	parser.set_defaults(submit=True)

	args = parser.parse_args()
	run_cases(args)
