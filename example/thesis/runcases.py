#! /bin/env python3

import os
import shutil
import numpy as np

FILEDIR=os.getcwd()

data=np.genfromtxt(os.path.join(FILEDIR,'cases.csv'),delimiter=',')

material='N06230'
diameter=22.4
thickness= 1.2
temperature=565
mydict={'N06230':'A230','N08810':'A800'}

MATDIR=os.path.join('/scratch/xa1/af5590/mdba_jobs','fatigue7') # CHANGE THIS FOLDER NAME
sf_vector='1.2 1.2 1.2 0.9725 0.8375 0.8475 0.8125 0.8125 0.7575'

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
		s+='#PBS -P xa1\n'
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
		s+='export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/xa1/software-package/om-inst-v1.14.2/lib\n'
		s+='export PATH=$PATH:/scratch/xa1/software-package/om-inst-v1.14.2/bin\n'
		s+='export PYTHONPATH=$PYTHONPATH:/scratch/xa1/af5590/st-inst-af/lib/python3.8/site-packages\n'
		s+='export OPENMODELICALIBRARY=/scratch/xa1/af5590/st-inst-af/lib/omlibrary:/scratch/xa1/software-package/om-inst-v1.14.2/lib/omlibrary/\n'
		s+='\n'
		s+='source /scratch/xa1/software-package/solstice-0.9.0/etc/solstice.profile\n'
		s+='\n'
		s+='cd %s\n'%CASEDIR
		s+='python3.8 %s/runGemasolar.py --material %s --D0 %s --WT %s --T %d --dnir %s --case %d  --sf_vector %s'%(FILEDIR,material,diameter,thickness,temperature,dnir,case,sf_vector)

		if not os.path.exists(CASEDIR):
			os.mkdir(CASEDIR)

		os.system('cp %s/pos_and_aiming_trimmed.csv %s/'%(FILEDIR, CASEDIR))
#		os.system('cp %s/%s_OD%.2f_WT%.2f_peakFlux.csv %s/'%(FILEDIR,material,diameter,thickness,CASEDIR))

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
s+='#PBS -P xa1\n'
s+='#PBS -q normal\n'
s+='#PBS -l walltime=02:00:00,mem=16GB,ncpus=8\n'
s+='#PBS -N runGemasolar\n'
s+='#PBS -o runGemasolar.out\n'
s+='#PBS -e runGemasolar.err\n'
s+='\n'
s+='echo ""\n'
s+='echo "start job script"\n'
s+='\n'
s+='module load python3/3.8.5\n'
s+='\n'
s+='export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/xa1/software-package/om-inst-v1.14.2/lib\n'
s+='export PATH=$PATH:/scratch/xa1/software-package/om-inst-v1.14.2/bin\n'
s+='export PYTHONPATH=$PYTHONPATH:/scratch/xa1/af5590/st-inst-af/lib/python3.8/site-packages\n'
s+='export OPENMODELICALIBRARY=/scratch/xa1/af5590/st-inst-af/lib/omlibrary:/scratch/xa1/software-package/om-inst-v1.14.2/lib/omlibrary/\n'
s+='\n'
s+='source /scratch/xa1/software-package/solstice-0.9.0/etc/solstice.profile\n'
s+='\n'
s+='cd %s\n'%CASEDIR
s+='python3.8 %s/runEquinox.py --material %s --D0 %s --WT %s --T %d --dnir 1.0 --case 0  --sf_vector %s'%(FILEDIR,material,diameter,thickness,temperature,sf_vector)

f=open(os.path.join(CASEDIR,'jobscript'),'w+')
f.write(s)
f.close()

for dni,dnir in enumerate([0.56,0.87,1.0,1.2,1.39]):
	DNIDIR=os.path.join(MATDIR,'dnir%s'%dnir)
	for case in range(data.shape[0]):
		CASEDIR=os.path.join(DNIDIR,'job%d'%case)
		os.chdir(CASEDIR)
		os.system('qsub -N %sT%dD%sJ%d jobscript'%(mydict[material],temperature,dni,case)) # CHANGE THIS JOB NAME

CASEDIR=os.path.join(MATDIR,'equinox')
os.chdir(CASEDIR)
os.system('qsub -N equinox jobscript') # CHANGE THIS JOB NAME

