#! /bin/env python3

import os
import shutil
import numpy as np
import colorama
import time

colorama.init()
def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

FILEDIR = os.getcwd()
data = np.genfromtxt(os.path.join(FILEDIR,'cases.csv'),delimiter=',')

#FILEDIR = '/scratch/xa1/af5590/mdba_jobs'
#MATDIR = os.path.join(FILEDIR,'N06230_OD22.40_WT1.20_565')
MATDIR=FILEDIR
for dni,dnir in enumerate([0.56,0.87,1.0,1.2,1.39]):
	DNIDIR = os.path.join(MATDIR,'dnir%s'%dnir)
	for case in range(data.shape[0]):
		CASEDIR = os.path.join(DNIDIR,'job%d'%case);print(yellow('%s'%CASEDIR))
		f=open(os.path.join(CASEDIR,'runGemasolar.err'),'r')
		l=f.readlines()
		for i in l:
			print(i)
#			if 'Walltime Used' in i:
#				print(yellow(i))#print(yellow('%s'%os.path.join(CASEDIR,'runGemasolar.err')))
#				os.system('cat %s'%os.path.join(CASEDIR,'runGemasolar.err'))

            #print(yellow('%s'%os.path.join(CASEDIR,'runGemasolar.err')))
            #os.system('cat %s'%os.path.join(CASEDIR,'runGemasolar.err'))
            #time.sleep(0.1)
