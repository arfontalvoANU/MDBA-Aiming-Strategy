#! /bin/env python3

import os
import shutil
import numpy as np
import colorama

colorama.init()
def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

FILEDIR = os.getcwd()
PROTODIR = os.path.join(FILEDIR,'opticsjob')

data = np.genfromtxt(os.path.join(PROTODIR,'cases.csv'),delimiter=',')
names = os.listdir(os.getcwd())

for name in names:
	MATDIR = os.path.join(FILEDIR,name)
	times = []
	s=0
	for dni,dnir in enumerate([0.56,0.87,1.0,1.2,1.39]):
		DNIDIR = os.path.join(MATDIR,'dnir%s'%dnir)
		for case in range(data.shape[0]):
			CASEDIR = os.path.join(DNIDIR,'job%d'%case)
			try:
				os.chdir(CASEDIR)
				f = open('runGemasolar.out','r')
				p = f.readlines()
				for i in p:
					if 'Walltime Used' in i:
						linev = i.split(' ')
						for jn,j in enumerate(linev):
							if 'Used' in j:
								h = linev[jn+1].split(':')
								h = float(h[0]) + float(h[1])/60 + float(h[2])/3600
								times.append(h)
								s+=h
			except:
				continue
	try:
		print(yellow(name),'Maximum time: {0:.1f} h'.format(max(times)))
	except:
		continue
	print(green(name),'Cummulative time: {0:.1f} h'.format(s))
