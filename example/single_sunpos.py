import os
import sys
import inspect
import datetime
from solsticepy.cal_sun import *
import colorama
colorama.init()

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import salt_aiming as aiming
from simulation import optical
from run_solstice import *
from python_postprocessing import *
from flux_reader import *
import pickle
import argparse
import matplotlib.pyplot as plt

def yellow(text):
	return colorama.Fore.YELLOW + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def green(text):
	return colorama.Fore.GREEN + colorama.Style.BRIGHT + text + colorama.Style.RESET_ALL

def simulate_dni_ratio(ratio,C_start,E_start,A_start):

	# define a unique case folder for the user
	basefolder = os.path.join(currentdir,'case-test-N08811')

	if not os.path.exists(basefolder):
		os.makedirs(basefolder)

	folder=basefolder
	source_path=parentdir

	# Inputs (Gemasolar)
	num_bundle=18         # Number of panels
	r_height=10.5         # Receiver height [m]
	r_diameter=8.5        # Receiver diameter [m]
	bins=50               # Vertical bins
	tower_h=114.75        # Tower height [m] (Floor to receiver bottom)
	azimuth=0.0           # Solar azimuth angle [degrees]
	elevation=55.15       # Solar elevation angle [degrees]
	DNI=930.0             # Beam irradiance [W/m2]
	num_fp=2              # Number of flow path
	D0=45.00              # Panel tube OD [mm]
	num_rays=int(5e6)     # Number of rays

	# Creating receiver model
	Model=aiming.one_key_start(
		folder,
		source_path,
		num_bundle,
		num_fp,
		r_height,
		r_diameter,
		bins,
		tower_h,
		azimuth,
		elevation,
		DNI,
		D0)
	Model.sweeping_algorithm(C_start,E_start,A_start)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Run representative sun positions of annual simulation for a specific DNI ratio')
	parser.add_argument('--ratio', type=int, default=1, help='DNI ratio to be simulated. Default=1')
	parser.add_argument('--fpath', type=int, default=1, help='Flowpath data to be compiled. Default=1')
	parser.add_argument('--C_start', type=float, default=0.5, help='The starting value of the aiming extent')
	parser.add_argument('--E_start', type=float, default=2.0, help='The starting value of the aiming exponent')
	parser.add_argument('--A_start', type=float, default=0.5, help='The starting value of the aiming asymetry factor')
	args = parser.parse_args()

	simulate_dni_ratio(args.ratio, args.C_start, args.E_start, args.A_start)

