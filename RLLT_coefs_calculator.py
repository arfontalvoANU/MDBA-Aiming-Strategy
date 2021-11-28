'''
This code uses a loss analysis method to calculate the receiver losses 
and effiency in an annual simulation.
'''

import os
from sys import path
import numpy as np
from pandas import DataFrame
from sklearn import linear_model
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

def Loss_analysis(fpath):
	# read from input file
	abs_t=0.98 
	eff_abs=abs_t/(2./np.pi*(1.-abs_t)+abs_t) # effective absorptivity of the pipe
	ems_t=0.91
	eff_emi=ems_t/(2./np.pi*(1.-ems_t)+ems_t) # effective emissivity of the pipe

	R_input=np.loadtxt('%s/RLLT_input.csv' % fpath   ,delimiter=',', skiprows=1)
	R_output_old=np.loadtxt('%s/RLLT.csv' % fpath   ,delimiter=',', skiprows=1)
#	R_input=np.genfromtxt(RLLT_input, delimiter=',', skip_header=1)
#	R_output_old=np.genfromtxt(RLLT, delimiter=',', skip_header=1)
	R_output=R_output_old[~np.isnan(R_output_old).any(axis=1)] # to remove the rows with nan
	R_input=R_input[~np.isnan(R_output_old).any(axis=1)]

	R_input=R_input[~np.isnan(R_output).any(axis=1)]
	Q_in=R_output[:,-2]/R_output[:,-1]/eff_abs/1e6  # MW
	T_amb=R_input[:,-2] # C
	V_wind=R_input[:,-1] # m/s
	R_eff=R_output[:,-1]*eff_abs

	A_rec=R_output[0,3]/(R_output[0,0]*(R_output[0,1]-(R_input[0,-2]+273.15))) # receiver surface area

	# the reflecive loss
	Q_ref = (1-eff_abs)*Q_in

	# linear regression for T_ext
	T_ext={'X2': Q_in, 'X3': (T_amb+273.15),'X4': V_wind,'T_ext': R_output[:,1]}
	df = DataFrame(T_ext,columns=['X2','X3','X4','T_ext'])
	X = df[['X2','X3','X4']] 
	Y = df['T_ext']
	regr = linear_model.LinearRegression()
	regr.fit(X, Y)
	C0=regr.intercept_
	C=regr.coef_
	T_ext_linear=C0+C[0]*Q_in+C[1]*(T_amb+273.15)+C[2]*V_wind

	# linear regression for T_ext_4_mean
	T_ext_4={'X2': Q_in, 'X3': (T_amb+273.15),'X4': V_wind,'T_ext_4': R_output[:,2]}
	df = DataFrame(T_ext_4,columns=['X2','X3','X4','T_ext_4'])
	X = df[['X2','X3','X4']] 
	Y = df['T_ext_4']
	regr = linear_model.LinearRegression()
	regr.fit(X, Y)
	C1=regr.intercept_
	C1_1=regr.coef_
	T_ext_4_linear=C1+C1_1[0]*Q_in+C1_1[1]*(T_amb+273.15)+C1_1[2]*V_wind

	coefs_T=[C0,C,C1,C1_1]

	# how to calculate h_conv
	coefs=np.polyfit(R_input[:,-1],R_output[:,0],4)
	h_conv=coefs[4]+coefs[3]*R_input[:,-1]+coefs[2]*(R_input[:,-1])**2+coefs[1]*(R_input[:,-1])**3+coefs[0]*(R_input[:,-1])**4 # h with polynominal fitting

	#print max(abs(h_conv-R_output[:,0])/R_output[:,0])
	Q_emi=eff_emi*5.67e-8*A_rec*(T_ext_4_linear**4-(T_amb+273.15)**4)/1e6

	#print abs(sum(Q_emi)-sum(R_output[:,4])/1e6)/sum(R_output[:,4]/1e6)

	Q_conv=h_conv*A_rec*(T_ext_linear-(T_amb+273.15))/1e6
	#print abs(sum(Q_conv)-sum(R_output[:,3])/1e6)/sum(R_output[:,3]/1e6)

	Qnet=eff_abs*Q_in-Q_conv-Q_emi
	Qnet_real= R_output[:,-2]/1e6
	Diff=(Qnet-Qnet_real)/Qnet_real
	#print max(abs(Qnet-Qnet_real)/Qnet_real)
	#print abs(sum(Qnet)-sum(Qnet_real))/sum(Qnet_real)

	return coefs_T,coefs,eff_abs,eff_emi,A_rec

if __name__=='__main__':
	currentdir=os.getcwd()

	fpath = '%s/example/case_H230_Sodium_290_565'%currentdir

	headers = \
		'coefs_T[0],coefs_T[1][0],coefs_T[1][1],coefs_T[1][2],'+\
		'coefs_T[2],coefs_T[3][0],coefs_T[3][1],coefs_T[3][2],'+\
		'coefs[0],coefs[1],coefs[2],coefs[3],coefs[4],eff_abs,eff_emi\n'

	RLLTS = open('%s/example/case_H230_Sodium_290_565/coefs.txt'%currentdir,'a')
	RLLTS.write(headers)

	coefs_T,coefs,eff_abs,eff_emi,A_rec=Loss_analysis(fpath)
	RLLTS.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(
		coefs_T[0],
		coefs_T[1][0],
		coefs_T[1][1],
		coefs_T[1][2],
		coefs_T[2],
		coefs_T[3][0],
		coefs_T[3][1],
		coefs_T[3][2],
		coefs[0],
		coefs[1],
		coefs[2],
		coefs[3],
		coefs[4],
		eff_abs,
		eff_emi)
		)
	RLLTS.close()
