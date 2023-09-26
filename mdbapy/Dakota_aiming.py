from onekey_aiming import one_key_start
from sys import path
import numpy as np

def dakota_opt(**kwargs):
	num_fns = kwargs['functions']  # equaling to one corresponds to the 
	x = kwargs['cv']
	ASV = kwargs['asv']
	retval = dict([])
	
	from sys import path
	folder=path[0]
	num_bundle_2=num_bundle=16
	r_height_2=r_height=24.0
	r_diameter_2=r_diameter=16.0
	bins=50
	tower_h_2=tower_h=175.0
	phi_1=phi=0.0
	elevation_1=elevation=55.15
	DNI_1=DNI=980.0
	num_fp_2=num_fp=8
	D0_2=D0=60.33
	Model=one_key_start(folder,num_bundle,num_fp,r_height,r_diameter,bins,tower_h,phi,elevation,DNI,D0)
	
	# initialization
	C_aiming=np.zeros(num_bundle)
	Exp=np.zeros(num_bundle)
	A_f=np.zeros(num_bundle)
	
	# flowpath 12
	C_aiming_12_1=C_aiming[0]=0.7000000000000002
	Exp_12_1=Exp[0]=1.5
	A_f_12_1=A_f[0]=0.67
	C_aiming_12_2=C_aiming[0]=0.7000000000000002
	Exp_12_2=Exp[0]=1.5
	A_f_12_2=A_f[0]=0.67
	
	# flowpath 11
	C_aiming_11_1=C_aiming[0]=0.7000000000000002
	Exp_11_1=Exp[0]=1.5
	A_f_11_1=A_f[0]=0.67
	C_aiming_11_2=C_aiming[0]=0.7000000000000002
	Exp_11_2=Exp[0]=1.5
	A_f_11_2=A_f[0]=0.67
	
	# flowpath 10
	C_aiming_10_1=C_aiming[0]=0.7000000000000002
	Exp_10_1=Exp[0]=1.5
	A_f_10_1=A_f[0]=0.67
	C_aiming_10_2=C_aiming[0]=0.7000000000000002
	Exp_10_2=Exp[0]=1.5
	A_f_10_2=A_f[0]=0.67
	
	# flowpath 9
	C_aiming_9_1=C_aiming[0]=0.7000000000000002
	Exp_9_1=Exp[0]=1.5
	A_f_9_1=A_f[0]=0.67
	C_aiming_9_2=C_aiming[0]=0.7000000000000002
	Exp_9_2=Exp[0]=1.5
	A_f_9_2=A_f[0]=0.67
		
	# flowpath 8
	C_aiming_8_1=C_aiming[4]=0.8000000000000003
	Exp_8_1=Exp[4]=1.5
	A_f_8_1=A_f[4]=0.33
	C_aiming_8_2=C_aiming[12]=0.9000000000000004
	Exp_8_2=Exp[12]=7.5000000000e-01
	A_f_8_2=A_f[12]=7.4500000000e-01
		
	# flowpath 7
	C_aiming_7_1=C_aiming[11]=0.8000000000000003
	Exp_7_1=Exp[11]=1.5
	A_f_7_1=A_f[11]=0.33
	C_aiming_7_2=C_aiming[3]=0.9000000000000004
	Exp_7_2=Exp[3]=7.5000000000e-01
	A_f_7_2=A_f[3]=7.4500000000e-01
		
	# flowpath 6
	C_aiming_6_1=C_aiming[5]=0.6500000000000001
	Exp_6_1=Exp[5]=1.5
	A_f_6_1=A_f[5]=0.33
	C_aiming_6_2=C_aiming[13]=1.0000000000000004
	Exp_6_2=Exp[13]=1.5
	A_f_6_2=A_f[13]=0.67
		
	# flowpath 5
	C_aiming_5_1=C_aiming[10]=0.6500000000000001
	Exp_5_1=Exp[10]=1.5
	A_f_5_1=A_f[10]=0.33
	C_aiming_5_2=C_aiming[2]=1.0000000000000004
	Exp_5_2=Exp[2]=1.5
	A_f_5_2=A_f[2]=0.67
		
	# flowpath 4
	C_aiming_4_1=C_aiming[6]=0.7500000000000002
	Exp_4_1=Exp[6]=1.5
	A_f_4_1=A_f[6]=0.33
	C_aiming_4_2=C_aiming[14]=0.8000000000000003
	Exp_4_2=Exp[14]=1.5
	A_f_4_2=A_f[14]=0.67
		
	# flowpath 3
	C_aiming_3_1=C_aiming[9]=0.7500000000000002
	Exp_3_1=Exp[9]=1.5
	A_f_3_1=A_f[9]=0.33
	C_aiming_3_2=C_aiming[1]=0.8000000000000003
	Exp_3_2=Exp[1]=1.5
	A_f_3_2=A_f[1]=0.67
		
	# flowpath 2
	C_aiming_2_1=C_aiming[7]=0.8000000000000003
	Exp_2_1=Exp[7]=1.5
	A_f_2_1=A_f[7]=0.33
	C_aiming_2_2=C_aiming[15]=0.7000000000000002
	Exp_2_2=Exp[15]=1.5
	A_f_2_2=A_f[15]=0.67
	
	# flowpath 1
	C_aiming_1_1=C_aiming[8]=0.8000000000000003
	Exp_1_1=Exp[8]=1.5
	A_f_1_1=A_f[8]=0.33
	C_aiming_1_2=C_aiming[0]=0.7000000000000002
	Exp_1_2=Exp[0]=1.5
	A_f_1_2=A_f[0]=0.67

	print C_aiming
	print Exp
	print A_f
	
	aiming_results,eff_interception,Strt=Model.aiming_loop(C_aiming,Exp,A_f)
	gx=aiming_results[2][15]
	print aiming_results[2]
	print eff_interception
	f=[gx]
	retval['fns'] = f
	print f
	return(retval)
