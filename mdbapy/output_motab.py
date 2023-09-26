import numpy as np
from datetime import datetime

def output_motab(table,savedir=None, title=None):
	'''
	output the .motab table fiel
	'''
	f=open(savedir, 'w')
	f.write('#1\n')

	if table.ndim==2:
		# size of the lookup table
		m=np.shape(table)[0]-2
		n=np.shape(table)[1]-2
		f.write('double optics(%s, %s)\n'%(m,n))

		hour_angle=table[2, 3:]
		declination=table[3:,2]

		for i in range(m):
			if i ==0:
				row_i=np.append(0, hour_angle)
			else:
				row_i=np.append(declination[i-1], table[2+i, 3:])

			#content=np.array2string(row_i, formatter={'float_kind':lambda x: "%.2f" % row_i})
			#content=np.array2string(row_i, precision=2, separator=' ', suppress_small=True)
			#f.write(content+'\n')
			f.write(" ".join(map(str, row_i)))
			f.write("\n")

	else:
		# 3D table, include the breakdown of the total energy
		a=len(table)
		m=np.shape(table[0])[0]-2
		n=np.shape(table[0])[1]-2

		hour_angle=table[0][2, 3:]
		declination=table[0][3:,2]

		for t in range(a):
			f.write('double %s(%s, %s)\n'%(title[t], m,n))

			for i in range(m):
				if i ==0:
					row_i=np.append(0, hour_angle)
				else:
					row_i=np.append(declination[i-1], table[t][2+i, 3:])
				f.write(" ".join(map(str, row_i)))
				f.write("\n")	
						
			f.write("\n")
			f.write("\n")

	f.close()


def output_matadata_motab(table, field_type, aiming, n_helios, A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,coefs_T, coefs, eff_abs, eff_emi,SM,savedir=None, details_en=None):
	"""Output the .motab file to work with the SolarTherm program

	``Arguments``
		* table (numpy array): the oelt that returned by design_crs.CRS.field_design_annual, listed by declination and hour angles
		* field_type (str): polar or surrounding
		* aiming (str): 'single' or 'isp' or others to specify the aiming strategy
		* n_helios (int): total number of heliostats
		* A_helio (float): area of each heliostat (m2)
		* eff_design (float): the optical efficiency at design point
		* d_receiver: receiver diameter (m)
		* h_receiver: receiver height (m)
		* H_tower (float), tower height (m)
		* eff_rec_design (float): the receiver efficiency at design point
		* savedir (str): the directory to save this .motab file
		* details_en (dict): the key of the dict is the name of the breakdown of energy (loss), the value of the dict is a numpy array that contains the value of the energy (loss) with the corresponding declination and hour angles
	
	``Return``
		write the table(s) to the .motab file
	"""
	f=open(savedir, 'w')
	f.write('#1\n')
	f.write('#Comments: Field type: %s, Aiming Strategy: %s, Date:%s\n'%(field_type, aiming, datetime.now()))
	f.write("#METALABELS,n_helios,A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,CoT1,CoT2,CoT3,CoT4,CoT'1,CoT'2,CoT'3,CoT'4,Coh1,Coh2,Coh3,Coh4,Coh5,eff_abs,eff_emi,SM\n")
	f.write('##METAUNITS,integer,m2,real,m,m,m,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real\n')
	f.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(int(n_helios), A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,
	  coefs_T[0],coefs_T[1][0],coefs_T[1][1],coefs_T[1][2],coefs_T[2],coefs_T[3][0],coefs_T[3][1],coefs_T[3][2],coefs[0],coefs[1],coefs[2],coefs[3],coefs[4],eff_abs,eff_emi,SM))

	# size of the lookup table  
	m=np.shape(table)[0]
	n=np.shape(table)[1]
	f.write('double optics(%s, %s)\n'%(m,n))

	hour_angle=table[0, 1:]
	declination=table[1:,0]

	for i in range(m):
		if i ==0:
			row_i=np.append(0, hour_angle)
		else:
			row_i=np.append(declination[i-1], table[i, 1:])
		f.write(" ".join(map(str, row_i)))
		f.write("\n")

	if details_en!=None:
		for key in details_en:
			breakdown=key
			table=details_en[key]
			# size of the lookup table  
			m=np.shape(table)[0]
			n=np.shape(table)[1]
			f.write('double %s(%s, %s)\n'%(key,m,n))
			for i in range(m):
				if i ==0:
					row_i=np.append(0, hour_angle)
				else:
					row_i=np.append(declination[i-1], table[i, 1:])
				f.write(" ".join(map(str, row_i)))
				f.write("\n")			
	f.close()

if __name__=='__main__':
	from sys import path
	table=np.loadtxt('%s/F_optic.csv'%path[0],delimiter=',')
	coefs_T=[945.7112573259491, np.array([ 0.02720568, -0.00172737,  0.07126733]), 953.7130902079241, np.array([ 0.02170311, -0.00196636,  0.08407119])]
	coefs=[7.61828573e-04,-3.54208032e-02,5.93470995e-01,-9.37379885e-01,9.26793247e+00]
	eff_abs=0.987174393114
	eff_emi=0.940767053132


	output_matadata_motab(table, field_type='surrounding', aiming='MDBA', n_helios=6766, A_helio=144.375, eff_design=0.65414, 
	  d_receiver=16., h_receiver=24., H_tower=175., eff_rec_design=0.87645, coefs_T=coefs_T, coefs=coefs, eff_abs=eff_abs, eff_emi= eff_emi,
	  savedir='%s/example.motab'%path[0], details_en=None)


def read_motab(filename):

	with open(filename) as f:
		content=f.read().splitlines()
	f.close()
	res=content[4].split(',')

	n_helios=float(res[1])
	A_helio=float(res[2])
	eff_des=float(res[3])
	Q_in_rcv=float(res[-2])
	A_land=float(res[-1])

	oelt=np.array([])
	solar_hour=np.array([])
	declination=np.array([])
	t=content[6].split(' ')

	for v in t[1:]:
		solar_hour=np.append(solar_hour, float(v))

	for t in content[7:]:
		v=t.split(' ')
		declination=np.append(declination, float(v[0]))
		oelt=np.append(oelt, np.array(v[1:],dtype=float))

	oelt=oelt.reshape(len(declination), len(solar_hour))

	return n_helios, A_helio, eff_des, Q_in_rcv, A_land, solar_hour, declination, oelt
	


