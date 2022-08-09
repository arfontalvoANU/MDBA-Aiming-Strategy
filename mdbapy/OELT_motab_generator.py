import numpy as np
import os
from sys import path
class to_motab:
	def __init__(self,path,filename,C):
		# Loading csv file
		self.data = np.loadtxt(filename, delimiter=',')
		# Extracting declination and hour angles
		nhra = self.data.shape[1]
		ndec = self.data.shape[0]
		dec = np.linspace(-23.45,23.45,ndec)
		hra = np.linspace(-180.0,180.0,nhra)
		# Creating motab text file
		table = open('%s/case_H230_Sodium_290_565.motab'%path,'w+')
		# Fixed
		n_helios = 6764
		A_helio = 12.20*12.20
		eff_design = 0.646894510992
		d_receiver = 16.0
		h_receiver = 24.0
		H_tower = 175
		eff_rec_design = 0.915313131417
		# Writting headers
		table.write('#1\n')
		table.write('#Comments\n')
		table.write("#METALABELS,n_helios,A_helio, eff_design, d_receiver, h_receiver, H_tower, eff_rec_design,CoT1,CoT2,CoT3,CoT4,CoT'1,CoT'2,CoT'3,CoT'4,Coh1,Coh2,Coh3,Coh4,Coh5,eff_abs,eff_emi\n")
		table.write('#METAUNITS,integer,m2,real,m,m,m,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real,real\n')
		table.write('#METADATA,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'%(
			int(n_helios),
			A_helio,
			eff_design,
			d_receiver,
			h_receiver,
			H_tower,
			eff_rec_design,
			C[0],C[1],C[2],
			C[3],C[4],C[5],
			C[6],C[7],C[8],
			C[9],C[10],C[11],
			C[12],C[13],C[14]
			)
			)
		table.write('double optics(%s,%s)\n'%(ndec+1,nhra+1))
		table.write('0.0,')
		for j in range(nhra):
			if j==nhra-1:
				table.write('%s\n'%hra[j])
			else:
				table.write('%s,'% hra[j])

		# writting values
		for i in range(ndec):
			table.write('%s,'% dec[i])
			for j in range(nhra):
				if j==nhra-1:
					table.write('%s\n'%self.data[i,j])
				else:
					table.write('%s,'%self.data[i,j])
		table.close()

if __name__=='__main__':
	fpath = '%s/example/case_H230_Sodium_290_565'%path[0]
	data=np.loadtxt('%s/coefs.txt'%fpath, delimiter=',', skiprows=1)
	to_motab(fpath,'%s/OELT.csv'%fpath,data)
