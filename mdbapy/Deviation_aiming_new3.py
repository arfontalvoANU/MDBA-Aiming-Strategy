'''
The aiming method is developed following the deviation-based multiple aiming.
Reference: Augsburger G. Thermo-economic optimisation of large solar tower power plants[R]. EPFL, 2013.
Basen on the method, this script improved the aiming by using different tube banks.
We call it localised deviation aiming strategy.
The input parameters include:
r_height: receiver height
r_diameter: receiver diameter
c_aiming: an aiming co-efficient array
csv: the field layout file
'''
import numpy as np
from sys import path

def return_metric(N_hst,A):
	Metric=np.zeros(N_hst) # the index metric
	if A>=0.5:
		N_hst_up=1
		N_hst_down=0
		Metric[0]=1
	else:
		N_hst_down=1
		N_hst_up=0
		Metric[0]=-1
	
	ratio=float(N_hst_up)/(N_hst_up+N_hst_down)
	#print 0,ratio,N_hst_up,N_hst_down,Metric[0]
	i=1
	while i<N_hst:
		ratio=float(N_hst_up)/(N_hst_up+N_hst_down)
		if ratio>A:
			Metric[i]=-1
			N_hst_down+=1
		elif ratio==A:
			Metric[i]=Metric[0]
			if Metric[i]==1:
				N_hst_up+=1
			else:
				N_hst_down+=1
		elif ratio<A:
			Metric[i]=1
			N_hst_up+=1
		#print i,ratio,N_hst_up,N_hst_down,Metric[i]
		i+=1
	return Metric

def aiming(folder,r_height,r_diameter,C_aiming,csv,tower_h,num_bundle,Exp,A_f,stand_by=True):
	#print C_aiming
	title=np.array(['x', 'y', 'z', 'foc', 'aim x', 'aim y', 'aim z', 'm', 'm', 'm', 'm', 'm', 'm', 'm'])
	r_radius=0.5*r_diameter
	hst_info=np.loadtxt(csv,delimiter=',', skiprows=2) 
	num_hst=int(hst_info.size/7)
	N_hst=num_hst
	num_bundle=int(num_bundle)
	Azimuth_boundary=np.zeros(num_bundle+1) # from -90 to 270
	for i in range(num_bundle+1):
		Azimuth_boundary[i]=-90.+360./num_bundle*i
	#print Azimuth_boundary
	
	Hst_info=[]  # the array to store all the hst_info according to the relevant tube bank
	for i in range(num_bundle):
		
		A=np.array([]) # a transition array
		for j in range(num_hst):
			# to calculate the azimuth angle
			x=hst_info[j,0]
			y=hst_info[j,1]
			if abs(x)<1e-6:
				azimuth=90.
			else:
				azimuth=np.arctan(y/x)/np.pi*180.				
			if x<0:
				azimuth+=180.
			# to compare with the boundaries
			if azimuth>=Azimuth_boundary[i] and azimuth<Azimuth_boundary[i+1]:
				A=np.append(A,hst_info[j,:])
		A=A.reshape(int(len(A)/7), 7)
		Hst_info.append(A)
		
	pos_and_aiming_new=np.array([])
	Hst_info_ranked=[]
	for i in range(num_bundle):
		hst_info=Hst_info[i]
		#print(len(hst_info))
		foc=hst_info[:,3]
		M=return_metric(N_hst,A_f[i])
		hst_info_ranked = hst_info[np.argsort(foc)[::1]]
		for j in range(int(len(hst_info))):
			lmax=np.max(hst_info_ranked[:,3])
			lmin=np.min(hst_info_ranked[:,3])
			li=hst_info_ranked[j,3]
			d=0.5*r_height*C_aiming[i]*((lmax-li)/(lmax-lmin))**Exp[i]
			hst_info_ranked[j,6]=tower_h+0.5*r_height+d*M[j]			
			hst_info_ranked[j,4]=hst_info_ranked[j,0]*r_radius/np.sqrt(hst_info_ranked[j,0]**2+hst_info_ranked[j,1]**2)
			hst_info_ranked[j,5]=hst_info_ranked[j,1]*r_radius/np.sqrt(hst_info_ranked[j,0]**2+hst_info_ranked[j,1]**2)
			hst_info_ranked[j,3]=np.sqrt((hst_info_ranked[j,0]-hst_info_ranked[j,4])**2+(hst_info_ranked[j,1]-hst_info_ranked[j,5])**2+(hst_info_ranked[j,2]-hst_info_ranked[j,6])**2)
		Hst_info_ranked.append(hst_info_ranked)
	
	# to move heliostats aiming outside to the stand-by point
	safety_aiming=[50.,50.,tower_h+r_height*0.5]
	Hst_stand=[]
	#if (not np.all(C_aiming<=1.)):
	for i in range(num_bundle):
		hst_stand=0
		J=[]
		for j in range(int(len(Hst_info_ranked[i]))):			
			if Hst_info_ranked[i][j,6]>(tower_h+r_height) or Hst_info_ranked[i][j,6]<tower_h:
				#print pos_and_aiming_new[i,6],i
				hst_stand+=1
				Hst_info_ranked[i][j,4]=safety_aiming[0]
				Hst_info_ranked[i][j,5]=safety_aiming[1]
				Hst_info_ranked[i][j,6]=safety_aiming[2]
				J.append(j)
		if stand_by!=True:
			Hst_info_ranked[i] = np.delete(Hst_info_ranked[i], J, axis=0)
		Hst_stand.append(hst_stand)
		pos_and_aiming_new=np.append(pos_and_aiming_new, Hst_info_ranked[i])
	#else:
		#for i in range(num_bundle):
			#pos_and_aiming_new=np.append(pos_and_aiming_new, Hst_info_ranked[i])
	
	pos_and_aiming_new=np.append(title,pos_and_aiming_new)
	pos_and_aiming_new=pos_and_aiming_new.reshape(int(len(pos_and_aiming_new)/7), 7)
	csv_new='%s/pos_and_aiming_new.csv' % folder # the output field file
	np.savetxt(csv_new, pos_and_aiming_new, fmt='%s', delimiter=',')
	return Hst_info_ranked,Hst_stand

if __name__=='__main__':
	r_height=18. # receiver height
	r_diameter=16.
	folder=path[0]
	csv='%s/1.15/pos_and_aiming_trimmed.csv' % folder # the input field file
	tower_h=175.
	num_bundle=16
	# C_aiming is an coefficienct array. Each coefficient corresponds to one tube bank.
	# The orientation is the same as the meshing of the cylinder: south-east-north-west
	C_aiming=np.zeros(num_bundle)
	C_aiming[:]=1.1
	A_f=np.zeros(num_bundle)
	A_f[:]=0.5 
	Exp=np.zeros(num_bundle)
	Exp[:]=1.4
	aiming(folder,r_height,r_diameter,C_aiming,csv,tower_h,num_bundle,Exp,A_f)
