import subprocess
import sys
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import math
import time
from sys import path

class SolsticeScene:
	def __init__(self,mainfolder,num_rays,dni,azimuth,zenith,att_factor,csv,tower_h,h_cyl,r_cyl,num_bundle,hst_w, hst_h, mirror_reflectivity=0.9, slope_error=0.0015, sunshape='buie', sunshape_param=0.02):

		'''
		Import and define the parameters

		Arguments:
		receiver - str, 'flat' or 'blade' or 'STL'
		folder - str, the directory of saving the case
		casename - str, define that name of the case, for saving files
		pmfile - str, the directory of the parameter input csv file
		annual - bool, do an annual performance simulation or not
		sunshape - str, 'buie' or 'pillbox'
		sunshape_param - float, the csr for buie sunshape, or the half angle for pillbox sunshape, note that the unit is deg for the pillbox angle 
		'''
		self.mainfolder=mainfolder
		self.folder=mainfolder
		self.casename='demo'
		
		# Solar Related parameters        
		self.azimuth=float(azimuth)
		self.zenith=float(zenith)
		self.num_rays=num_rays
		self.dni=float(dni)
		# heliostat related parameters
		self.mirror_reflectivity=mirror_reflectivity
		#self.slope=0.0014020
		self.slope=slope_error
		self.hst_dir=csv
		self.hst_w=hst_w
		self.hst_h=hst_h
		self.att_factor=att_factor
		self.tower_h=tower_h
		self.sunshape=sunshape
		self.sunshape_param=sunshape_param

		# receiver related parameters
		self.absorptivity=1. # absorptivity of the pipes
		self.rec_x=0.
		self.rec_y=0.
		self.rec_z=tower_h+h_cyl*0.5
		self.rec_slice=100
		self.h_cyl=h_cyl
		self.r_cyl=r_cyl
		self.num_bundle=num_bundle
		
	def gen_YAML(self):

		'''
		Generate YAML file 
		'''
		#TODO try pyYAML?

		# OPEN the input YAML file
		fd = open('%s/%s.yaml'%(self.folder, self.casename), 'w')

		# ----------------------------- CREATION OF THE SUN ------------------------
		#
		# CREATE The sun
		#hal_angle = 4.65e-3*180./np.pi
		#st_dev = 2.51e-3*180./np.pi
		#fd.write('- sun: {dni:  %s, pillbox: {half_angle: %s}}\n' % (self.dni,hal_angle))
		if self.sunshape=='buie':
			sunshape_par_n='csr'
		elif self.sunshape=='pillbox':
			sunshape_par_n='half_angle'

		fd.write('- sun: {dni:  %s, %s: {%s: %s}}\n' % (self.dni, self.sunshape, sunshape_par_n, self.sunshape_param))
		#fd.write('- sun: {dni:  %s, gaussian: {std_dev: %s}}\n' % (self.dni,st_dev))

		# ----------------------------- CREATION OF THE ATMOSPHERE ------------------
		#
		#print self.att_factor
		fd.write('- atmosphere: {extinction: %s}\n' % self.att_factor)

		# ------------------------------CREATION OF MATERIALS ----------------------
		#
		# CREATE an occultant material
		r_f = 0. # front
		r_b = 0. # and back reflectivity
		fd.write('- material: &%s\n' % 'material_black')
		fd.write('   front:\n')
		fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_f))    
		fd.write('   back:\n')
		fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b))
		fd.write('\n')
		#
		# CREATE a material for the target

		
		r_f = 1.-self.absorptivity # front
		r_b = 1. # and back reflectivity
		fd.write('- material: &%s\n' % 'material_rec')
		fd.write('   front:\n')
		fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_f))    
		fd.write('   back:\n')
		fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b))
		fd.write('\n')

		# CREATE a specular material
		slope_error = self.slope 
		r_b = 0. # and back reflectivity
		fd.write('- material: &%s\n' % 'material_mirror')
		fd.write('   front:\n')
		fd.write('     mirror: {reflectivity: %6.4f, slope_error: %15.8e }\n' % (self.mirror_reflectivity,float(slope_error) ) ) 
		fd.write('   back:\n')
		fd.write('     matte: {reflectivity: %6.4f }\n' % float(r_b))
		fd.write('\n')
		#
		# CREATE a material for the Large target to compute spillage
		fd.write('- material: &%s\n' % 'material_virtual')
		fd.write('   virtual:\n')
		fd.write('\n')


		#----------------------- CREATION OF GEOMETRIES -----------------------------
		#
		

		#    Receiver Geometry
		#

		# IMPORT CSV file for the heliostat positions
		hst_info=np.loadtxt(self.hst_dir,delimiter=',', skiprows=2)
		hst_x=hst_info[:,0]
		hst_y=hst_info[:,1]
		hst_z=hst_info[:,2]

		foc=hst_info[:,3]

		aim_x=hst_info[:,4]
		aim_y=hst_info[:,5]
		aim_z=hst_info[:,6]

		#foc = np.sqrt((hst_x-rec_x)**2+ (hst_y-rec_y)**2+(hst_z-rec_z)**2) # ideal focus

		h_hst = self.hst_h # heliostat height
		w_hst = self.hst_w # heliostat width
		slices = 4 # slices for the envelop circle
		pts_hst = [ [-w_hst*0.5, -h_hst*0.5], [-w_hst*0.5, h_hst*0.5], [w_hst*0.5, h_hst*0.5], [w_hst*0.5,-h_hst*0.5] ]

		# CREATE a reflective facet (mirror) used in the  heliostat template
		for i in range(0,len(foc)):
			name_hst_g = 'hst_g_'+str(i)
			fd.write('- geometry: &%s\n' % name_hst_g )
			fd.write('  - material: *%s\n' % 'material_mirror' )
		#    fd.write('    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) )
			fd.write('    parabol: \n')
			fd.write('      focal: %s\n' % foc[i]) 
			fd.write('      clip: \n')    
			fd.write('      - operation: AND \n')
			fd.write('        vertices: %s\n' % pts_hst)
			fd.write('      slices: %d\n' % slices )  

		# CREATE the pylon "pylon_g" geometry cylindrical shape
		h_pyl = 0.001 # pylon height
		r_pyl = 0.2 # pylon radius
		slices = 4 # slices for the envelop circle
		fd.write('- geometry: &%s\n' % 'pylon_g' )
		fd.write('  - material: *%s\n' % 'material_black' )
		fd.write('    transform: { translation: %s, rotation: %s }\n' % ([0, 0, -h_pyl*3], [0, 90, 0]) )
		fd.write('    cylinder: {height: %7.3f, radius: %7.3f, slices: %d }\n' % (h_pyl,r_pyl,slices) )
		#
			
		# -------------------------------------- CREATE THE TEMPLATES using the geometries

		# CREATE the heliostat templates

		for i in range(0,len(foc)):    
			name_hst_t = 'hst_t_'+str(i)
			fd.write('- template: &%s\n' % name_hst_t )
			name_hst_n = 'hst_'+ str(i)
			fd.write('    name: %s\n' % name_hst_n )
			fd.write('    primary: 0\n' )    
			fd.write('    geometry: *pylon_g\n')
			fd.write('    children: \n' )
			fd.write('    - name: pivot\n')
			fd.write('      zx_pivot: {target: {position: %s}} \n' % ([aim_x[i],aim_y[i],aim_z[i]]) )
			fd.write('      children: \n')
			fd.write('      - name: reflect_surface\n')
			fd.write('        primary: 1\n')
			fd.write('        transform: {rotation: [-90,0,0]} \n')    
			name_hst_g = 'hst_g_'+str(i)
			fd.write('        geometry: *%s\n' % name_hst_g )



        #
        #
        # -------------------------------------- CREATE THE ENTITIES using the geometries or the templates
        #
        ### CREATE entities from geometries: specifying if primary reflector
        #
		#   Tower Geometry 
		#- cylindrical shape
		fd.write("- entity: \n")
		fd.write("    name: tower_g\n")
		fd.write("    primary: 0\n")
		fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.tower_h*0.5], [0., 0, 0]))
		fd.write("    geometry:\n")
		fd.write("    - material: *material_rec\n")
		fd.write('      cylinder: \n')
		fd.write('        height: %s\n' % (1.))  
		fd.write('        radius: %s\n' % (self.r_cyl))
		fd.write('\n')
		#
        
		# the pipes
		
		fd.write("- entity: \n")
		fd.write("    name: cylinder\n")
		fd.write("    primary: 0\n")
		fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.rec_z], [0., 0, 0]))
		fd.write("    geometry:\n")
		fd.write("    - material: *material_rec\n")
		fd.write('      cylinder: \n')
		fd.write('        height: %s\n' % (self.h_cyl))  
		fd.write('        radius: %s\n' % (self.r_cyl))
		fd.write('        slices: %s\n' % self.num_bundle)
		fd.write('        stacks: %s\n' % 50)
		fd.write('\n')
		
		'''
		# build the central plate, which is not a round
		fd.write('\n- entity:\n')
		fd.write('    name: "round_plate"\n')
		fd.write('    primary: 0\n')
		fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.rec_z-0.5*self.h_cyl], [0., 180., 0]))
		fd.write('    geometry: \n')
		fd.write('    - material: *material_rec\n' )
		fd.write('      plane: \n')
		fd.write('        clip: \n')    
		fd.write('        - operation: AND \n')
		fd.write('          circle: {radius: %s, center: [0, 0]}\n' % self.r_cyl)
		fd.write('        slices: 10\n')
		fd.write('\n')
		'''
		
		# the virtual plate to calculate spillage
		width_vir=1000.
		pts = [[-width_vir, -width_vir], [-width_vir, width_vir], [width_vir, width_vir], [width_vir,-width_vir]]
		fd.write('\n- entity:\n')
		fd.write('    name: "virtual_target_e"\n')
		fd.write('    primary: 0\n')
		fd.write('    transform: { translation: %s, rotation: %s }\n' % ([self.rec_x, self.rec_y, self.rec_z-self.h_cyl*0.5-1.], [0, 0, 0]))
		fd.write('    geometry: \n')
		fd.write('    - material: *%s\n' % 'material_virtual' )
		fd.write('      plane: \n')
		fd.write('        clip: \n')
		
		fd.write('        - operation: AND \n')
		fd.write('          vertices: %s\n' % pts)
		fd.write('\n')
		

		#
		#
		### CREATE entities from templates
		#
		#
		for i in range(0,len(foc)):
			name_e ='H_'+str(i)
			name_hst_t = 'hst_t_'+str(i)
			fd.write('\n- entity:\n')
			fd.write('    name: %s\n' % name_e)
			fd.write('    transform: { translation: %s, rotation: %s }\n' % ([hst_x[i], hst_y[i], hst_z[i]], [0, 0, 0]) )
			fd.write('    children: [ *%s ]\n' % name_hst_t)

		# END OF THE YAML INPUT FILE
		fd.close() 


		# WRITE THE YAML FILE FOR THE RECEIVERS
		fd = open('%s/%s-rcv.yaml'%(self.folder,self.casename), 'w')

		#/* Receivers */
		fd.write("- {name: cylinder, side: FRONT, per_primitive: ABSORBED}\n")
		#fd.write("- {name: round_plate, side: FRONT, per_primitive: ABSORBED}\n")
		fd.write('- name: virtual_target_e \n' )
		fd.write('  side: %s \n' % 'BACK')
		fd.write('  per_primitive: %s \n' % 'INCOMING')
		
		
		# -------------------------------------- END OF THE YAML PARSER WRITTING
		fd.close()

	def runSOLSTICE(self, savefile, azi=0., zenith=0., view=False):
		'''
		run SOLSTICE with the corresponding sun position  
		and postprocessing the result   
		azi: from East to North
		zenith: 0 is the horizontal (in Solstice)
		view - if check it in paraview   
		'''
		azi=self.azimuth
		zenith=self.zenith

		os.system('solstice -D%s,%s -t 4 -v -n %s -R %s/%s-rcv.yaml -fo %s/simul %s/%s.yaml'%(azi,zenith, self.num_rays, self.folder, self.casename,savefile,self.folder,self.casename))


