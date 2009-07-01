#!/usr/bin/env python

from __future__ import division
from math import *
import numpy
import sys
from glob import glob
from zgoubi_constants import *
# use these to convert things to metres and tesla
m = 1
cm = 0.01
mm = 0.001
T = 1

# define things in metres, and tesla
# then use these when setting up elements to convert
# to a zgoubi unit
m_ = 1
cm_ = 100
mm_ = 1000
kgauss_ = 10


def coords_grid(min=[0,0], max=[1,1], step=[1,1], type=float):
	"""make a list of coordinate aranged in a grid in a generalised space
	eg 
	coords_grid(min=[0,0],max=[1,1],step=[10,10])
	will give you 100 points in 2 dimensions between 0 and 1

	"""
	assert(len(min)==len(max)==len(step))
	#assert(len(min)==2)
	dimention = len(min)

	stride=[(max[x]-min[x])/step[x] for x in xrange(dimention)]

	indexes = numpy.zeros(step)

	coords = []

	for ind, dummy in numpy.ndenumerate(indexes):
		#print ind
		this_coord = []
		for x in xrange(len(ind)):
			this_coord.append(stride[x] * ind[x] + min[x])
		coords.append(this_coord)

	return numpy.array(coords, dtype=type)


def search_pattern(step, range, start=0):
	current = start
	while True:
		yield current
		if (current > 0):
			current = 0-current
		elif (current < 0):
			current = (0 - current) + step
		else: # zero case
			current += step
		if (abs(current)>range):
			break


def get_cmd_param(key, default=None):
	"get things formated like key=value from the commandline. returns the requested value"
	#print sys.argv
	params = {}
	for item in sys.argv:
		if '=' in item:
			k,e,v = item.partition('=')
			if (k == str(key)):
				return v
	if default != None:
		return default
	print key, "not given on the command line"
	print "please run with "+ str(key) +"=x as an argument"
	raise ValueError


def get_cmd_param_bool(key, default=None):
	"""get things formated like key=value from the commandline.
	returns true for: yes, on, true, 1
	returns false for: no, off, false, 0
	case insensitve
	"""
	#print sys.argv
	params = {}
	for item in sys.argv:
		if '=' in item:
			k,e,v = item.partition('=')
			if (k == str(key)):
				if v.lower() in ['yes', 'on', 'true', '1', 'oui']:
					return True
				if v.lower() in ['no', 'off', 'false', '0', 'non']:
					return False
				
	if default != None:
		if type(default) == type(True):
			return default
		else:
			print "default must be a bool (True or False)"
			raise ValueError
	print key, "not given on the command line"
	print "please run with "+ str(key) +"=x as an argument"
	raise ValueError

def readArray(filename, skipchar = '#', dtype=float):
	"""
	# PURPOSE: read an array from a file. Skip empty lines or lines
	#          starting with the comment delimiter (defaulted to '#').
	#
	# OUTPUT: a float numpy array
	#
	# EXAMPLE: >>> data = readArray("/home/albert/einstein.dat")
	#          >>> x = data[:,0]        # first column
	#          >>> y = data[:,1]        # second column
	from http://www.scipy.org/Cookbook/InputOutput

	"""

	myfile = open(filename, "r")
	contents = myfile.readlines()
	myfile.close()

	data = []
	for line in contents:
		stripped_line = line.lstrip()
		if (len(stripped_line) != 0):
			if (stripped_line[0] != skipchar):
				items = stripped_line.split()
				data.append(map(float, items))

	a = numpy.array(data, dtype=dtype)
	(Nrow,Ncol) = a.shape
	if ((Nrow == 1) or (Ncol == 1)): a = ravel(a)

	#print >> sys.stderr, "Read",Nrow,"rows in",Ncol,"columns in file", filename
	return(a)


def find_centre(ellipse):
	"find centre by using mean of coords, this assumes that point are evenly distributed around ellipse"
	try:
		cx = ellipse[:,0].mean()
		cy = ellipse[:,1].mean()
		return (cx, cy)
	except TypeError:
		pass
			

	sx = 0
	sy = 0
	for point in ellipse:
			x, y = point
			sx += x
			sy += y

	cx = sx / len(ellipse)
	cy = sy / len(ellipse)
	return (cx, cy)


def calc_area_simple(ellipse, centre=(0,0)):
	"can't handle noise. finds closest and furthest point from centre, assumes they are a and b"
	if len(ellipse) < 2:
		raise NoTrackError, "Not enough data points to measure area"
	distances = []
	for point in ellipse:
		x, y = point
		x -= centre[0]
		y -= centre[1]
		dist = sqrt(x**2 + y**2)
		distances.append(dist)

	a = min(distances)
	b = max(distances)

	return a * b * pi

def ke_to_rigidity(ke, mass):
	"""ke in eV, mass in eV/c^2
	gives result in kGauss cm, which is handy for OBJET
	"""
	te = ke + mass
	mom = sqrt(te**2 - mass**2) # in eV/c

	rigidity = mom / SPEED_OF_LIGHT * 1000 # to give kG cm

	return rigidity

def mom_to_rigidity(mom):
	"""mom in eV/c
	gives result in kGauss cm, which is handy for OBJET
	"""

	rigidity = mom / SPEED_OF_LIGHT * 1000 # to give kG cm

	return rigidity

def ke_to_relativistic_beta_gamma(ke, mass):
	"""ke in eV, mass in eV/c^2
	"""
	te = ke + mass
	mom = sqrt(te**2 - mass**2) # in eV/c

	beta_gamma_rel = mom / mass

	return beta_gamma_rel

def show_file(file_path,mode):
	"""shows a file. mode can be:
	cat - dumps to screen
	less - scrollable view
	win - scollable new win (less in an xterm)
	
	"""
	if mode == "cat":
		fh = open(file_path)
		for line in fh:
			print line,

	elif mode == "less":
		from os import system
		command = "less %s"% file_path
		system(command)
	
	elif mode == "win":
		from os import system
		command = 'xterm -e "less %s"'% file_path
		system(command)


def find_closed_orbit(line, init_YTZP=[0,0,0,0], max_iterations=100, tol = 1e-8, D=1):
	#check line has an objet2
	for e in line.element_list:
		if ("OBJET2" in str(type(e)).split("'")[1]):
			objet = e
			break
	else:
		raise ValueError, "Line has no OBJET2 element"
	
	for e in line.element_list:
		if ("FAISCNL" in str(type(e)).split("'")[1]):
			break
	else:
		raise ValueError, "Line has no FAISCNL element"

	current_YTZP = init_YTZP
	areas = []
	coords = []
	tracks = []
	close_orbit_found = False
	for iteration in xrange(max_iterations):
		coords.append(current_YTZP)
		objet.clear()	# remove existing particles
		objet.add(Y=current_YTZP[0], T=current_YTZP[1], Z=current_YTZP[2], P=current_YTZP[3], LET='A', D=D)

		r = line.run(xterm=False)
		#r = Results(line)

		track = r.get_track('fai', ['Y','T','Z','P'])
		if not r.run_success():
			print "No stable orbit"
			return None
		track.insert(0, current_YTZP)
		track_a = numpy.array(track)
		tracks.append(track_a)
		
		centre_h = find_centre(track_a[:,0:2])
		area_h = calc_area_simple(track_a[:,0:2], centre=centre_h)
		centre_v = find_centre(track_a[:,2:4])
		area_v = calc_area_simple(track_a[:,2:4], centre=centre_v)


		areas.append([area_h , area_v])
		current_YTZP = [centre_h[0], centre_h[1], centre_v[0], centre_v[1]]
		if area_h < tol and area_v < tol:
			close_orbit_found = True
			close_orbit = current_YTZP
			break


	
	for x in xrange(len(areas)):
		print "it:",x, "\tYTZP", coords[x], "\tarea", areas[x]


	if close_orbit_found:
		print "found closed orbit"
		print "Y=%s, T=%s, Z=%s, P=%s"%tuple(current_YTZP)
		return current_YTZP
	else:
		print "Iterations did not converge, no closed orbit found"
		return None


def get_twiss_profiles(line,file_result):
	""" Calculates the twiss parameters at all points written to zgoubi.plt. 11 particle trajectories are used to calculate the
	transfer matrix, just as is done in zgoubi. The code mirrors that found in mat1.f, mksa.f. The twiss parameters are first
	calculated at the end of the cell using get_twiss_parameters and then mapped to the magnets using the transfer matrix at
	each point. The results are stored in list twiss_profiles where each row represents a point in the zgoubi.plt file and
	has format [s_coord, label,  mu_y, beta_y, alpha_y, gamma_y, mu_z, beta_z, alpha_z, gamma_z]

	Requires an OBJET type 5, and a MATRIX element.

        Note - This calculation uses trajectories as measured in the local coordinate system of the magnet."""


	has_object5 = False
       	has_matrix = False
       	for e in line.element_list:
       		t = str(type(e)).split("'")[1].rpartition(".")[2]
       		if t == 'OBJET5':
       			has_object5 = True
       		if t == 'MATRIX':
       			has_matrix = True
       	if not (has_object5 and has_matrix):
       		raise BadLineError, "beamline need to have an OBJET with kobj=5 (OBJET5), and a MATRIX elementi to get tune"

	#run Zgoubi
	r = line.run(xterm = False)
	#-----------------------------------------------------------------------------------------
	#  PREPARE DATA FOR CALCULATION
	#-----------------------------------------------------------------------------------------
	plt_track = r.get_track('plt', ['LET','D0','Y0','T0','Z0','P0','X0','D','Y','T','Z','P','S','X'])
	transpose_plt_track = map(list,zip(*plt_track))
	track_tag = transpose_plt_track[0]
	D0 = transpose_plt_track[1]
	Y0 = transpose_plt_track[2]
	T0 = transpose_plt_track[3]
	Z0 = transpose_plt_track[4]			
	P0 = transpose_plt_track[5]
	X0 = transpose_plt_track[6]
	D = transpose_plt_track[7]
	Y = transpose_plt_track[8]
	T = transpose_plt_track[9]
	Z = transpose_plt_track[10]
	P = transpose_plt_track[11]	
	S = transpose_plt_track[12]
	X = transpose_plt_track[13]

	#read labels
	labelraw = flatten(r.get_track('plt', ['element_label1']))
        #strip out trailing blanks in label string
	label = []
	for lab in labelraw:
		label.append(lab.rstrip())

	#sort out individual tracks
	alphabet = "A,B,C,D,E,F,G,H,I,J,K,L,M,N,P,Q,R,S,T,U,V,W,X,Y,Z"
	alphabet = alphabet.split(",")

	D0_alltracks = []
	Y_alltracks = []
	Y0_alltracks = []
	T_alltracks = []
	T0_alltracks = []
	Z_alltracks = []
	Z0_alltracks = []
	P_alltracks = []
	P0_alltracks = []		
	S_alltracks = []
	label_ref = []
	Y_track = []
	T_track = []
	Z_track = []
	P_track = []
	S_track = []
	#first add track coords for reference track
	ref_indices = find_indices(track_tag,'O')
	D0_alltracks.append(D0[ref_indices[0]])
	Y0_alltracks.append(Y0[ref_indices[0]])
	T0_alltracks.append(T0[ref_indices[0]])
	Z0_alltracks.append(Z0[ref_indices[0]])
	P0_alltracks.append(P0[ref_indices[0]])
	for i in ref_indices:
		Y_track.append(Y[i])
		T_track.append(T[i])
		Z_track.append(Y[i])
		P_track.append(T[i])
		S_track.append(S[i])
		label_ref.append(label[i])
	Y_alltracks.append(Y_track)
	T_alltracks.append(T_track)
	Z_alltracks.append(Z_track)
	P_alltracks.append(P_track)
	S_alltracks.append(S_track)
	#add other tracks
	for track_tag_name in alphabet:
		Y_track = []
		T_track = []
		Z_track = []
		P_track = []
		S_track = []
		indices = find_indices(track_tag,track_tag_name)
		if len(indices) == 0:
			break
		else:
			D0_alltracks.append(D0[indices[0]])
			Y0_alltracks.append(Y0[indices[0]])
			T0_alltracks.append(T0[indices[0]])
			Z0_alltracks.append(Z0[indices[0]])
			P0_alltracks.append(P0[indices[0]])
			for i in indices:
				Y_track.append(Y[i])
				T_track.append(T[i])
				Z_track.append(Z[i])
				P_track.append(P[i])
				S_track.append(S[i])
			Y_alltracks.append(Y_track)
			T_alltracks.append(T_track)
			Z_alltracks.append(Z_track)
			P_alltracks.append(P_track)
			S_alltracks.append(S_track)	
       				

	#11 coordinate in Y0_alltracks,T0_alltracks etc correspond to 11 starting conditions required by MATRIX

	#-----------------------------------------------------------------------------------------
	#  TRANSFER MATRIX CALCULATION AT ALL POINTS in PLT file (as in mat1.f)
	#-----------------------------------------------------------------------------------------

	#initialise transfer matrix
	#transfer_matrix = zeros((6, 6))

	#momentum ratio DP
	# FORTRAN (mat1.f)   DP = ( FO(1,I10) - FO(1,I11) ) / (.5D0*( FO(1,I10) + FO(1,I11) ) )
	DP = ((1+D0_alltracks[9])-(1+D0_alltracks[10]))/(0.5*((1+D0_alltracks[9])+(1+D0_alltracks[10])))

        #Transfer matrix elements R(J-1,6), J in range (2,6)
	# FORTRAN (mat1.f)   R(J-1,6)  = (F(J,I10) - F(J,I11)) /DP
	R16_list = [x/DP for x in map(numpy.subtract,Y_alltracks[9],Y_alltracks[10])]
	R26_list = [x/DP for x in map(numpy.subtract,T_alltracks[9],T_alltracks[10])]
	R36_list = [x/DP for x in map(numpy.subtract,Z_alltracks[9],Z_alltracks[10])]
	R46_list = [x/DP for x in map(numpy.subtract,P_alltracks[9],P_alltracks[10])]
	R56_list = [x/DP for x in map(numpy.subtract,S_alltracks[9],S_alltracks[10])]
	#-----------------------------------------------------------------------------------------

	#assume it1=1
	it1 = 1
	for i in range(1,5):
		i2 = 2*i +it1-1
		i3 = i2 + 1
		#adjust indices since python lists are numbered from zero
		i2 = i2 - 1
		i3 = i3 - 1
		if i == 1:
			UO = Y0_alltracks[i2]-Y0_alltracks[i3]
			R11_list =  [x/UO for x in map(numpy.subtract,Y_alltracks[i2],Y_alltracks[i3])]
			R21_list =  [x/UO for x in map(numpy.subtract,T_alltracks[i2],T_alltracks[i3])]
			R31_list =  [x/UO for x in map(numpy.subtract,Z_alltracks[i2],Z_alltracks[i3])]
			R41_list =  [x/UO for x in map(numpy.subtract,P_alltracks[i2],P_alltracks[i3])]
			R51_list =  [x/UO for x in map(numpy.subtract,S_alltracks[i2],S_alltracks[i3])]
		elif i == 2:
			UO = T0_alltracks[i2]-T0_alltracks[i3]
			R12_list =  [x/UO for x in map(numpy.subtract,Y_alltracks[i2],Y_alltracks[i3])]
			R22_list =  [x/UO for x in map(numpy.subtract,T_alltracks[i2],T_alltracks[i3])]
			R32_list =  [x/UO for x in map(numpy.subtract,Z_alltracks[i2],Z_alltracks[i3])]
			R42_list =  [x/UO for x in map(numpy.subtract,P_alltracks[i2],P_alltracks[i3])]
			R52_list =  [x/UO for x in map(numpy.subtract,S_alltracks[i2],S_alltracks[i3])]
		elif i == 3:
			UO = Z0_alltracks[i2]-Z0_alltracks[i3]
			R13_list =  [x/UO for x in map(numpy.subtract,Y_alltracks[i2],Y_alltracks[i3])]
			R23_list =  [x/UO for x in map(numpy.subtract,T_alltracks[i2],T_alltracks[i3])]
			R33_list =  [x/UO for x in map(numpy.subtract,Z_alltracks[i2],Z_alltracks[i3])]
			R43_list =  [x/UO for x in map(numpy.subtract,P_alltracks[i2],P_alltracks[i3])]
			R53_list =  [x/UO for x in map(numpy.subtract,S_alltracks[i2],S_alltracks[i3])]
		elif i == 4:
			UO = P0_alltracks[i2]-P0_alltracks[i3]
			R14_list =  [x/UO for x in map(numpy.subtract,Y_alltracks[i2],Y_alltracks[i3])]
			R24_list =  [x/UO for x in map(numpy.subtract,T_alltracks[i2],T_alltracks[i3])]
			R34_list =  [x/UO for x in map(numpy.subtract,Z_alltracks[i2],Z_alltracks[i3])]
			R44_list =  [x/UO for x in map(numpy.subtract,P_alltracks[i2],P_alltracks[i3])]
			R54_list =  [x/UO for x in map(numpy.subtract,S_alltracks[i2],S_alltracks[i3])]
		
	#----------------------------------
	# Adjust units to SI (as in mksa.f)
	#----------------------------------
	unit_list = [1e-2,1e-3,1e-2,1e-3,1e-2,1,1e-6]
	R21_list = [x*unit_list[1]/unit_list[0] for x in R21_list]
	R31_list = [x*unit_list[2]/unit_list[0] for x in R31_list]
	R41_list = [x*unit_list[3]/unit_list[0] for x in R41_list]
	R51_list = [x*unit_list[4]/unit_list[0] for x in R51_list]

	R12_list = [x*unit_list[0]/unit_list[1] for x in R12_list]
	R32_list = [x*unit_list[2]/unit_list[1] for x in R32_list]
	R42_list = [x*unit_list[3]/unit_list[1] for x in R42_list]
	R52_list = [x*unit_list[4]/unit_list[1] for x in R52_list]

	R13_list = [x*unit_list[0]/unit_list[2] for x in R13_list]
	R23_list = [x*unit_list[1]/unit_list[2] for x in R23_list]
	R43_list = [x*unit_list[3]/unit_list[2] for x in R43_list]
	R53_list = [x*unit_list[4]/unit_list[2] for x in R53_list]

	R14_list = [x*unit_list[0]/unit_list[3] for x in R14_list]
	R24_list = [x*unit_list[1]/unit_list[3] for x in R24_list]
	R34_list = [x*unit_list[2]/unit_list[3] for x in R34_list]
	R54_list = [x*unit_list[4]/unit_list[3] for x in R54_list]

	R16_list = [x*unit_list[0]/unit_list[5] for x in R16_list]
	R26_list = [x*unit_list[1]/unit_list[5] for x in R26_list]
	R36_list = [x*unit_list[2]/unit_list[5] for x in R36_list]
	R46_list = [x*unit_list[3]/unit_list[5] for x in R46_list]
	R56_list = [x*unit_list[4]/unit_list[5] for x in R56_list]
	
        #------------------------------------------------------
        # Get inital twiss paramters
	#------------------------------------------------------
	twissparam = r.get_twiss_parameters()
	beta_y_0 = twissparam[0]
	alpha_y_0 = twissparam[1]
        gamma_y_0 = twissparam[2]
	beta_z_0 = twissparam[3]
	alpha_z_0 = twissparam[4]
	gamma_z_0 = twissparam[5]

        #----------------------------------------------------------------------------------------------------------------
        # Calculate twiss parameters at all points in plt file 
	#  - The twiss parameters at the end of the periodic cell have been calculated by call to get_twiss_parameters 
	#  - Twiss parameters at other points in the cell may calculated knowing these values and the transfer matrix.
	#  - e.g. S.Y.Lee Accelerator physics equation 2.54
	#  - The phase advance can be found by applying a Floquet transformation to the transfer matrix (S.Y.Lee eqn 2.65)
        #-----------------------------------------------------------------------------------------------------------------
 
	fresults= open(file_result, 'w')
              
	mu_y_list = []
	beta_y_list= []
	alpha_y_list = []
	gamma_y_list = []
	mu_z_list = []
	beta_z_list= []
	alpha_z_list = []
	gamma_z_list = []
	sign_sine_y_old = 1.0
	sign_sine_z_old = 1.0
	n_pi_y = 0
	n_pi_z = 0
	for i in range(len(R11_list)):
		# Horizontal plane
		#-----------------
		beta_y = (R11_list[i]**2)*beta_y_0 - 2.0*R11_list[i]*R12_list[i]*alpha_y_0 + (R12_list[i]**2)*gamma_y_0
		beta_y_list.append(beta_y)
		# Horizontal phase advance calculation
		# To account for range of acos=(0,Pi), check sign of sin(angle) to know when to change from angle to (2*Pi-angle) etc.	
		sine_angle = R12_list[i]/sqrt(beta_y*beta_y_0)
		if abs(sine_angle) > 1: 
			mu_y_list.append(mu_y_list[i-1])
		else:
			sign_sine_y = numpy.sign(asin(sine_angle))
			if sign_sine_y - sign_sine_y_old == -2:
				n_pi_y = n_pi_y + 1
			sign_sine_y_old = sign_sine_y
			if sign_sine_y >= 0:
				mu_y_list.append(n_pi_y*2*pi+acos(sqrt(beta_y_0/beta_y)*R11_list[i]-alpha_y_0*R12_list[i]/(beta_y*beta_y_0)**0.5))
			else:
				mu_y_list.append(n_pi_y*2*pi-acos(sqrt(beta_y_0/beta_y)*R11_list[i]-alpha_y_0*R12_list[i]/(beta_y*beta_y_0)**0.5))
		alpha_y_list.append(-R11_list[i]*R21_list[i]*beta_y_0 + (R11_list[i]*R22_list[i]+R12_list[i]*R21_list[i])*alpha_y_0 \
			- R12_list[i]*R22_list[i]*gamma_y_0)
		gamma_y_list.append((R21_list[i]**2)*beta_y_0-2*R21_list[i]*R22_list[i]*alpha_y_0+(R22_list[i]**2)*gamma_y_0)
		# Vertical plane
		#---------------
		beta_z = (R33_list[i]**2)*beta_z_0 - 2.0*R33_list[i]*R34_list[i]*alpha_z_0 + (R34_list[i]**2)*gamma_z_0
		beta_z_list.append(beta_z)
		# Vertical phase advance calculation
		sine_angle = R34_list[i]/sqrt(beta_z*beta_z_0)
		if abs(sine_angle) > 1:
			mu_z_list.append(mu_z_list[i-1])
		else:
			sign_sine_z = numpy.sign(asin(sine_angle))
			if sign_sine_z - sign_sine_z_old == -2:
				n_pi_z = n_pi_z + 1
			sign_sine_z_old = sign_sine_z
			if sign_sine_z >= 0:
				mu_z_list.append(n_pi_z*2*pi+acos(sqrt(beta_z_0/beta_z)*R33_list[i]-alpha_z_0*R34_list[i]/(beta_z*beta_z_0)**0.5))
			else: 
				mu_z_list.append(n_pi_z*2*pi-acos(sqrt(beta_z_0/beta_z)*R33_list[i]-alpha_z_0*R34_list[i]/(beta_z*beta_z_0)**0.5))
		alpha_z_list.append(-R33_list[i]*R43_list[i]*beta_z_0 + (R33_list[i]*R44_list[i]+R34_list[i]*R43_list[i])*alpha_z_0 \
			- R34_list[i]*R44_list[i]*gamma_z_0)
		gamma_z_list.append((R43_list[i]**2)*beta_z_0-2*R43_list[i]*R44_list[i]*alpha_z_0+(R44_list[i]**2)*gamma_z_0)
		print >>fresults, '%2f %2s %2f %2f %2f %2f %2f %2f %2f %2f' % (S_alltracks[0][i], label_ref[i], mu_y_list[i],beta_y,alpha_y_list[i],\
			gamma_y_list[i],mu_z_list[i],beta_z,alpha_z_list[i],gamma_z_list[i])

	#put twiss parameters together. Each row has format [s_coord, mu_y,beta_y, alpha_y, gamma_y, mu_z,beta_z, alpha_z, gamma_z]. Units are SI
	twiss_profiles = numpy.transpose([[s*cm for s in S_alltracks[0]],label_ref,mu_y_list,beta_y_list,alpha_y_list,gamma_y_list,mu_z_list,beta_z_list,alpha_z_list,gamma_z_list])

        return twiss_profiles



def find_indices(list,list_element):
	""" Find all occurences of list_element in list. Return the indices.
	"""
	indices = []
	i = -1
	try:
		while 1:
			i = list.index(list_element, i+1)
			indices.append(i)
	except ValueError:
		pass
	
	return indices


def plot_data_xy(data, filename, labels=["","",""]):
	import pylab
	data_a = numpy.array(data)
	print data_a
	pylab.hold(False)
	pylab.plot(data_a[:,0], data_a[:,1])
	pylab.title(labels[0])
	pylab.xlabel(labels[1])
	pylab.ylabel(labels[2])
	pylab.savefig(filename)

