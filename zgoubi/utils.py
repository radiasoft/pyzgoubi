#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from math import *
import numpy
import sys
from glob import glob
from constants import *
# use these to convert things to metres and tesla
m = 1
cm = 0.01
mm = 0.001
T = 1
kgauss = 0.1

# define things in metres, and tesla
# then use these when setting up elements to convert
# to a zgoubi unit
m_ = 1
cm_ = 100
mm_ = 1000
kgauss_ = 10
T_ = 1


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


def flatten(x):
	"""flatten(sequence) -> list

	Returns a single, flat list which contains all elements retrieved
	from the sequence and all recursively contained sub-sequences
	(iterables).

	http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

	Examples:
	>>> [1, 2, [3,4], (5,6)]
	[1, 2, [3, 4], (5, 6)]
	>>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
	[1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

	result = []
	for el in x:
	    #if isinstance(el, (list, tuple)):
	    if hasattr(el, "__iter__") and not isinstance(el, basestring):
		result.extend(flatten(el))
	    else:
		result.append(el)
	return result


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


def mom_to_ke(mom,mass):
	"""mom in eV/c
	gives results in eV
	"""
	ke = sqrt(mom**2 + mass**2) - mass

	return ke


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

		#clean up tmp directory
		line.clean()
		
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


def get_twiss_profiles(line, file_result,input_twiss_parameters = [0,0,0,0,0,0]):
	""" Calculates the twiss parameters at all points written to zgoubi.plt. 11 particle trajectories are used to calculate the
	transfer matrix, just as is done in zgoubi. The code mirrors that found in mat1.f, mksa.f. The twiss parameters are first
	calculated at the end of the cell using either input_twiss_parameters (format [beta_y, alpha_y, gamma_y, beta_z, alpha_z, gamma_z]) 
	or, if this is not supplied, it assumes the cell is periodic and uses get_twiss_parameters to find the boundary condition. 

	The twiss parameters are then mapped to all points in the magnets using the transfer matrix calculated at each point. The results are stored
	in list twiss_profiles where each row represents a point in the zgoubi.plt file and has format 

	[s_coord, label,  mu_y, beta_y, alpha_y, gamma_y, mu_z, beta_z, alpha_z, gamma_z]

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
#!-----------------------------------------------------------------------------------------
#!  PREPARE DATA FOR CALCULATION
#!-----------------------------------------------------------------------------------------
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
	label = [x[0].strip() for x in r.get_track('plt', ['element_label1'])]

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

#!-----------------------------------------------------------------------------------------
#!  TRANSFER MATRIX CALCULATION AT ALL POINTS in PLT file (as in zgoubi source file mat1.f)
#!-----------------------------------------------------------------------------------------

#pyreport latex
#$ The linear transfer matrix relates the starting coordinates to the final coordinates \newline
#$ \begin{pmatrix}
#$ R11 & R12 & R13 & R14 & R15 & R16 \\
#$ R21 & R22 & R23 & R24 & R25 & R26 \\
#$ R31 & R32 & R33 & R34 & R35 & R36 \\
#$ R41 & R42 & R43 & R44 & R45 & R46 \\
#$ R51 & R52 & R53 & R54 & R55 & R56 \\
#$ R61 & R62 & R63 & R64 & R65 & R66 \end{pmatrix}*
#$\begin{pmatrix}
#$ Y0 \\ T0 \\ Z0 \\ P0 \\ X0 \\ D0 \end{pmatrix}=
#$\begin{pmatrix}
#$ Y \\ T \\ Z \\ P \\ X \\ D \end{pmatrix}

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
		
#!----------------------------------
#! Adjust units to SI (as in mksa.f)
#!----------------------------------
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
	

#! Get inital twiss paramters. If no input_twiss_parameters supplied, assume cell is periodic and find results using get_twiss_parameters

	if input_twiss_parameters == [0,0,0,0,0,0]:
		twissparam = r.get_twiss_parameters()
		beta_y_0 = twissparam[0]
		alpha_y_0 = twissparam[1]
		gamma_y_0 = twissparam[2]
		beta_z_0 = twissparam[3]
		alpha_z_0 = twissparam[4]
		gamma_z_0 = twissparam[5]
	else:
		beta_y_0 = input_twiss_parameters[0]
		alpha_y_0 = input_twiss_parameters[1]
		gamma_y_0 = input_twiss_parameters[2]
		beta_z_0 = input_twiss_parameters[3]
		alpha_z_0 = input_twiss_parameters[4]
		gamma_z_0 = input_twiss_parameters[5]


#! Calculate twiss parameters at all points in plt file 
#!######################################################

#! The twiss parameters at the end of the periodic cell have been calculated by call to get_twiss_parameters, otherwise twiss parameters are input 
#!
#! - e.g. S.Y.Lee Accelerator physics equation 2.54
#! - The phase advance can be found by applying a Floquet transformation to the transfer matrix (S.Y.Lee eqn 2.65)

 
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

	#put twiss parameters together. Format [s_coord, mu_y, beta_y, alpha_y, gamma_y, mu_z, beta_z, alpha_z, gamma_z]. Units are SI
	twiss_profiles = [[s*cm for s in S_alltracks[0]],label_ref,mu_y_list,beta_y_list,alpha_y_list,gamma_y_list,mu_z_list,beta_z_list,alpha_z_list,gamma_z_list]

	return twiss_profiles


def fourier_tune(line,initial_YTZP,D_in,nfourierturns, coords = []):
	"""Calculate tune using FFT. nfourierturns determines the number of passes through the lattice.
       Can supply set of horizontal and vertical coordinates in coords = [ycoords,zcoords], otherwise
         routine will calculate coordinates
	"""
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


	#check line has a REBELOTE
	for e in line.element_list:
		if ("REBELOTE" in str(type(e)).split("'")[1]):
			reb=e
			break
	else:
		raise ValueError, "Line has no REBELOTE element"

	if coords == []:
		objet.clear()	# remove existing particles
		#start at closed orbit with a small amplitude added to vertical component to get tune readings
		objet.add(Y=initial_YTZP[0], T=initial_YTZP[1], Z=initial_YTZP[2]+1e-5, P=initial_YTZP[3], LET='A', D=D_in)

		reb.set(NPASS=nfourierturns-1)	
		r = line.run(xterm = False)
		YZ = r.get_track('fai', ['Y','Z'])
		ycoords = numpy.transpose(YZ)[0]
		zcoords = numpy.transpose(YZ)[1]
	else:
		#use input coords
		ycoords = coords[0]
		zcoords = coords[1]



	#perform FFT
	import pylab
	yfft = pylab.fft(ycoords)
	zfft = pylab.fft(zcoords)
	#extract amplitudes
	yampfreq = numpy.abs(yfft)
	zampfreq = numpy.abs(zfft)

	import operator
	ysortfreqamp = sorted(enumerate(yampfreq), key=operator.itemgetter(1))
	ysortfreqamp.reverse()
	zsortfreqamp = sorted(enumerate(zampfreq), key=operator.itemgetter(1))
	zsortfreqamp.reverse()

	#drop peak at index 0
	if ysortfreqamp[0][0] == 0:
		del ysortfreqamp[0]
	if zsortfreqamp[0][0] == 0:
		del zsortfreqamp[0]

	ypeaksloc = sorted([ysortfreqamp[0][0],ysortfreqamp[1][0]])
	zpeaksloc = sorted([zsortfreqamp[0][0],zsortfreqamp[1][0]])

	yfouriertune = ypeaksloc[0]/len(ycoords)
	zfouriertune = zpeaksloc[0]/len(zcoords)

	#plot tunes if desired
	plot_fourier = False
	if(plot_fourier):
		pylab.plot(yampfreq, 'k-')
		pylab.plot(zampfreq, 'b-')
		pylab.savefig('fspectrum')

	return yfouriertune, zfouriertune

def scan_dynamic_aperture(line, emit_list, closedorb_YTZP, npass, D_mom, beta_gamma_input = 1, ellipse_coords = 1, plot_data = False):
        """ Check a list of emittances (emit_list, units Pi m rad) to see if tracking succeeds. Can be used to establish the dynamic aperture. If the elements in emit_list are increasing then will stop tracking once it finds the lowest emittance where a particle is lost. On the other hand, if the elements in emit_list are decreasing then will stop tracking once it reaches the first point where tracking succeeds without loss. 

			Required input:
				emit_list - List of emittances to check
				closedorb_YTZP - Closed orbit coordinates as returned by find_closed_orbit
				npass - Number of passes through lattice
				D_mom - Momentum factor p/pref
			Optional input
				beta_gamma_input - If this is supplied then the emittances supplied in emit_list are assumed to be normalised and a conversion to geometrical emittance is made. Otherwise the emittances are assumed to be geometrical.
				ellipse_coords - If greater than 1 will test uniformly distributed set of coordinates around around phase space 			ellipse. Otherwise will take single point where phase space cuts y (or z) axis.
						Can also specify ellipse coords = [n,m] -- distribute n points around ellipse but only test point m of these. 
				plot_data - If True, creates phase space plots at all emittances scanned in both transverse planes.
	
		If a particle is lost, returns index_lost in emit_list where the loss occurs. Otherwise index_lost remains 0.
        """

	import numpy
	import os
	import zgoubi.core as zg

	#create objet so initial phase settings can be set
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


	#create reb so that number of passes can be set
	for e in line.element_list:
			if ("REBELOTE" in str(type(e)).split("'")[1]):
					reb = e
					break
	else:
		raise ValueError, "Line has no REBELOTE element"

	rigidity = objet.BORO

	reb.set(NPASS=npass-1)

	coord_pick = None
	if type(ellipse_coords) == list:
	    coord_pick = ellipse_coords[1]
	    ellipse_coords = ellipse_coords[0]

	#calculate optical parameters on closed orbit. This is required to convert emittances into a coordinate.
	objet5 = zg.OBJET5()
	line.replace(objet, objet5)
	objet5.set(BORO=rigidity)
	objet5.set(PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
	objet5.set(YR=closedorb_YTZP[0],TR=closedorb_YTZP[1],ZR=closedorb_YTZP[2],PR=closedorb_YTZP[3], DR=D_mom)
	matrix= zg.MATRIX(IORD=1,IFOC=11)
	line.replace(reb,matrix)
	r = line.run(xterm = False)
	twissparam = r.get_twiss_parameters()
	#alphayz = [twissparam[1],twissparam[4]]
	betayz = [twissparam[0],twissparam[3]]
	gammayz = [twissparam[2],twissparam[5]]

	#revert to objet2 mode with rebelote
	line.replace(objet5, objet)
	objet.clear()	# remove existing particles
	objet.add(Y=closedorb_YTZP[0], T=closedorb_YTZP[1], Z=closedorb_YTZP[2], P=closedorb_YTZP[3], LET='A', D=D_mom)
	line.replace(matrix,reb)


	index_lost = None
	coord_index_lost = None
	from operator import add
	if plot_data:
	    YTZP_list = []

	fourier_tune_emit = []

	#if emit_list is a list of decreasing emittances, looking for the first point where particle is not lost
	reverse_search = False
	if emit_list[0] > emit_list[-1]:
		print "decreasing list "
		reverse_search = True
		

	for emit in emit_list:
		print "check emit ",emit
		if coord_pick == None:
			#obtain coordinates on phase space ellipses. Use beta, gamma values found at closed orbit
			coords_YTZP_ini = emittance_to_coords(emit, emit, gammayz, betayz, beta_gamma_input, ncoords = ellipse_coords)
		else:
			#select one point on phase space ellipse. Use beta, gamma values found at closed orbit
			coords_YTZP_ini = emittance_to_coords(emit, emit, gammayz, betayz, beta_gamma_input, ncoords = ellipse_coords)[coord_pick]

		try:
			l = len(coords_YTZP_ini[0])
			coords_YTZP_ini = [map(add,closedorb_YTZP,coords) for coords in coords_YTZP_ini]
		except TypeError:
			coords_YTZP_ini = [map(add,closedorb_YTZP,coords_YTZP_ini)] 
	
		for coord_index, current_YTZP in enumerate(coords_YTZP_ini):

			print "ellipse coord index ",coord_index
			objet.clear()	# remove existing particles
			objet.add(Y=current_YTZP[0], T=current_YTZP[1], Z=current_YTZP[2], P=current_YTZP[3], LET='A', D=D_mom)
				
			#run Zgoubi
			r = line.run(xterm= False)
    
			rebelote_completed = r.test_rebelote()

			lost = not rebelote_completed
    
			if(lost):
				print "tracking failed at emit = ",emit
				coord_index_lost = coord_index
				#break
				#return [index_lost, coord_index],fourier_tune_emit

			if plot_data and not lost:
				YTZP_list.append(r.get_track('fai', ['Y','T','Z','P']))
				coords_in = [flatten(r.get_track('fai', ['Y'])),flatten(r.get_track('fai', ['Z']))]
				fourier_tune_result = fourier_tune(line,[],1,1, coords = coords_in)
				fourier_tune_emit.append(fourier_tune_result)

		#condition to end scan
		if reverse_search and not lost:
			index_lost = emit_list.index(emit)+1
			break
		elif not reverse_search and lost:
			index_lost = emit_list.index(emit)
			break			


	else:
		print "Successful tracking in all test cases"

	if plot_data:
		#obtain coordinates on phase space ellipse using closed orbit twiss parameters
		coords_YTZP_full = emittance_to_coords(emit, emit, gammayz, betayz, beta_gamma_input, ncoords = 100)
		coords_YTZP_full = [map(add,closedorb_YTZP,coords) for coords in coords_YTZP_full]
		coords_YTZP_full.append(coords_YTZP_full[0])

		#add many coordinates to draw phase space ellipse
		Y_data = [numpy.transpose(coords_YTZP_full)[0]]
		T_data = [numpy.transpose(coords_YTZP_full)[1]]
		Z_data = [numpy.transpose(coords_YTZP_full)[2]]
		P_data = [numpy.transpose(coords_YTZP_full)[3]]


		#add points on phase space ellipse actually used to initialise scan above
		Y_data.append([numpy.transpose(coords_YTZP_ini)[0]])
		T_data.append([numpy.transpose(coords_YTZP_ini)[1]])
		Z_data.append([numpy.transpose(coords_YTZP_ini)[2]])
		P_data.append([numpy.transpose(coords_YTZP_ini)[3]])
    
		for index in range(len(YTZP_list)):
			Y_data.append(numpy.transpose(YTZP_list[index])[0])
			T_data.append(numpy.transpose(YTZP_list[index])[1])
			Z_data.append(numpy.transpose(YTZP_list[index])[2])
			P_data.append(numpy.transpose(YTZP_list[index])[3])


		plot_data_xy_multi(Y_data,T_data,'yt_phasespace', labels=["Horizontal phase space","y [cm]","y' [mrad]"],style = ['k-','ro','b+','r+','g+','m+','y+'],xlim=[-32.5,-29.5],ylim=[-40,40])
		plot_data_xy_multi(Z_data,P_data,'zp_phasespace', labels=["Vertical phase space","z [cm]","z' [mrad]"],style = ['k-','ro','b+','r+','g+','m+','y+'],xlim=[-2.0,2.0],ylim=[-30,30])

		#y turn-by-turn, including initial point
		y_turnbyturn = numpy.transpose(YTZP_list[index])[0]
		y_turnbyturn = [y for y in y_turnbyturn]
		y_turnbyturn.reverse()
		y_turnbyturn.append(coords_YTZP_ini[0][0])
		y_turnbyturn.reverse()
		y_turnbyturn = [y-closedorb_YTZP[0] for y in y_turnbyturn]
		yrange = [i+1 for i in range(len(y_turnbyturn))]
		#z turn-by-turn, including initial point
		z_turnbyturn = [z for z in numpy.transpose(YTZP_list[index])[2]]
		z_turnbyturn.reverse()
		z_turnbyturn.append(coords_YTZP_ini[0][2])
		z_turnbyturn.reverse()
		#z_turnbyturn = [z-closedorb_YTZP[2] for z in z_turnbyturn]
		zrange = [i+1 for i in range(len(z_turnbyturn))]

		#plot_data_xy_multi(yrange,y_turnbyturn,'y_turnbyturn', labels=["Y turn-by-turn","turn","y [cm]"],style = ['k+'],xlim = [0,len(yrange)+1])
		#plot_data_xy_multi(zrange,z_turnbyturn,'z_turnbyturn', labels=["Z turn-by-turn","turn","z [cm]"],style = ['k+'],xlim = [0,len(zrange)+1])
		plot_data_xy_multi(y_turnbyturn,z_turnbyturn,'yz_space', labels=["YZ coords","y [cm]","z [cm]"],style = ['k+'])
	    

	return [index_lost, coord_index],fourier_tune_emit



def emittance_to_coords(emit_horizontal, emit_vertical, gammayz, betayz, beta_gamma_input = 1, ncoords = 1):
	"""Given some emittance in horizonal and vertical space

	If ncoords >= 1 return points where phase space ellipse crosses the y,y' and z,z' axis.

	If ncoords > 1, will instead give a distribution of points (y,t) around the phase space ellipse uniform in angle theta where
	    y = a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi)
	    t = a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi)
	where (a,b) are major and minor radii and phi is the tilt of the major radius w.r.t horizontal axis. A (z,p) distribution is 
	similarly calculated.

	Can use these points, returned in coords_YTZP to start tracking at the desired emittance (assuming that the optical functions don't change with amplitude).

	Emittances in both the horizontal and vertical planes may be supplied. Twiss parameters beta and gamma in 
	both places may be determined calling get_twiss_parameters beforehand i.e.
			twissparam = r.get_twiss_parameters()
			betayz = [twissparam[0],twissparam[3]]
			gammayz = [twissparam[2],twissparam[5]]"""
	import numpy
	
	coords_YTZP = []

	if ncoords <= 1:
	    #crosses y or z axis at sqrt(emit/gamma). Convert to cm.
	    y_axis=cm_*sqrt(emit_horizontal/(gammayz[0]*beta_gamma_input))       
	    z_axis=cm_*sqrt(emit_vertical/(gammayz[1]*beta_gamma_input))
	    YTZP=[y_axis,0,z_axis,0]
	    coords_YTZP.append(YTZP)
    
	    #crosses y' or z' axis at sqrt(emit/beta) and convert to mrad (same as scaling by mm)
	    t_axis =mm_*sqrt(emit_horizontal/(betayz[0]*beta_gamma_input))       
	    p_axis =mm_*sqrt(emit_vertical/(betayz[1]*beta_gamma_input))
	    YTZP=[0,t_axis,0,p_axis]
	    #coords_YTZP.append(YTZP)

	else:
	    start_angle = 0.0
	    end_angle = 2*pi
	    angle_list = [start_angle+i*(end_angle-start_angle)/(ncoords) for i in range(ncoords)]
	    emityz = [emit_horizontal/beta_gamma_input, emit_vertical/beta_gamma_input]
	    ydat = []
	    tdat = []
	    zdat = []
	    pdat = []
	    for index in range(2):
		#calculate twiss parameter alpha
		alpha = (abs(betayz[index]*gammayz[index]-1))**0.5
		#calculate major and minor radii in each plane
		h = 0.5*(betayz[index] + gammayz[index])
		major_radius = ((emityz[index]/2)**0.5)*( (h+1)**0.5 + (h-1)**0.5 )
		minor_radius = ((emityz[index]/2)**0.5)*( (h+1)**0.5 - (h-1)**0.5 )
		phi = 0.5*numpy.arctan(-2*alpha/(betayz[index]-gammayz[index]))

		#decide which axis is the major one
		if betayz[index] >= gammayz[index]:
		    horiz_radius = major_radius
		    vert_radius = minor_radius
		else:
		    horiz_radius = minor_radius
		    vert_radius = major_radius		

		for theta in angle_list:
		    y_z = horiz_radius*cos(theta)*cos(phi) - vert_radius*sin(theta)*sin(phi)
		    t_p = horiz_radius*cos(theta)*sin(phi) + vert_radius*sin(theta)*cos(phi)
		    if index == 0:
			ydat.append(cm_*y_z)
			tdat.append(mm_*t_p)
		    else:
			zdat.append(cm_*y_z)
			pdat.append(mm_*t_p)

	    #plot_data_xy_multi(ydat,tdat,'ytdat', style = ['k+','b+','r+','g+'])
	    #plot_data_xy_multi(zdat,pdat,'zpdat', style = ['k+','b+','r+','g+'])

	    #put coords_YTZP together
	    for index in range(ncoords):
		coords_YTZP.append([ydat[index],tdat[index],zdat[index],pdat[index]])

	return coords_YTZP


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


def scaling_to_dipole(k, r0, b0_in, d_r0=0, scale_factor=1.0, terms=4):
	"""
		A scaling FFAG has field given by B/B0 = (r/r0)^k where k is the scaling factor.
		The Taylor expansion about r0 yields
		
		B(r)/B0 = (k!/(n!(k-n)!)) * (r-r0)^n/r0^n for each multipole n

		In DIPOLES the magnetic field (ignoring fringe field factor) is given by
	
		B(r)/B0 = bcoef[n] (r-r0)^n/r0^n for each multipole n

		It follows that setting bcoef[n] = (k!/(n!(k-n)!)) will create a DIPOLE magnet whose field approximates the scaling law.
		Routine returns bcoef up to number of terms specified. 

		If r0 is scaled by scale_factor, the coefficients are scaled accordingly. 
		If r0 is shifted by d_r0, the routine will update B0 (b0_in) according to the scaling law"""


	#calculate k*(k-1)*(k-2)*.../n!
	bcoef = [k]
	kfac = k
	for i in range(terms-1):
		kfac = kfac*(k-i-1)/(i+2)
		bcoef.append(kfac)
		
	#scale coefficients to be consistent with radius scaling
	for index in range(len(bcoef)):
		bcoef[index] = bcoef[index]*scale_factor**(index+1)

	#adjust B0 if r0 is shifted by d_r0
	b0 = b0_in*(1 + d_r0/r0)**k

	return bcoef, b0


def get_enclosing_circle(ellipse_data):
	""" Find smallest circle that encloses a set of ellipses centred on the midplane. The ellipses are defined by their horizontal and vertical radii and by their centre along the horizontal axis (a,b,c). The algorithm optimised both the centre of the enclosing circle and its radius. 

	Ellipse algorithm by J. Scott Berg (BNL)
	
	Reference - J. Scott Berg. 'Finding the Circular Magnet Aperture which Encloses and Arbitrary number of midplace-centered Beam Ellipses', Proc. EPAC04, 5-9 July, Lucerne Congress Centre, Lucerne, Switzerland

	""" 

	import ellipse

	bc = ellipse.BestCircle()

	for edata in ellipse_data:
		a = edata[0]
		b = edata[1]
		c = edata[2]
		bc.append((a,b,c))

	(centre,radius) = bc.get_circle()
	print 'centre of enclosing circle =', centre, ", with radius =", radius

	return centre, radius


def misalign_element(line, element_indices, mean, sigma, sigma_cutoff, misalign_dist = [], seed = 123456):
	""" Misalign the set of elements identified by element_indices. The misalignment is Gaussian
	with sigma, sigma_cutoff and mean among the input parameters. If sigma = 0.0, the misalignment is constant 
	and given by mean

	Alternatively, misalignments specified by list misalign_dist

	Returns the misailgnment distribution (misalign_dist), unit meters

	element_indices can be determined using line.find_elements(element_name)

	mean and sigma should be input in meters """

	import random
	import zgoubi.core as zg
	random.seed(seed)

	#generator function
	def surround(l):
		update_index = 0
		for e in l:
			yield e+update_index
			update_index = update_index+2
	
	#indice of entrance changref
	changref_indices = list(surround(element_indices))

	#generate misalignment distribution (microns)
	if misalign_dist == []:
		while len(misalign_dist) < len(element_indices):
			if sigma > 0.0:
				rand = random.gauss(mean,sigma)
				if abs(rand) < sigma_cutoff*sigma:
					misalign_dist.append(rand)
			else:
				misalign_dist.append(mean)

	#insert CHANGREF with misalignments given by misalign_dist.
	for index, pos in enumerate(changref_indices):
		misalign_elem = zg.CHANGREF('misalign',YCE= misalign_dist[index]*cm_)
		misalign_elem_rev = zg.CHANGREF('misalign',YCE= -misalign_dist[index]*cm_)
		line.insert(pos, misalign_elem)
		line.insert(pos+2, misalign_elem_rev)

	return misalign_dist


def plot_data_xy(data, filename, labels=["","",""], style='b-', xlim = [0,0], ylim = [0,0]):
	import pylab
	data_a = numpy.array(data)
	pylab.hold(False)
	pylab.plot(data_a[:,0], data_a[:,1], style)
	pylab.title(labels[0])
	if xlim != [0,0]:
		pylab.xlim( (xlim[0] ,xlim[1]) )
	if ylim != [0,0]:
		pylab.ylim( (ylim[0] ,ylim[1]) )
	pylab.xlabel(labels[1])
	pylab.ylabel(labels[2])
	pylab.savefig(filename)


def plot_data_xy_multi(data_x_list, data_y_list, filename, labels=["","",""], style='', legend = ' ', legend_location = 'best',xlim = [0,0], ylim = [0,0]):
	import pylab
	""" Plots multiple sets of data where the X and Y coordinates are each specified in a list of lists. Should also
	    work if a single set of X, Y data is specified or if one X is supplied with multiple Y data points (as long 
	    the dimensions of Y equals that of X in all cases). 

		Required input -
		data_x_list - Set of X data. May be a single list or a list of lists each with the same dimension. 
		data_y_list - Set of Y data. May be a single list or a list of lists each with the same dimension. Can have a single list of X data and multiple lists of data_y_list.
		filename - Plot saved to filename.png

		Optional parameters -
		labels = ["Plot Title","X axis label","Y axis label"]
		style - The plot style can be supplied as a list or as a single value. If no style is specified,  will cycle over a predefined style list.
		legend - Legends can be supplied as a list of strings or as a single string. 
		legend_location - Default of location of legend box is 'best', otherwise can select 'upper right' or numerical code for position. See matplotlib documentation. 
		xlim = [lower_value, upper_value] - Limit horizontal axis
		ylim = [lower_value, upper_value] - Limit vertical axis  """
                                                                  

	#check if data is a single list or a list of lists
	single_x_data = False
	single_y_data = False

	try:
		l = len(data_x_list[0])
	except TypeError:
		single_x_data = True

	try:
		l = len(data_y_list[0])
	except TypeError:
		single_y_data = True

	pylab.hold(True)

	if type(style) != list:
	    style = [style]
	if style == ['']:
	    style = ['k-','b-','r-','g-','m-','y-']

	if type(legend) != list:
	    legend = [legend]


	if single_y_data:
	    pylab.plot(data_x_list, data_y_list, style[0])
	else:
	    for index, data_y in enumerate(data_y_list):
			if single_x_data:
				pylab.plot(data_x_list, data_y, style[index%len(style)], label=legend[index%len(legend)])
			else:
				pylab.plot(data_x_list[index], data_y, style[index%len(style)], label= legend[index%len(legend)])

	pylab.title(labels[0])
	if xlim != [0,0]:
		pylab.xlim( (xlim[0] ,xlim[1]) )
	if ylim != [0,0]:
		pylab.ylim( (ylim[0] ,ylim[1]) )
	pylab.xlabel(labels[1])
	pylab.ylabel(labels[2])

	print "legend ",legend
	if legend != [' ']:
		pylab.legend(loc=legend_location)

	pylab.savefig(filename)
	pylab.cla()

