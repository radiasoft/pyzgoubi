#!/usr/bin/env python
# -*- coding: utf-8 -*-

"Various useful functions and utilities"

from __future__ import division
from math import *
import numpy
import sys
import os
import warnings
from glob import glob
from zgoubi.common import *
from zgoubi.constants import *
from zgoubi.exceptions import *
from zgoubi.rel_conv import *
import itertools
from zgoubi.core import zlog, dep_warn
from StringIO import StringIO

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


def coords_grid(min_c=None, max_c=None, step=None, dtype=float):
	"""make a list of coordinate aranged in a grid in a generalised space eg::

		coords_grid(min_c=[0,0],max_c=[1,1],step=[10,10])
	
	will give you 100 points in 2 dimensions between 0 and 1
	"""

	if min_c == None:
		min_c = [0, 0]
	if max_c == None:
		max_c = [1, 1]
	if step == None:
		step = [1, 1]
	
	assert(len(min_c)==len(max_c)==len(step))
	#assert(len(min)==2)
	dimention = len(min_c)

	stride = [(max_c[x]-min_c[x])/step[x] for x in xrange(dimention)]

	indexes = numpy.zeros(step)

	coords = []

	for ind, dummy in numpy.ndenumerate(indexes):
		#print ind
		this_coord = []
		for x in xrange(len(ind)):
			this_coord.append(stride[x] * ind[x] + min[x])
		coords.append(this_coord)

	return numpy.array(coords, dtype=dtype)


def search_pattern(step, range, start=0):
	# this has been in the code since pre-sourceforge, i dont think it is used, and i dont think it works right if start!=0
	warnings.warn("utils.search_pattern()"+dep_warn, DeprecationWarning)
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
	for item in sys.argv:
		if '=' in item:
			k, e, v = item.partition('=')
			if (k == str(key)):
				return v
	if default != None:
		return default
	message = "'" + str(key) + "' not given on the command line\nplease run with "+ str(key) +"=x as an argument"
	zlog.error(message)
	raise ValueError


def get_cmd_param_bool(key, default=None):
	"""get things formated like key=value from the commandline.
	returns true for: yes, on, true, 1
	returns false for: no, off, false, 0
	case insensitve
	"""
	#print sys.argv
	for item in sys.argv:
		if '=' in item:
			k, e, v = item.partition('=')
			if (k == str(key)):
				if v.lower() in ['yes', 'on', 'true', '1', 'oui']:
					return True
				if v.lower() in ['no', 'off', 'false', '0', 'non']:
					return False
				zlog.error("Value for "+str(key)+" must be a bool (True or False)")
				raise ValueError
				
	if default != None:
		if type(default) == type(True):
			return default
		else:
			zlog.error("default must be a bool (True or False)")
			raise ValueError
	message = "'" + str(key) + "' not given on the command line\nplease run with "+ str(key) +"=x as an argument"
	zlog.error(message)
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

	Note: please switch to using numpy functions: numpy.savetxt(), numpy.loadtxt() or numpy.genfromtxt()

	"""

	warnings.warn("utils.readArray()"+dep_warn, DeprecationWarning)

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
	(Nrow, Ncol) = a.shape
	if ((Nrow == 1) or (Ncol == 1)): a = numpy.ravel(a)

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
		cx = ellipse[:, 0].mean()
		cy = ellipse[:, 1].mean()
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


def calc_area_simple(ellipse, centre=(0, 0)):
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


def mom_to_ke(mom, mass):
	"""mom in eV/c
	gives results in eV
	"""
	ke = sqrt(mom**2 + mass**2) - mass

	return ke


def ke_to_relativistic_beta(ke, mass):
	"""ke in eV, mass in eV/c^2
	"""
	te = ke + mass

	mom = sqrt(te**2 - mass**2) # in eV/c

	beta_rel = mom/te

	return beta_rel


def ke_to_relativistic_beta_gamma(ke, mass):
	"""ke in eV, mass in eV/c^2
	"""
	te = ke + mass
	mom = sqrt(te**2 - mass**2) # in eV/c

	beta_gamma_rel = mom / mass

	return beta_gamma_rel


def show_file(file_path, mode):
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
		command = "less %s"% file_path
		os.system(command)
	
	elif mode == "win":
		command = 'xterm -e "less %s"'% file_path
		os.system(command)

def find_closed_orbit_range(line, init_YTZP=None, max_iterations=100, fai_label = None, tol = 1e-6, D=1, range_YTZP=None, count_YTZP=None, record_fname=None):
	"""Same as find_closed_orbit, but if init_YTZP is unstable, will generate a bunch or particles, with a spread of range_YTZP, and if any of those are stable will do a closed orbit search with that particle::

		init_YTZP=[5,0,0,0], range_YTZP=[10,0,0,0], count_YTZP[50,0,0,0]

	will vary the Y coordinate from -5 to 15 cm in 50 steps::

		range_YTZP=[10,5,0,0], count_YTZP=[10,10,0,0]

	would create 100 particles in a grid.
	
	"""
	zlog.debug("enter function")
	if range_YTZP == None:
		range_YTZP = [10, 10, 10, 10]
	if init_YTZP == None:
		init_YTZP = [0, 0, 0, 0]
	if count_YTZP == None:
		count_YTZP = [10, 10, 10, 10]
	
	if record_fname:
		record_fh = open_file_or_name(record_fname, "w", mkdir=True)
	else:
		record_fh = None
	
	#first check center
	if record_fname:
		record_fh.write("#center\n")
	result = find_closed_orbit(line=line, init_YTZP=init_YTZP, max_iterations=max_iterations, fai_label=fai_label, tol=tol, D=D, record_fname=record_fh)
	if result != None:
		zlog.debug("Found a closed orbit without needing range")
		return result

	# if the init_YTZP is not stable, then send a bunch, and see if any of those are not lost
	
	for e in line.element_list:
		if ("OBJET2" in str(type(e)).split("'")[1]):
			objet = e
			break
	else:
		raise ValueError, "Line has no OBJET2 element"
	
	for e in line.element_list:
		if ("FAISCNL" or "FAISTORE" in str(type(e)).split("'")[1]):
			break
	else:
		raise ValueError, "Line has no FAISCNL element"

	current_YTZP = numpy.array(init_YTZP)
	objet.clear()	# remove existing particles
	
	# generate bunch
	ranges = []
	for x in xrange(4):
		if range_YTZP[x] == 0 or count_YTZP[0] <= 1:
			ranges.append([0])
		else:
			ranges.append(numpy.linspace(-range_YTZP[x], range_YTZP[x], count_YTZP[x]))

	n_coords = 0
	for coord in itertools.product(*ranges):
		n_coords += 1
		objet.add(Y=current_YTZP[0]+coord[0], T=current_YTZP[1]+coord[1], Z=current_YTZP[2]+coord[2], P=current_YTZP[3]+coord[3], LET='A', D=D)

	zlog.debug("Search for a stable orbit with %s particles"%n_coords)
	r = line.run(xterm=False)

	if not r.run_success():
		zlog.debug("No stable orbit within range: "+str(range_YTZP))
		return None
	else:
		# should have a track, so can find a stable particle
		track = r.get_all('fai')
		del r
		final_lap_n = track['PASS'].max()
		# select one particle which survived (IEX is positive), and made it to last lap
		surviving_particles = track[ numpy.logical_and(track['PASS'] == final_lap_n , track['IEX']>0)  ]
		zlog.debug("%d particles survived"%len(surviving_particles))

		#measure lengths and sort
		surviving_particles2 = [ [ sqrt((p['Y']-p['Y0'])**2 + (p['T']-p['T0'])**2 + (p['Z']-p['Z0'])**2 + (p['P']-p['P0'])**2),p] for p in surviving_particles ]
		surviving_particles2.sort()

		if record_fname:
			record_fh.write("#survivors %s\n"%len(surviving_particles))
			for length, surviving_particle in surviving_particles2:
				surviving_coords = (surviving_particle['Y0'], surviving_particle['T0'], surviving_particle['Z0'], surviving_particle['P0'], surviving_particle['Y'], surviving_particle['T'], surviving_particle['Z'], surviving_particle['P'])
				record_fh.write("%s %s %s %s %s %s %s %s\n"%surviving_coords)


		for length, surviving_particle in surviving_particles2[:10]:
			surviving_init_coord = [surviving_particle['Y0'], surviving_particle['T0'], surviving_particle['Z0'], surviving_particle['P0']]
		
			# use stable particle to find closed orbit
			result = find_closed_orbit(line=line, init_YTZP=surviving_init_coord, max_iterations=max_iterations, fai_label=fai_label, tol=tol, D=D, record_fname=record_fh)
			if result != None:
				return result
		zlog.warning("Despite finding surviving particles, none were stable")	
		return None


def find_closed_orbit(line, init_YTZP=None, max_iterations=100, fai_label = None, tol = 1e-6, D=1, record_fname=None):
	"""Find a closed orbit for the line. can optionally give a list of initial coordinates, init_YTZP, eg:
	find_closed_orbit(line, init_YTZP=[1.2,2.1,0,0])
	otherwise [0,0,0,0] are used.

	The line is expected to have a OBJET2 and either a FAISCNL at the end or, using FAISTORE, a MARKER at the end of the cell identified by fai_label. The latter option allows the zgoubi.fai file to contain coordinates throughout the lattice while using just those points with fai_label in this function.

	It is recommend to have a REBELOTE, with several laps. The area of the phase space ellipse is approximated from the coordinates from the FAISCNL (or MARKER with fai_label), and the center is used for the new coordinates. Once the relative variation between iterations is less that the tolerance, the function returns the closed orbit coordinates. If a coordinate is close to zero (less than the tolerance) then it is compared absolutely instead of relatively.
	
	record_fname is used to record search details to a file that can be used with plot_find_closed_orbit().
	"""
	zlog.debug("enter function")
	if init_YTZP == None:
		init_YTZP = [0, 0, 0, 0]

	if record_fname:
		record_fh = open_file_or_name(record_fname, "w", mkdir=True)
	#check line has an objet2
	for e in line.element_list:
		if ("OBJET2" in str(type(e)).split("'")[1]):
			objet = e
			break
	else:
		raise ValueError, "Line has no OBJET2 element"
	
	for e in line.element_list:
		if ("FAISCNL" or "FAISTORE" in str(type(e)).split("'")[1]):
			break
	else:
		raise ValueError, "Line has no FAISCNL element"

	current_YTZP = numpy.array(init_YTZP)
	areas = []
	coords = []
	tracks = []
	close_orbit_found = False
	for iteration in xrange(max_iterations):
		zlog.debug("start iteration: "+str(iteration)+ " with coords "+str(current_YTZP))
		coords.append(current_YTZP)
		objet.clear()	# remove existing particles
		objet.add(Y=current_YTZP[0], T=current_YTZP[1], Z=current_YTZP[2], P=current_YTZP[3], LET='A', D=D)

		r = line.run(xterm=False)

		if not r.run_success():
			if iteration == 0:
				zlog.warning("Not a stable orbit")
				return None
			else:
				zlog.warning("Center of orbit from last iteration, not stable")
				return None
		else:
			track = r.get_track('fai', ['Y', 'T', 'Z', 'P'])
			
			#use track data at location given by fai_label
			if fai_label != None:
				if iteration == 0:
					label = flatten(r.get_track('fai', ['element_label1']))
					label = [lab.rstrip() for lab in label]

					if len(find_indices(label, fai_label)) < 1:
						zlog.error("number of instances of label "+ str(fai_label) + " not equal 1")
						return

					fai_index = find_indices(label, fai_label)
				try:
					track = numpy.array([track[fi] for fi in fai_index])
				except:
					track = [track[fai_index]]

		track_a = numpy.zeros([len(track)+1, len(track[0])])
		track_a[0] = current_YTZP
		track_a[1:] = track
		
		if record_fname:
			record_fh.write("#track %s\n"%track_a.shape[0])
			numpy.savetxt(record_fh, track_a)

		tracks.append(track_a)

		#clean up tmp directory
		line.clean()
		
		centre_h = find_centre(track_a[:, 0:2])
		area_h = calc_area_simple(track_a[:, 0:2], centre=centre_h)
		centre_v = find_centre(track_a[:, 2:4])
		area_v = calc_area_simple(track_a[:, 2:4], centre=centre_v)

		zlog.debug("End iteration: "+str(iteration)+ " final coords "+str(track_a[-1])+", center "+str([centre_h, centre_v])+", area "+str([area_h,area_v]))


		areas.append([area_h , area_v])
		prev_YTZP = current_YTZP
		current_YTZP = numpy.array([centre_h[0], centre_h[1], centre_v[0], centre_v[1]])
		
		difs = numpy.zeros(4)
		for x in xrange(4):
			if abs(prev_YTZP[x]) < tol or abs(current_YTZP[x]) < tol:
				difs[x] = abs(prev_YTZP[x] - current_YTZP[x])
			else:
				difs[x] = abs((prev_YTZP[x] - current_YTZP[x])/ prev_YTZP[x])

		if difs.max() < tol:
			close_orbit_found = True
			close_orbit = current_YTZP
			break

		#if area_h < tol and area_v < tol:
		#	close_orbit_found = True
		#	close_orbit = current_YTZP
		#	break


	
	for x in xrange(len(areas)):
		print "it:", x, "\tYTZP", coords[x], "\tarea", areas[x]


	if close_orbit_found:
		print "found closed orbit"
		print "Y=%s, T=%s, Z=%s, P=%s" % tuple(current_YTZP)
		return current_YTZP
	else:
		zlog.warn("Iterations did not converge, no closed orbit found")
		return None


def plot_find_closed_orbit(data_fname, outfile=None):
	"""When the closed orbit search fails it can be useful to see what happened. find_closed_orbit() and find_closed_orbit_range() can take an optional argument record_fname. This causes them to write a log file, which can be read by this function and plotted. ::

		closed_orbit =  find_closed_orbit(tline, init_YTZP=[1,20,0,0], record_fname="res/closedorbit.log")
		plot_find_closed_orbit(data_fname="res/closedorbit.log", outfile="res/closedorbit.png")
	
	This will show you the succession of guesses. If using plot_find_closed_orbit_range() it will also show a plot of all the coordinates that survive being tracked through the lattice.

	From the plots it may be clear whether the lattice is unstable, or if the convergence is just too slow. In the latter case it may help to repeat the cell using REBELOTE, or to increase the number of iterations allowed.

	"""
	import pylab
	all_data = open(data_fname).readlines()
	outfile_pre, dummy, outfile_ext = outfile.rpartition(".")

	line_n = 0
	tracks = []
	survivors = None
	while line_n < len(all_data):
		line = all_data[line_n]
		if line.startswith("#track"):
			dummy, track_len = line.split()
			track_len = int(track_len)

			track_str = "\n".join(all_data[line_n+1: line_n+track_len+1])
			#print track_str
			track = numpy.genfromtxt(StringIO(track_str))
			#print track

			line_n = line_n+track_len
			tracks.append(track)
		if line.startswith("#survivors"):
			dummy, surv_len = line.split()
			surv_len = int(surv_len)

			surv_str = "\n".join(all_data[line_n+1: line_n+surv_len+1])
			surv = numpy.genfromtxt(StringIO(surv_str))

			line_n = line_n+surv_len
			survivors = surv

		line_n += 1

	for axis in ['h', 'v']:
		pylab.clf()
		pylab.title("closed orbit searches")

		for track in tracks[:]:
			if axis=='h':
				x, y = track[:,0], track[:,1]
				pylab.xlabel("x")
				pylab.ylabel("x'")
			if axis=='v':
				x, y = track[:,2], track[:,3]
				pylab.xlabel("y")
				pylab.ylabel("y'")
			#print
			#print x, y
			pylab.plot(x,y, 'x')
		
		
		if outfile != None:
			pylab.savefig(outfile_pre+"track_"+axis+"."+outfile_ext)
		else:
			pylab.show()

		pylab.clf()
		pylab.title("surviving particle start and ends")

		if survivors != None:
			for surv in survivors[:10]:
				if axis=='h':
					x, y, x1, y1 = surv[[0,1,4,5]]
					pylab.xlabel("x")
					pylab.ylabel("x'")
				if axis=='v':
					x, y, x1, y1 = surv[[2,3,6,7]]
					pylab.xlabel("y")
					pylab.ylabel("y'")
				pylab.plot([x], [y], 'x')
				pylab.plot([x,x1], [y,y1])


			if outfile != None:
				pylab.savefig(outfile_pre+"survivor_"+axis+"."+outfile_ext)
			else:
				pylab.show()


def get_twiss_profiles(line, file_result, input_twiss_parameters=None):
	""" Calculates the twiss parameters at all points written to zgoubi.plt. 11 particle trajectories are used to calculate the
	transfer matrix, just as is done in zgoubi. The code mirrors that found in mat1.f, mksa.f. The twiss parameters are first
	calculated at the end of the cell using either input_twiss_parameters (format [beta_y, alpha_y, gamma_y, beta_z, alpha_z, gamma_z]) 
	or, if this is not supplied, it assumes the cell is periodic and uses get_twiss_parameters to find the boundary condition. 

	The twiss parameters are then mapped to all points in the magnets using the transfer matrix calculated at each point. The results are stored
	in list twiss_profiles where each row represents a point in the zgoubi.plt file and has format 

	[s_coord, label,  mu_y, beta_y, alpha_y, gamma_y, mu_z, beta_z, alpha_z, gamma_z]

	Requires an OBJET type 5, and a MATRIX element.

	Note - This calculation uses trajectories as measured in the local coordinate system of the magnet.
	
	The table is also writen to file_result"""

	if input_twiss_parameters == None:
		input_twiss_parameters = [0, 0, 0, 0, 0, 0]

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
	try:
		try:
			plt_track = r.get_track('plt', ['LET', 'D0', 'Y0', 'T0', 'Z0', 'P0', 'X0', 'D', 'Y', 'T', 'Z', 'P', 'S', 'X'])
			track_type = 'plt'
		except IOError:
			plt_track = r.get_track('fai', ['LET', 'D0', 'Y0', 'T0', 'Z0', 'P0', 'X0', 'D', 'Y', 'T', 'Z', 'P', 'S'])
			track_type = 'fai'
	except ValueError:
		# new file format has no 'X0', now called 'S0' to match zgoubi output
		# also 'D' is now know as 'D-1'
		try:
			plt_track = r.get_track('plt', ['LET', 'D0-1', 'Y0', 'T0', 'Z0', 'P0', 'S0', 'D-1', 'Y', 'T', 'Z', 'P', 'X', 'S'])
			track_type = 'plt'
		except IOError:
			plt_track = r.get_track('fai', ['LET', 'D0-1', 'Y0', 'T0', 'Z0', 'P0', 'S0', 'D-1', 'Y', 'T', 'Z', 'P', 'S'])
			track_type = 'fai'

	transpose_plt_track = map(list, zip(*plt_track))
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
	if track_type == 'plt':
		X = transpose_plt_track[13]
		#read labels
		label = [x[0].strip() for x in r.get_track('plt', ['element_label1'])]
	else:
		X = [0 for x in xrange(len(S))]
		label = [x[0].strip() for x in r.get_track('fai', ['element_label1'])]


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
	ref_indices = find_indices(track_tag, 'O')
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
		indices = find_indices(track_tag, track_tag_name)
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
	R16_list = [x/DP for x in map(numpy.subtract, Y_alltracks[9], Y_alltracks[10])]
	R26_list = [x/DP for x in map(numpy.subtract, T_alltracks[9], T_alltracks[10])]
	R36_list = [x/DP for x in map(numpy.subtract, Z_alltracks[9], Z_alltracks[10])]
	R46_list = [x/DP for x in map(numpy.subtract, P_alltracks[9], P_alltracks[10])]
	R56_list = [x/DP for x in map(numpy.subtract, S_alltracks[9], S_alltracks[10])]
	#-----------------------------------------------------------------------------------------

	#assume it1=1
	it1 = 1
	for i in range(1, 5):
		i2 = 2*i +it1-1
		i3 = i2 + 1
		#adjust indices since python lists are numbered from zero
		i2 = i2 - 1
		i3 = i3 - 1
		if i == 1:
			UO = Y0_alltracks[i2]-Y0_alltracks[i3]
			R11_list =  [x/UO for x in map(numpy.subtract, Y_alltracks[i2], Y_alltracks[i3])]
			R21_list =  [x/UO for x in map(numpy.subtract, T_alltracks[i2], T_alltracks[i3])]
			R31_list =  [x/UO for x in map(numpy.subtract, Z_alltracks[i2], Z_alltracks[i3])]
			R41_list =  [x/UO for x in map(numpy.subtract, P_alltracks[i2], P_alltracks[i3])]
			R51_list =  [x/UO for x in map(numpy.subtract, S_alltracks[i2], S_alltracks[i3])]
		elif i == 2:
			UO = T0_alltracks[i2]-T0_alltracks[i3]
			R12_list =  [x/UO for x in map(numpy.subtract, Y_alltracks[i2], Y_alltracks[i3])]
			R22_list =  [x/UO for x in map(numpy.subtract, T_alltracks[i2], T_alltracks[i3])]
			R32_list =  [x/UO for x in map(numpy.subtract, Z_alltracks[i2], Z_alltracks[i3])]
			R42_list =  [x/UO for x in map(numpy.subtract, P_alltracks[i2], P_alltracks[i3])]
			R52_list =  [x/UO for x in map(numpy.subtract, S_alltracks[i2], S_alltracks[i3])]
		elif i == 3:
			UO = Z0_alltracks[i2]-Z0_alltracks[i3]
			R13_list =  [x/UO for x in map(numpy.subtract, Y_alltracks[i2], Y_alltracks[i3])]
			R23_list =  [x/UO for x in map(numpy.subtract, T_alltracks[i2], T_alltracks[i3])]
			R33_list =  [x/UO for x in map(numpy.subtract, Z_alltracks[i2], Z_alltracks[i3])]
			R43_list =  [x/UO for x in map(numpy.subtract, P_alltracks[i2], P_alltracks[i3])]
			R53_list =  [x/UO for x in map(numpy.subtract, S_alltracks[i2], S_alltracks[i3])]
		elif i == 4:
			UO = P0_alltracks[i2]-P0_alltracks[i3]
			R14_list =  [x/UO for x in map(numpy.subtract, Y_alltracks[i2], Y_alltracks[i3])]
			R24_list =  [x/UO for x in map(numpy.subtract, T_alltracks[i2], T_alltracks[i3])]
			R34_list =  [x/UO for x in map(numpy.subtract, Z_alltracks[i2], Z_alltracks[i3])]
			R44_list =  [x/UO for x in map(numpy.subtract, P_alltracks[i2], P_alltracks[i3])]
			R54_list =  [x/UO for x in map(numpy.subtract, S_alltracks[i2], S_alltracks[i3])]
		
#!----------------------------------
#! Adjust units to SI (as in mksa.f)
#!----------------------------------
	unit_list = [1e-2, 1e-3, 1e-2, 1e-3, 1e-2, 1, 1e-6]
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

	if input_twiss_parameters == [0, 0, 0, 0, 0, 0]:
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
	
	zlog.debug("Initial parameters:\nbeta_y_0, alpha_y_0, gamma_y_0, beta_z_0, alpha_z_0, gamma_z_0\n%s, %s, %s, %s, %s, %s" % (beta_y_0, alpha_y_0, gamma_y_0, beta_z_0, alpha_z_0, gamma_z_0))
	if beta_y_0 == 0 or beta_z_0 == 0:
		zlog.error("Beam is unstable")


#! Calculate twiss parameters at all points in plt file 
#!######################################################

#! The twiss parameters at the end of the periodic cell have been calculated by call to get_twiss_parameters, otherwise twiss parameters are input 
#!
#! - e.g. S.Y.Lee Accelerator physics equation 2.54
#! - The phase advance can be found by applying a Floquet transformation to the transfer matrix (S.Y.Lee eqn 2.65)

 
	fresults = open(file_result, 'w')
              
	mu_y_list = []
	beta_y_list = []
	alpha_y_list = []
	gamma_y_list = []
	mu_z_list = []
	beta_z_list = []
	alpha_z_list = []
	gamma_z_list = []
	sign_sine_y_old = 1.0
	sign_sine_z_old = 1.0
	n_pi_y = 0
	n_pi_z = 0
	try:
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
			print >> fresults, '%2f %2s %2f %2f %2f %2f %2f %2f %2f %2f' % (S_alltracks[0][i], label_ref[i], mu_y_list[i], beta_y, alpha_y_list[i], \
				gamma_y_list[i], mu_z_list[i], beta_z, alpha_z_list[i], gamma_z_list[i])
	except (IndexError, ZeroDivisionError, ValueError):
		print "Error calculating twiss parameters from twiss matrix"
		print "i=", i, "(", label[i], ")"
		print "Matrix:"
		for ai in range(1,7):
			for aj in range(1,7):
				if aj == 5 or ai == 6:
					print "-\t",
					continue
				print eval("R%s%s_list[i]"%(ai,aj)),"\t",
			print
		raise
	#put twiss parameters together. Format [s_coord, mu_y, beta_y, alpha_y, gamma_y, mu_z, beta_z, alpha_z, gamma_z]. Units are SI
	twiss_profiles = [[s*cm for s in S_alltracks[0]], label_ref, mu_y_list, beta_y_list, alpha_y_list, gamma_y_list, mu_z_list, beta_z_list, alpha_z_list, gamma_z_list]

	return twiss_profiles


def fourier_tune(line, initial_YTZP, D_in, nfourierturns, plot_fourier=False, coords=None):
	"""Calculate tune using FFT. nfourierturns determines the number of passes through the lattice.
	Can supply set of horizontal and vertical coordinates in coords = [ycoords,zcoords], otherwise
	routine will calculate coordinates

	Set plot_fourier= True to show Fourier spectrum. Default is False
	"""
	if coords == None:
		coords = []

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
			reb = e
			break
	else:
		raise ValueError, "Line has no REBELOTE element"

	if coords == []:
		objet.clear()	# remove existing particles
		#start at closed orbit with a small amplitude added to vertical component to get tune readings
		objet.add(Y=initial_YTZP[0], T=initial_YTZP[1], Z=initial_YTZP[2]+1e-5, P=initial_YTZP[3], LET='A', D=D_in)

		reb.set(NPASS=nfourierturns-1)	
		r = line.run(xterm = False)
		YZ = r.get_track('fai', ['Y', 'Z'])
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

	ypeaksloc = sorted([ysortfreqamp[0][0], ysortfreqamp[1][0]])
	zpeaksloc = sorted([zsortfreqamp[0][0], zsortfreqamp[1][0]])

	yfouriertune = ypeaksloc[0]/len(ycoords)
	zfouriertune = zpeaksloc[0]/len(zcoords)

	#plot tunes if desired	
	if(plot_fourier):
		pylab.subplot(211)
		pylab.plot(yampfreq, 'k-')
		pylab.ylim((0, max(yampfreq[1:-1])))
		pylab.ylabel("FFT(y)")
		pylab.subplot(212)
		pylab.plot(zampfreq, 'b-')
		pylab.ylim((0, max(zampfreq[1:-1])))
		pylab.ylabel("FFT(z)")
		pylab.savefig('fspectrum')

	return yfouriertune, zfouriertune

def scan_dynamic_aperture(line, emit_list_h, emit_list_v, closedorb_YTZP, npass, D_mom, beta_gamma_input = 1, ellipse_coords = 1, coord_pick = 0, twiss_parameters = [], plot_data = False):
	""" Check a list of emittances (emit_list, units Pi m rad) to see if tracking succeeds. Can be used to establish the dynamic aperture. If the elements in emit_list are increasing then will stop tracking once it finds the lowest emittance where a particle is lost. On the other hand, if the elements in emit_list are decreasing then will stop tracking once it reaches the first point where tracking succeeds without loss. 

	Required input:
		emit_list_h - List of emittances in horizontal plane to check
		emit_list_v - List of emittances in vertical plane to check

		closedorb_YTZP - Closed orbit coordinates as returned by find_closed_orbit

		npass - Number of passes through lattice

		D_mom - Momentum factor p/pref

	Optional input
		beta_gamma_input - If this is supplied then the emittances supplied in emit_list are assumed to be normalised and a conversion to geometrical emittance is made. Otherwise the emittances supplied are assumed to be geometrical.


		ellipse_coords - If greater than 1 will test uniformly distributed set of coordinates around around phase space
		ellipse. The points are generated by emittance_to_coords. ellipse_coords = 1 (or -1) will take single point where phase space cuts either the y, t, z or p axes (decided by coord_pick). ellipse_coords = -1 selects the negative point in the phase space ellipse. The default settings ellipse_coords = 1, coord_pick = 0 selects the point where the phase space ellipse cuts the y, z axis to the outside of the closed orbit point.

		coord_pick - Choose one of ellipse coords, i.e. ellipse_coords[coord_pick] is selected. If ellipse_coords equal 1 or -1, then coords_pick = 0 selects the points where the phase space ellipse crosses the Y, Z axes. If coord_pick = 1 then selects the point where the phase space ellipse crosses the T, P axes. If coord_pick = None, will select all points generated by emittance_to_coords.

		Can also specify ellipse coords = [n,m] -- distribute n points around ellipse but only test point m of these.

		twiss_parameters - Can supply twiss parameters in format [beta_y, alpha_y, gamma_y, disp_y, disp_py, beta_z, alpha_z, gamma_z, disp_z, disp_pz], i.e the output of r.get_twiss_parameters(). If not supplied an attempt will be made to calculate the twiss parameters in this function.

		plot_data - If True, creates phase space plots at all emittances scanned in both transverse planes.
	
		Returns [index_lost, coord_index], fourier_tune_emit, coords_YTZP_ini]  where 
		index_lost is the index in emit_list where particle is lost (0 if no loss). 
		coord_index indicates at which coord in ellipse_coords the particle is lost. 
		fourier_tune_emit is a list of tunes found at each emittance using Fourier analysis. 
		coords_YTZP_ini_list is a list of all initial coordinates tested. 

		If a particle is lost, returns index_lost in emit_list where the loss occurs. Otherwise index_lost remains 0.
	"""

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


	if twiss_parameters != []:
		betayz = [twiss_parameters[0], twiss_parameters[5]]
		gammayz = [twiss_parameters[2], twiss_parameters[7]]
	else:
		#calculate optical parameters on closed orbit. This is required to convert emittances into a coordinate.
		objet5 = zg.OBJET5()
		line.replace(objet, objet5)
		objet5.set(BORO=rigidity)
		objet5.set(PY=1e-4, PT=1e-3, PZ=1e-4, PP=1e-3, PX=1e-3, PD=1e-3)
		objet5.set(YR=closedorb_YTZP[0], TR=closedorb_YTZP[1], ZR=closedorb_YTZP[2], PR=closedorb_YTZP[3], DR=D_mom)
		matrix = zg.MATRIX(IORD=1, IFOC=11)
		line.replace(reb, matrix)
		r = line.run(xterm = False)
		twissparam = r.get_twiss_parameters()
		#alphayz = [twissparam[1],twissparam[4]]
		betayz = [twissparam[0], twissparam[5]]
		gammayz = [twissparam[2], twissparam[7]]
	
		#revert to objet2 mode with rebelote
		line.replace(objet5, objet)
		objet.clear()	# remove existing particles
		objet.add(Y=closedorb_YTZP[0], T=closedorb_YTZP[1], Z=closedorb_YTZP[2], P=closedorb_YTZP[3], LET='A', D=D_mom)
		line.replace(matrix, reb)


	index_lost = None
	coord_index_lost = None
	from operator import add
	if plot_data:
		YTZP_list = []

	fourier_tune_emit = []

	#if emit_list is a list of decreasing emittances, looking for the first point where particle is not lost
	reverse_search = False
	if emit_list_h[0] > emit_list_h[-1]:
		print "decreasing list (emit_list_h) "
		reverse_search = True
		

	coords_YTZP_ini_list = []
	for emit_h, emit_v in zip(emit_list_h, emit_list_v):
		print "check emit (h/v) ", emit_h, emit_v
		print "ellipse_coords, coord_pick ",ellipse_coords, coord_pick

		if coord_pick == None:
			#obtain coordinates on phase space ellipses. Use beta, gamma values found at closed orbit
			coords_YTZP_ini = emittance_to_coords(emit_h, emit_v, gammayz, betayz, beta_gamma_input, ncoords = abs(ellipse_coords))
		else:
			#select one point on phase space ellipse. Use beta, gamma values found at closed orbit
			coords_YTZP_ini = emittance_to_coords(emit_h, emit_v, gammayz, betayz, beta_gamma_input, ncoords = abs(ellipse_coords))[coord_pick]

		try:
			l = len(coords_YTZP_ini[0])
			if ellipse_coords != -1:
				coords_YTZP_ini = [map(add, closedorb_YTZP, coords) for coords in coords_YTZP_ini]
			else:
				coords_YTZP_ini = [[a-b for a,b in zip(closedorb_YTZP, coords)] for coords in coords_YTZP_ini]
		except TypeError:
			if ellipse_coords != -1:
				coords_YTZP_ini = [map(add, closedorb_YTZP, coords_YTZP_ini)]
			else:
				coords_YTZP_ini = [[a-b for a,b in zip(closedorb_YTZP,coords_YTZP_ini)]]


		coords_YTZP_ini_list.append(coords_YTZP_ini)

		for coord_index, current_YTZP in enumerate(coords_YTZP_ini):

			print "ellipse coord index ", coord_index
			print "current coords to test ",current_YTZP

			objet.clear()	# remove existing particles
			objet.add(Y=current_YTZP[0], T=current_YTZP[1], Z=current_YTZP[2], P=current_YTZP[3], LET='A', D=D_mom)

			#pure-x
			#objet.add(Y=current_YTZP[0], T=current_YTZP[1], Z= 0.0, P= 0.0, LET='A', D=D_mom)
				
			#run Zgoubi
			r = line.run(xterm= False)
    
			rebelote_completed = r.test_rebelote()

			lost = not rebelote_completed
    
			if(lost):
				print "tracking failed at emit_h, emit_v = ", emit_h, emit_v
				coord_index_lost = coord_index
				#break
				#return [index_lost, coord_index],fourier_tune_emit

			if plot_data and not lost:
				YTZP_list.append(r.get_track('fai', ['Y', 'T', 'Z', 'P']))
				coords_in = [flatten(r.get_track('fai', ['Y'])), flatten(r.get_track('fai', ['Z']))]
				fourier_tune_result = fourier_tune(line, [], 1, 1, coords = coords_in)
				fourier_tune_emit.append(fourier_tune_result)

		#condition to end scan
		if reverse_search and not lost:
			index_lost = emit_list_h.index(emit_h)+1
			break
		elif not reverse_search and lost:
			index_lost = emit_list_h.index(emit_h)
			break			


	else:
		print "Successful tracking in all test cases"

	if plot_data:

		Y_data = []
		T_data = []
		Z_data = []
		P_data = []
    
		for index in range(len(YTZP_list)):
			Y_data.append(numpy.transpose(YTZP_list[index])[0])
			T_data.append(numpy.transpose(YTZP_list[index])[1])
			Z_data.append(numpy.transpose(YTZP_list[index])[2])
			P_data.append(numpy.transpose(YTZP_list[index])[3])

		plot_data_xy_multi(Y_data, Z_data, 'yz_space', labels=["YZ coords", "y [cm]", "z [cm]"], style=['k+'])

		#obtain coordinates on phase space ellipse using closed orbit twiss parameters
		if lost:
			emit_plot_h = emit_list_h[index_lost-1]
			emit_plot_v = emit_list_v[index_lost-1]
		else: 
			emit_plot_h = emit_list_h[-1]
			emit_plot_v = emit_list_v[-1]

		coords_YTZP_full = emittance_to_coords(emit_plot_h, emit_plot_v, gammayz, betayz, beta_gamma_input, ncoords = 100)
		coords_YTZP_full = [map(add, closedorb_YTZP, coords) for coords in coords_YTZP_full]
		coords_YTZP_full.append(coords_YTZP_full[0])


		if coord_pick == None:
			#obtain coordinates on phase space ellipses. Use beta, gamma values found at closed orbit
			coords_YTZP_lim = emittance_to_coords(emit_plot_h, emit_plot_v, gammayz, betayz, beta_gamma_input, ncoords = abs(ellipse_coords))
		else:
			#select one point on phase space ellipse. Use beta, gamma values found at closed orbit
			coords_YTZP_lim = emittance_to_coords(emit_plot_h, emit_plot_v, gammayz, betayz, beta_gamma_input, ncoords = abs(ellipse_coords))[coord_pick]

		try:
			l = len(coords_YTZP_lim[0])
			if ellipse_coords != -1:
				coords_YTZP_lim = [map(add, closedorb_YTZP, coords) for coords in coords_YTZP_lim]
			else:
				coords_YTZP_lim = [[a-b for a,b in zip(closedorb_YTZP, coords)] for coords in coords_YTZP_lim]
		except TypeError:
			if ellipse_coords != -1:
				coords_YTZP_lim = [map(add, closedorb_YTZP, coords_YTZP_lim)]
			else:
				coords_YTZP_lim = [[a-b for a,b in zip(closedorb_YTZP,coords_YTZP_lim)]]


		#coords_YTZP_ini = [map(add, closedorb_YTZP, coords) for coords in coords_YTZP_ini]

		#add points on phase space ellipse actually used to initialise scan above
		Y_data.insert(0,numpy.transpose(coords_YTZP_lim)[0])
		T_data.insert(0,numpy.transpose(coords_YTZP_lim)[1])
		Z_data.insert(0,numpy.transpose(coords_YTZP_lim)[2])
		P_data.insert(0,numpy.transpose(coords_YTZP_lim)[3])

		#add many coordinates to draw phase space ellipse
		Y_data.insert(0,numpy.transpose(coords_YTZP_full)[0])
		T_data.insert(0,numpy.transpose(coords_YTZP_full)[1])
		Z_data.insert(0,numpy.transpose(coords_YTZP_full)[2])
		P_data.insert(0,numpy.transpose(coords_YTZP_full)[3])

		style_list =  ['k-', 'ro']
		style2 = ['b+','r+', 'g+', 'm+', 'y+']
		for i in range(len(Y_data)-2):
			style_list.append(style2[i%(len(style2))])

		plot_data_xy_multi(Y_data, T_data, 'yt_phasespace', labels=["Horizontal phase space", "y [cm]", "y' [mrad]"], style=style_list)
		plot_data_xy_multi(Z_data, P_data, 'zp_phasespace', labels=["Vertical phase space", "z [cm]", "z' [mrad]"], style=style_list)

	return [index_lost, coord_index], fourier_tune_emit, coords_YTZP_ini_list



def emittance_to_coords(emit_horizontal, emit_vertical, gammayz, betayz, beta_gamma_input = 1, ncoords = 1):
	"""Given some emittance in horizonal and vertical space
	
	If ncoords <= 1 return points where phase space ellipse crosses the y,y' and z,z' axis.
	
	If ncoords > 1, will instead give a distribution of points (y,t) around the phase space ellipse uniform in angle theta where::
	
		y = a*cos(theta)*cos(phi) - b*sin(theta)*sin(phi)
		t = a*cos(theta)*sin(phi) + b*sin(theta)*cos(phi)
	
	where (a,b) are major and minor radii and phi is the tilt of the major radius w.r.t horizontal axis. A (z,p) distribution is 
	similarly calculated.
	
	Can use these points, returned in coords_YTZP to start tracking at the desired emittance (assuming that the optical functions don't change with amplitude).
	
	Emittances in both the horizontal and vertical planes may be supplied. Twiss parameters beta and gamma in 
	both places may be determined calling get_twiss_parameters beforehand i.e.::
	
		twissparam = r.get_twiss_parameters()
		betayz = [twissparam[0],twissparam[5]]
		gammayz = [twissparam[2],twissparam[7]]
	"""
	
	coords_YTZP = []

	if ncoords <= 1:
		#crosses y or z axis at sqrt(emit/gamma). Convert to cm.
		y_axis = cm_*sqrt(emit_horizontal/(gammayz[0]*beta_gamma_input))       
		z_axis = cm_*sqrt(emit_vertical/(gammayz[1]*beta_gamma_input))
		coords_YTZP.append([y_axis, 0, z_axis, 0])

		#crosses y' or z' axis at sqrt(emit/beta) and convert to mrad (same as scaling by mm)
		t_axis = mm_*sqrt(emit_horizontal/(betayz[0]*beta_gamma_input))       
		p_axis = mm_*sqrt(emit_vertical/(betayz[1]*beta_gamma_input))
		coords_YTZP.append([0, t_axis, 0, p_axis])

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
			coords_YTZP.append([ydat[index], tdat[index], zdat[index], pdat[index]])

	return coords_YTZP


def find_indices(a_list, list_element):
	""" Find all occurences of list_element in list. Return the indices.
	"""
	indices = []
	i = -1
	try:
		while 1:
			i = a_list.index(list_element, i+1)
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

def scaling_to_poly(b0, k, r0, rmin, rmax, step, order=5):
	"""In a scaling FFAG, the magnetic field follows the scaling law given by B=B0*(r/r0)^k where r is the radial coordinate,
	   B0 is the field at r=r0 and k is the scaling factor.
	   This def fits a polynomial of a given order to the scaling field at points in the region rmin < r < rmax, where the number of points is determined by the increment step. The coefficients of 
	   polynomial are returned.

	   Required input - 
	   b0 - the field at r0
	   r0 - the reference radius
	   rmin,rmax - region of fit
	   step - step size in (rmin,rmax) at which scaling field is evaluated

	   Optional input-
	   order - Order of polynomial fit (default 5)
	   
	   Author -  S. Sheehy """

	x = numpy.arange(rmin, rmax, step)
	def scalelaw(a): return b0*pow((r0+a)/r0, k)
	By = map(scalelaw, x)
	polyval = numpy.polyfit(x, By, order)
	vals = polyval.tolist()
	vals.reverse()
	while len(vals)<6:
		vals.append(0)
	return vals

def get_enclosing_circle(ellipse_data):
	""" Find smallest circle that encloses a set of ellipses centred on the midplane. The ellipses are defined by their horizontal and vertical radii and by their centre along the horizontal axis (a,b,c). The algorithm optimised both the centre of the enclosing circle and its radius. 

	Ellipse algorithm by J. Scott Berg (BNL)
	
	Reference - J. Scott Berg. 'Finding the Circular Magnet Aperture which Encloses and Arbitrary number of midplace-centered Beam Ellipses', Proc. EPAC04, 5-9 July, Lucerne Congress Centre, Lucerne, Switzerland

	""" 

	import zgoubi.ellipse

	bc = zgoubi.ellipse.BestCircle()

	for edata in ellipse_data:
		a = edata[0]
		b = edata[1]
		c = edata[2]
		bc.append((a, b, c))

	(centre, radius) = bc.get_circle()
	print 'centre of enclosing circle =', centre, ", with radius =", radius

	return centre, radius


def misalign_element(line, element_indices, mean, sigma, sigma_cutoff, misalign_dist=None, seed = None):
	""" Misalign the set of elements identified by element_indices. The misalignment is Gaussian
	with sigma, sigma_cutoff and mean among the input parameters. If sigma = 0.0, the misalignment is constant 
	and given by mean

	Alternatively, misalignments specified by list misalign_dist

	Returns the misailgnment distribution (misalign_dist), unit meters

	element_indices can be determined using line.find_elements(element_name)

	mean and sigma should be input in meters """
	if misalign_dist == None:
		misalign_dist = []

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
		misalign_dist = gaussian_cutoff(len(element_indices), mean, sigma, sigma_cutoff, seed=seed)

	#insert CHANGREF with misalignments given by misalign_dist.
	for index, pos in enumerate(changref_indices):
		misalign_elem = zg.CHANGREF('misalign', YCE=misalign_dist[index]*cm_)
		misalign_elem_rev = zg.CHANGREF('misalign', YCE=-misalign_dist[index]*cm_)
		line.insert(pos, misalign_elem)
		line.insert(pos+2, misalign_elem_rev)

	return misalign_dist


def gaussian_cutoff(npoints, mean, sigma, sigma_cutoff, seed = None):
	""" Generate Gaussian distribution about mean with standard deviation sigma
          npoints - number of points in distribution required
          sigma_cutoff - the cutoff factor i.e. points outside sigma*sigma_cutoff eliminated

          Seed - Optional parameter to set random seed"""

	import random

	if seed != None:
		random.seed(seed)

	dist = []
	while len(dist) < npoints:
		if sigma > 0.0:
			rand = random.gauss(0, sigma)
			if abs(rand) < sigma_cutoff*sigma or sigma_cutoff == 0:
				dist.append(mean+rand)
		else:
			dist.append(mean)

	return dist


def tune_diagram(tune_list, order=3, xlim=None, ylim=None):
	"""Plot a list of tunes on a tune diagram with resonance line up to a given order.
	
	Required input 
		tune_list - format [[list of horizontal tunes],[list of vertical tunes]]
	
	Optional
		order - Integer denoting order up to which resonance lines are drawn. Default value 3 (third order).
		xlim = [lower_value, upper_value] - Limit horizontal axis
		ylim = [lower_value, upper_value] - Limit vertical axis  
	
	At the moment, since the tune diagram covers the range [0,1] just fractional tunes can be shown.
	The default value of order is 3. 
	Written by Y. Giboudet
	"""
	if xlim == None:
		xlim = [0, 1]
	if ylim == None:
		ylim = [0, 1]

	import pylab

	XL = xlim[0]
	XR = xlim[1]
	YB = ylim[0]
	YT = ylim[1]
	NY = 0
	LOC = [0, 0]
	col = ['b', 'r', 'g', 'm', 'y', 'k', 'c']*3
	for I in xrange(1, order+1):
		for J in xrange(-I, I+1):
			NX = J
			if (NX>=0) and (NY>=0):
				NY =  I-NX
			elif (NX>=0) and (NY<0):
				NY = -I+NX
			elif (NX<0) and (NY>=0):
				NY =  I+NX
			elif (NX<0) and (NY<0):
				NY = -I-NX
			for K in xrange(-I+1, I):
				X1 = -1
				X2 = -1
				Y1 = -1
				Y2 = -1
				
				if NX != 0:
					if (float(K)/NX > 0) and (float(K)/NX <= 1):
						X1 = K/float(NX)
					if ((K-NY)/float(NX) >= 0) and ((K-NY)/float(NX) < 1):
						X2 = (K-NY)/float(NX)
				if NY != 0:
					if (float(K)/NY >= 0) and (float(K)/NY < 1):
						Y1 = float(K)/NY
					if ((K-NX)/float(NY) > 0) and ((K-NX)/float(NY) <= 1):
						Y2 = (K-NX)/float(NY)
				XM = 0
				if X1 != -1:
					if Y2 != -1 :
						pylab.plot([X1, 1], [0, Y2], color=col[I])    

					elif Y1 != -1:
						pylab.plot([X1, 0], [0, Y1], color=col[I])
						
					elif X2 != -1:
						pylab.plot([X1, X2], [0, 1], color=col[I])

				if Y1 != -1:
					if Y2 != -1 :
						pylab.plot([0, 1], [Y1, Y2], color=col[I])
						

					elif X2 != -1:
						pylab.plot([0, X2], [Y1, 1], color=col[I])
						
				if X2 != -1 :
					if Y2 != -1:
						pylab.plot([X2, 1], [1, Y2], color=col[I])
						
						
	#plot box around area of resonance lines
	pylab.plot([XL,    XL+XL], [YB,    YB+YT], color='k')
	pylab.plot([XL+XL, XL+XR], [YB+YT, YB+YT], color='k')
	pylab.plot([XL+XR, XL+XR], [YB+YT, YB+YB], color='k')
	pylab.plot([XL+XR, XL+XL], [YB+YB, YB+YB], color='k')

	#plot list of tunes
	pylab.plot(tune_list[0], tune_list[1], 'k+')

	pylab.axis([xlim[0], xlim[1], ylim[0], ylim[1]])

	pylab.savefig('tune_diagram')
	pylab.cla()
	


 
def plot_data_xy(data, filename, labels=None, style='b-', xlim=None, ylim=None):
	if labels == None:
		labels = ["", "", ""]

	import pylab
	data_a = numpy.array(data)
	pylab.hold(False)
	pylab.plot(data_a[:, 0], data_a[:, 1], style)
	pylab.title(labels[0])
	if xlim != None:
		pylab.xlim( (xlim[0], xlim[1]) )
	if ylim != None:
		pylab.ylim( (ylim[0], ylim[1]) )
	pylab.xlabel(labels[1])
	pylab.ylabel(labels[2])
	pylab.savefig(filename)

	pylab.cla()


def plot_data_xy_multi(data_x_list, data_y_list, filename, labels=None, style='', legend=' ', legend_location='best', xlim=None, ylim=None, tick_multiple = None):
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
		ylim = [lower_value, upper_value] - Limit vertical axis
		majorLocator = horizontal axis tick marks defined as multiple of majorLocator"""

	import pylab

	if labels == None:
		labels = ["", "", ""]
	#check if data is a single list or a list of lists
	single_x_data = False
	single_y_data = False

	try:
		dummy = len(data_x_list[0])
	except TypeError:
		single_x_data = True

	try:
		dummy = len(data_y_list[0])
	except TypeError:
		single_y_data = True

	pylab.hold(True)

	if type(style) != list:
		style = [style]
	if style == ['']:
		style = ['k-', 'b-', 'r-', 'g-', 'm-', 'y-']

	if type(legend) != list:
		legend = [legend]


	if tick_multiple != None:
		from matplotlib.ticker import MultipleLocator, FormatStrFormatter
		majorLocator   = MultipleLocator(tick_multiple)
		ax = pylab.subplot(111)

	if single_y_data:
		pylab.plot(data_x_list, data_y_list, style[0])
	else:
		for index, data_y in enumerate(data_y_list):
			if single_x_data:
				pylab.plot(data_x_list, data_y, style[index%len(style)], label=legend[index%len(legend)])
			else:
				pylab.plot(data_x_list[index], data_y, style[index%len(style)], label= legend[index%len(legend)])

	pylab.title(labels[0])
	if xlim != None:
		pylab.xlim( (xlim[0], xlim[1]) )
	if ylim != None:
		pylab.ylim( (ylim[0], ylim[1]) )
	pylab.xlabel(labels[1])
	pylab.ylabel(labels[2])

	if legend != [' ']:
		pylab.legend(loc=legend_location)

	if tick_multiple != None:
		ax.xaxis.set_major_locator(majorLocator)

	pylab.savefig(filename)
	pylab.cla()



def calc_transfer_matrix(start_bunch, end_bunch):
	"""Track a bunch generated with OBJET5 through a line, and pass the start and end bunch to this function to calculate the twiss matrix. No unit conversion is done, so units match the passed bunch.
	NOT COMPLETE:
	only use cells in top 4 rows
	"""
	try:
		start = start_bunch.particles()
	except AttributeError:
		start = start_bunch
	try:
		end = end_bunch.particles()
	except AttributeError:
		end = end_bunch

	co = list(" DYTZPS") # coordinates

	tm = numpy.zeros([6,6])
	tm[4,4] = 1
	tm[5,5] = 1

	IT1 = 0 # offset of particles A to J
	I10 = IT1+9
	I11 = IT1+10
	# FO(1,I10) => start['D'][I10] 
	DP = (start['D'][I10] - start['D'][I11] ) / 0.5 /( start['D'][I10] + start['D'][I11])

	for j in xrange(2,6):
		tm[j-2, 5] = (end[co[j]][I10] - end[co[j]][I11]) / DP
		#print co[j], I10, I11, (end[co[j]][I10] - end[co[j]][I11]) / DP

		for i in xrange(1,5):
			i2 = 2*i + IT1-1
			i3 = i2 + 1
			u0 = start[co[i+1]][i2] - start[co[i+1]][i3]
			tm[j-2, i-1] = (end[co[j]][i2] - end[co[j]][i3]) / u0
			#print co[j],i2, i3, (end[co[j]][i2] - end[co[j]][i3]) / u0
			if (j == 5):
				tm[4,i-1] = (end[co[6]][i2] - end[co[6]][i3]) /u0
				#print co[6], i2, i3, (end[co[6]][i2] - end[co[6]][i3]) /u0
	tm[4,5] = (end[co[6]][I10] - end[co[6]][I11] ) /DP

	if (tm[0,0] + tm[1,1] > 2) or (tm[2,2] + tm[3,3] > 2):
		zlog.warning("Lattice is unstable")
	
	return tm


def calc_twiss_from_matrix(trans_matrix):
	"Calculate the twiss parameters (beta_y, alpha_y, gamma_y, beta_z, alpha_z, gamma_z) from a transfer matrix. Either use Results.get_transfer_matrix() or calc_transfer_matrix() to get matrix."
	tm = trans_matrix

	mu_y = acos(0.5 * (tm[0,0]+tm[1,1]))
	beta_y = tm[0,1]/sin(mu_y)
	alpha_y = (tm[0,0]-cos(mu_y))/sin(mu_y)
	gamma_y = - tm[1,0]/sin(mu_y)

	mu_z = acos(0.5 * (tm[2,2]+tm[3,3]))
	beta_z = tm[2,3]/sin(mu_z)
	alpha_z = (tm[2,2]-cos(mu_z))/sin(mu_z)
	gamma_z = - tm[3,2]/sin(mu_z)
	return (beta_y, alpha_y, gamma_y, beta_z, alpha_z, gamma_z)

def calc_phase_ad_from_matrix(trans_matrix):
	"Calculate the phase advance (mu_y,mu_z) from a transfer matrix. Either use Results.get_transfer_matrix() or calc_transfer_matrix() to get matrix."
	tm = trans_matrix
	mu_y = acos(0.5 * (tm[0,0]+tm[1,1]))
	mu_z = acos(0.5 * (tm[2,2]+tm[3,3]))
	return (mu_y,mu_z)
