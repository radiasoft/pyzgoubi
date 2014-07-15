from __future__ import division
import numpy
import scipy.interpolate
from zgoubi.core import *
from zgoubi.constants import *
from zgoubi.common import *
from zgoubi.utils import *
from zgoubi.exceptions import *

"""This module contains high level functions for analysing periodic cells and rings

Primary functions are
get_cell_properties()
get_cell_tracks()

"""

data_def = [
('KE',numpy.float64),
('stable',numpy.bool), # has closed orbit and stable tunes around it
('found_co',numpy.bool), # has closed orbit (may be unstable vertically)
('stable_tm_YT',numpy.bool), # has a stable transfer matrix in YT
('stable_tm_ZP',numpy.bool), # has a stable transfer matrix in ZP
('Y',numpy.float64),
('T',numpy.float64),
('Z',numpy.float64),
('P',numpy.float64),
('BETA_Y',numpy.float64),
('BETA_Z',numpy.float64),
('ALPHA_Y',numpy.float64),
('ALPHA_Z',numpy.float64),
('GAMMA_Y',numpy.float64),
('GAMMA_Z',numpy.float64),
('DISP_Y',numpy.float64),
('DISP_Z',numpy.float64),
('DISP_PY',numpy.float64),
('DISP_PZ',numpy.float64),
('NU_Y',numpy.float64),
('NU_Z',numpy.float64),
('matrix', numpy.float64, (6,6)),
('matrix_trace_YT', numpy.float64),
('matrix_trace_ZP', numpy.float64),
('tof',numpy.float64),
('S',numpy.float64),
('MAX_BY', numpy.float64),
('MIN_BY', numpy.float64),
('MAX_BZ', numpy.float64),
('MIN_BZ', numpy.float64),
('twiss_profile', numpy.object),
('full_twiss_profile', numpy.object),
('DA', numpy.object),
('DA_angles', numpy.object),
('ftrack', numpy.object),
('ptrack', numpy.object),
]




def part_info(particle):
	"Look up a particle by name and return the zgoubi PARTICUL, mass and charge sign"
	if particle == "p":
		part_ob = PROTON()
		mass = PROTON_MASS
		charge_sign = 1
	elif particle == "e":
		part_ob = ELECTRON()
		mass = ELECTRON_MASS
		charge_sign = -1
	elif particle == "mu-":
		part_ob = IMMORTAL_MUON()
		mass = MUON_MASS
		charge_sign = -1
	elif particle == "mu+":
		part_ob = IMMORTAL_MUON()
		mass = MUON_MASS
		charge_sign = +1
	elif particle == "Bismuth":
		part_ob = PARTICUL(M=ATOMIC_MASS_UNIT*238/1e6, Q=PROTON_CHARGE*4)
		mass = ATOMIC_MASS_UNIT*238
		charge_sign = 1
	else:
		raise ValueError("Unknown particle")
	return part_ob, mass, charge_sign



def get_cell_properties(cell, min_ke, max_ke=None, ke_steps=1, particle=None, tol=1e-6, stop_at_first_unstable=False, closed_orbit_range=None, closed_orbit_range_count=None, closed_orbit_init_YTZP=None, reuse_co_coords=True, closed_orbit_debug=False, full_tracking=False):
	"""Get the closed orbits and basic properties of a periodic cell.

	cell: A PyZgoubi Line object containing the beamline elements
	min_ke, max_ke, ke_steps: kinetic energy in eV. For a single step just set min_ke
	particle: "p", "e", "mu-", "mu+", or "Bismuth"
	tol: tolerance when finding closed orbit
	stop_at_first_unstable: If an energy is found to be unstable, give up
	closed_orbit_range, closed_orbit_range_count, closed_orbit_init_YTZ: see zgoubi.utils.find_closed_orbit_range()
	reuse_co_coords: use closed orbit from previous energy to start search for next energy
	closed_orbit_debug: output debuging information
	full_tracking=True is required in order get minimum and maximum magnetic fields along orbit

	returns orbit_data, an array with ke_steps elements, with the following data
	KE : particle KE in eV
	stable, found_co, stable_tm_YT, stable_tm_ZP : boolean stability flags
	Y, T, Z, P : closed orbit at start of cell, in cm and mrad
	BETA_Y, BETA_Z, ALPHA_Y, ALPHA_Z, GAMMA_Y, GAMMA_Z, DISP_Y, DISP_Z, DISP_PY, DISP_PZ, NU_Y, NU_Z: periodic twiss functions
	tof, S: time of flight and path length along closed orbit
	matrix, matrix_trace_YT, matrix_trace_ZP: transfer matrix and traces
	MAX_BY, MIN_BY, MAX_BZ, MIN_BZ: minimum and maximum fields seen along closed orbit (requires full_tracking=True)
	"""

	if max_ke==None: max_ke = min_ke
	ke_list = numpy.linspace(min_ke, max_ke, ke_steps)
	orbit_data =  numpy.zeros(ke_steps, data_def)
	for n, particle_ke in enumerate(ke_list):
		orbit_data['KE'][n] = particle_ke

	# get closed orbits
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)

	tline.add(DRIFT("fco", XL=0* cm_))

	tline.add(cell)

	tline.add(DRIFT("end", XL=0* cm_))
	tline.add(FAISCNL("end", FNAME='zgoubi.fai'))
	tline.add(REBELOTE(NPASS=20, K=99))
	tline.add(END())

	search_coords = [0,0,0,0]
	if closed_orbit_init_YTZP:
		search_coords = closed_orbit_init_YTZP
	if closed_orbit_range is not None and closed_orbit_range_count is None:
		# if range but no count, then do 10 steps in dimensions that are used	
		closed_orbit_range_count = []
		for r in closed_orbit_range:
			if r == 0: closed_orbit_range_count.append(0)
			else: closed_orbit_range_count.append(10)

	for n, particle_ke in enumerate(ke_list):
		orbit_data['found_co'][n] = False
		orbit_data['stable_tm_YT'][n] = False
		orbit_data['stable_tm_ZP'][n] = False
		orbit_data['stable'][n] = False
		print "closed orbit, energy = ", particle_ke
		rigidity = ke_to_rigidity(particle_ke,mass) * charge_sign
		ob.set(BORO=rigidity)

		if closed_orbit_debug:
			record_fname = "closedorbit_%s.log"%n
		else:
			record_fname = None

		if closed_orbit_range is None:
			closed_orbit =  find_closed_orbit(tline, init_YTZP=search_coords, tol=tol, max_iterations=50, record_fname=record_fname)
		else:
			closed_orbit =  find_closed_orbit_range(tline, range_YTZP=closed_orbit_range, count_YTZP=closed_orbit_range_count, init_YTZP=search_coords, tol=tol, max_iterations=50, record_fname=record_fname)
	
		search_coords = closed_orbit
		if closed_orbit != None:
			orbit_data['Y'][n],orbit_data['T'][n],orbit_data['Z'][n], orbit_data['P'][n]= closed_orbit
			orbit_data['found_co'][n] = True
		else:
			zlog.warn("No closed orbit at: %s"% particle_ke)
			if closed_orbit_debug:
				plot_find_closed_orbit(data_fname=record_fname, outfile=record_fname+".pdf")
				print "Search plots writen to ", record_fname+".pdf"
			if stop_at_first_unstable: break

	# get tunes and twiss
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET5()
	tline.add(ob)
	tline.add(part_ob)
	tline.add(DRIFT("fco", XL=0* cm_))
	tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(cell)

	tline.add(DRIFT("end", XL=0* cm_))
	tline.add(FAISCNL("end",FNAME='zgoubi.fai',))
	tline.add(MATRIX(IORD=1, IFOC=11))
	tline.add(END())

	for n, particle_ke in enumerate(ke_list):
		print "twiss, energy = ", particle_ke
		
		rigidity = ke_to_rigidity(particle_ke,mass) * charge_sign
		ob.set(BORO=rigidity)
		
		ref_Y,ref_T,ref_Z,ref_P = orbit_data['Y'][n],orbit_data['T'][n],orbit_data['Z'][n], orbit_data['P'][n]
		#print "ref_Y,ref_T,ref_Z,ref_P", ref_Y,ref_T,ref_Z,ref_P
		ob.set(YR=ref_Y, TR=ref_T, ZR=ref_Z, PR=ref_P, XR=0, DR=1)
		step_disp = 0.001 #cm
		step_ang =  0.01 #mrad
		ob.set(PY=step_disp, PT=step_ang, PZ=step_disp, PP=step_ang, PX=step_disp, PD=0.001)
		tline.full_tracking(full_tracking)

		res =tline.run()
		try:
			orbit_data['matrix'][n] = res.get_transfer_matrix()
		except BadLineError:
			continue
		
		orbit_data['matrix_trace_YT'][n] = orbit_data['matrix'][n][0,0] + orbit_data['matrix'][n][1,1]
		orbit_data['matrix_trace_ZP'][n] = orbit_data['matrix'][n][2,2] + orbit_data['matrix'][n][3,3]
		
		orbit_data['stable_tm_YT'][n] = abs(orbit_data['matrix_trace_YT'][n]) < 2
		orbit_data['stable_tm_ZP'][n] = abs(orbit_data['matrix_trace_ZP'][n]) < 2

		if orbit_data['stable_tm_YT'][n] and orbit_data['stable_tm_ZP'][n] and orbit_data['found_co'][n]:
			orbit_data['stable'][n] = True
		
		tune = res.get_tune()
		orbit_data['NU_Y'][n],orbit_data['NU_Z'][n] = tune
		
		if tune[0] == 0 or tune[1] == 0:
			zlog.warn("tune[0] == 0 or tune[1] == 0")
			#tline.run(xterm=1)
		
		twiss = res.get_twiss_parameters()
		orbit_data['BETA_Y'][n] = twiss['beta_y']
		orbit_data['ALPHA_Y'][n] = twiss['alpha_y']
		orbit_data['GAMMA_Y'][n] = twiss['gamma_y']
		orbit_data['DISP_Y'][n] = twiss['disp_y']
		orbit_data['DISP_PY'][n] = twiss['disp_py']
		orbit_data['BETA_Z'][n] = twiss['beta_z']
		orbit_data['ALPHA_Z'][n] = twiss['alpha_z']
		orbit_data['GAMMA_Z'][n] = twiss['gamma_z']
		orbit_data['DISP_Z'][n] = twiss['disp_z']
		orbit_data['DISP_PZ'][n] = twiss['disp_pz']

		if not orbit_data['found_co'][n]: continue
		if full_tracking:
			ptrack = res.get_all('plt')
			by = ptrack['BY']
			bz = ptrack['BZ']
			orbit_data['MAX_BY'][n] = by.max()
			orbit_data['MIN_BY'][n] = by.min()
			orbit_data['MAX_BZ'][n] = bz.max()
			orbit_data['MIN_BZ'][n] = bz.min()

		ftrack = res.get_all('fai')
		ftrack = ftrack[ftrack['ID']==1]
		orbit_data['tof'][n] = ftrack['tof'][-1]
		orbit_data['S'][n] = ftrack['S'][-1]
		res.clean()

		if orbit_data['NU_Y'][n] == -1 or orbit_data['NU_Z'][n] == -1:
			orbit_data['stable'][n] = False
			if stop_at_first_unstable: break
		
	return orbit_data


def get_cell_tracks(cell, data, particle, full_tracking=False):
	"""Get tracks along the closed orbit for values in data
	cell: periodic cell
	data: the data structure return from get_cell_properties()
	particle: see get_cell_properties()
	full_tracking: record all steps in the magnet

	tracks are added in the data structures ftrack and ptrack fields
	
	"""
	for n, particle_ke in enumerate(data['KE']):
		ref_Y,ref_T,ref_Z,ref_P = data[n]['Y'], data[n]['T'], data[n]['Z'], data[n]['P']
		tracks = get_tracks(cell=cell, start_YTZP=[ref_Y,ref_T,ref_Z,ref_P],
		                    particle=particle, ke=particle_ke, full_tracking=full_tracking)
		data[n]['ftrack'] = tracks['ftrack']
		data[n]['ptrack'] = tracks['ptrack']

	

def get_tracks(cell, start_YTZP, particle, ke, full_tracking=False, return_zgoubi_files=False, xterm=False):
	"""Run a particle through a cell from a give starting point and return track from fai and plt files.
	
	This is mostly used by other functions in this module, but can be useful for debugging a lattice
	
	"""
	split = 1
	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	
	tline.add(DRIFT('start', XL=0))

	for e in cell.elements():
		if split == 1:
			tline.add(e)
			if hasattr(e, "XL") or hasattr(e, "AT"):
				tline.add(FAISCNL(FNAME='zgoubi.fai',))
		else:
			for es in split_element(e, split):
				tline.add(es)
				if hasattr(e, "XL") or hasattr(e, "AT"):
					tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(DRIFT('end', XL=0))
	tline.add(END())

	rigidity = ke_to_rigidity(ke,mass) * charge_sign
	ob.set(BORO=rigidity)

	ref_Y,ref_T,ref_Z,ref_P = start_YTZP
	ob.clear()
	ob.add(Y=ref_Y, T=ref_T, Z=ref_Z, P=ref_P, X=0, D=1)
	if full_tracking:
		tline.full_tracking(True, drift_to_multi=True)

	res = tline.run(xterm=xterm)
	try:
		ftrack = res.get_all('fai')
	except IOError:
		ftrack = None
	if full_tracking:
		try:
			ptrack = res.get_all('plt')
		except EmptyFileError:
			zlog.warn("Empty plt file")
			ptrack = None
	else:
		ptrack = None
	ret_data = dict(ftrack=ftrack, ptrack=ptrack)

	if return_zgoubi_files:
		ret_data['dat_file'] = res.dat()
		ret_data['res_file'] = res.res()

	return ret_data



def plot_cell_properties(data, output_prefix="results/cell_", file_fmt=".pdf", ncells=1):
	"""Produce a set of standard plots from the data structure
	
	data: the data structure returned by get_cell_properties()

	"""
	import os
	from matplotlib import pyplot
	pyplot.clf()

	mkdir_p(os.path.dirname(output_prefix))

	#data = numpy.load(data_file_path)
	n_vars = len(data.dtype.names)

	# drop rows that dont start with zero
	stable_data =  data[data['stable']]

	n_energies =  stable_data.shape[0]
	if n_energies == 0:
		zlog.warn("No stable properties to plot")
		return

	off_x, off_y = 0.1, -3
	# horizontal
	for n in xrange(n_energies):
		pyplot.plot([stable_data['Y'][n]], [stable_data['T'][n]], 'bx')
		if (n_energies < 5) or ((n%int(n_energies/5))==0):
			pyplot.annotate("%.3f"%(stable_data['KE'][n]/1e6), xy=(stable_data['Y'][n], stable_data['T'][n]),
			               xycoords='data', xytext=(10, 10), textcoords='offset points',
						   arrowprops=dict(arrowstyle="->", color="grey"))
	pyplot.title("Horizontal closed orbit")
	pyplot.xlabel("Y (cm)")
	pyplot.ylabel("T (mrad)")
	if (stable_data['T'].max() - stable_data['T'].min()) < 1e-3:
		y_min = stable_data['T'].min()
		y_max = stable_data['T'].max()
		pyplot.ylim([ y_min-(y_max-y_min)*1 , y_max+(y_max-y_min)*1 ])
	pyplot.tight_layout()
	pyplot.savefig("%s_hor"%(output_prefix)+file_fmt)
	label = "Min = %s, Max = %s" % (stable_data['Y'].min(), stable_data['Y'].max())

	# vert
	pyplot.clf()
	off_x, off_y = 0, 0
	for n in xrange(n_energies):
		pyplot.plot([stable_data['Z'][n]], [stable_data['P'][n]], 'bx')
		if (n_energies < 5) or ((n%int(n_energies/5))==0):
			pyplot.annotate("%.3f"%(stable_data['KE'][n]/1e6), xy=(stable_data['Z'][n], stable_data['P'][n]),
			               xycoords='data', xytext=(10, 10), textcoords='offset points',
						   arrowprops=dict(arrowstyle="->", color="grey"))
	pyplot.title("Vertical closed orbit")
	pyplot.xlabel("Z (cm)")
	pyplot.ylabel("P (mrad)")
	pyplot.tight_layout()
	pyplot.savefig("%s_ver"%(output_prefix)+file_fmt)
	label = "Min = %s, Max = %s" % (stable_data['Z'].min(), stable_data['Z'].max())

	#tune
	pyplot.clf()
	#pyplot.title("Cell Tune")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell Tune")

	pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Y'], "b-x", label='Horizontal')
	pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Z'], "g-x", label='Vertical')
	pyplot.legend(loc='best')

	pyplot.tight_layout()
	pyplot.savefig("%s_tune"%(output_prefix)+file_fmt)

	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Ring Tune")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring Tune")

		pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Y']*ncells, "b-x", label='Horizontal')
		pyplot.plot(stable_data['KE']/1e6,stable_data['NU_Z']*ncells, "g-x", label='Vertical')
		pyplot.legend(loc='best')

		pyplot.savefig("%s_ring_tune"%(output_prefix)+file_fmt)

	#tof
	pyplot.clf()
	#pyplot.title("Cell time of flight")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell Time of flight ($\mu$s)")
	pyplot.plot(stable_data['KE']/1e6,stable_data['tof'], "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_tof"%(output_prefix)+file_fmt)

	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Cell time of flight")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring Time of flight ($\mu$s)")
		pyplot.plot(stable_data['KE']/1e6, stable_data['tof']*ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_tof"%(output_prefix)+file_fmt)
		label = "Mean time of flight = %s microseconds" % (stable_data['tof']*ncells).mean()

	#freq
	pyplot.clf()
	#pyplot.title("Revolution Frequency")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell revolution frequence (MHz)")
	pyplot.plot(stable_data['KE']/1e6,1/(stable_data['tof']), "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_freq"%(output_prefix)+file_fmt)
	
	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Revolution Frequency")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Revolution frequence (MHz)")
		pyplot.plot(stable_data['KE']/1e6,1/(stable_data['tof'])/ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_freq"%(output_prefix)+file_fmt)
		label = "Min = %s, Max = %s" % ((1/(stable_data['tof'])/ncells).min(), (1/(stable_data['tof'])/ncells).max())

	#path length
	pyplot.clf()
	#pyplot.title("Path length")
	pyplot.xlabel("KE (MeV)")
	pyplot.ylabel("Cell path length")
	pyplot.plot(stable_data['KE']/1e6,(stable_data['S']), "bx")
	pyplot.tight_layout()
	pyplot.savefig("%s_pathlen"%(output_prefix)+file_fmt)
	
	if(ncells>1):
		pyplot.clf()
		#pyplot.title("Path length")
		pyplot.xlabel("KE (MeV)")
		pyplot.ylabel("Ring path length")
		pyplot.plot(stable_data['KE']/1e6,(stable_data['S'])*ncells, "bx")
		pyplot.tight_layout()
		pyplot.savefig("%s_ring_pathlen"%(output_prefix)+file_fmt)
		label = "Min = %s, Max = %s" % (((stable_data['S'])*ncells).min(), ((stable_data['S'])*ncells).max())
	page_images = []


def cell_properties_table(data, keys, sep="\t"):
	"""Quick function for printing parameter tables, eg:
	
	    print gcp.cell_properties_table(data_arc, ["KE", "stable", "Y", "T", "NU_Y", "NU_Z"])

	"""
	out = ""
	out += sep.join(keys) + "\n"
	for d in data:
		row = [str(d[k]) for k in keys]
		out += sep.join(row) + "\n"
	return out



def get_twiss_params(cell, data, particle, output_prefix="results/twiss_profiles_", full_tracking=False):
	"""Get the periodic Twiss (Courant and Snyder) parameters for the cell. Tracks a bunch of particles containing pairs offset in each plane.

	Must pass in data returned from get_cell_properties() to give initial conditions. The profile is added to the 'twiss_profile' key in the data structure. If full_tracking is enabled, the Twiss parameters are every integration step are recorded, otherwise parameters are recorded at the end of each magnetic element.


	"""
	mkdir_p(os.path.dirname(output_prefix))

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET5(DR=1, PY=1e-4,PT=1e-3,PZ=1e-4,PP=1e-3,PX=1e-3,PD=1e-3)
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	
	tline.add(DRIFT('start', XL=0))

	for e in cell.elements():
		isphysical = False
		try:
			dummy = e.XL
			isphysical = True
		except AttributeError:pass
		try:
			dummy = e.AT
			isphysical = True
		except AttributeError:pass
		tline.add(e)
		if isphysical:
			tline.add(FAISCNL(FNAME='zgoubi.fai',))

	tline.add(DRIFT('end', XL=0))
	tline.add(MATRIX(IORD=1,IFOC=11))	
	tline.add(END())

	for n, particle_ke in enumerate(data['KE']):
		if not data[n]['stable']: continue
		print "energy = ", particle_ke
		rigidity = ke_to_rigidity(particle_ke,mass) * charge_sign
		ob.set(BORO=rigidity)

		ref_Y,ref_T,ref_Z,ref_P = data[n]['Y'], data[n]['T'], data[n]['Z'], data[n]['P']
		ob.set(YR=ref_Y, TR=ref_T, ZR=0, PR=0, XR=0)
		init_twiss = numpy.zeros(1,dtype=[('beta_y','f8'),('alpha_y','f8'),('gamma_y','f8'),('disp_y','f8'),('disp_py','f8'),
		                         ('beta_z','f8'),('alpha_z','f8'),('gamma_z','f8'),('disp_z','f8'),('disp_pz','f8')])
		init_twiss["beta_y"] = data[n]["BETA_Y"]
		init_twiss['alpha_y'] = data[n]['ALPHA_Y']
		init_twiss['gamma_y'] = data[n]['GAMMA_Y']
		init_twiss['disp_y'] = data[n]['DISP_Y']
		init_twiss['disp_py'] = data[n]['DISP_PY']
		init_twiss['beta_z'] = data[n]['BETA_Z']
		init_twiss['alpha_z'] = data[n]['ALPHA_Z']
		init_twiss['gamma_z'] = data[n]['GAMMA_Z']
		init_twiss['disp_z'] = data[n]['DISP_Z']
		init_twiss['disp_pz'] = data[n]['DISP_PZ']
		#data[n][["BETA_Y", "ALPHA_Y", "GAMMA_Y", "DISP_Y", "DISP_PY", "BETA_Z", "ALPHA_Z", "GAMMA_Z", "DISP_Z", "DISP_PZ"]]
	
		tline.full_tracking(False)
		twiss_profiles = get_twiss_profiles(tline,'%s%s.txt'%(output_prefix, particle_ke), input_twiss_parameters=init_twiss, calc_dispersion=0)
		data[n]['twiss_profile'] = twiss_profiles
		if full_tracking:
			tline.full_tracking(True)
			twiss_profiles = get_twiss_profiles(tline, '%s%s_full.txt'%(output_prefix, particle_ke), calc_dispersion=False, input_twiss_parameters=init_twiss)
			data[n]['full_twiss_profile'] = twiss_profiles


def plot_twiss_params(data, output_prefix="results/twiss_profiles_"):
	"""Plot the Twiss profiles found with get_twiss_params()

	"""
	import pylab
	stable_data =  numpy.extract(data['stable'], data)
	full_tracking_message = 0
	for n, particle_ke in enumerate(stable_data['KE']):
		if stable_data[n]['full_twiss_profile'] != 0:
			twiss_profiles = stable_data[n]['full_twiss_profile']
		elif stable_data[n]['twiss_profile'] != 0:
			if full_tracking_message == 0:
				zlog.warn("get_twiss_params() called without full_tracking=True, so only plotting twiss at element ends")
				full_tracking_message = 1
			twiss_profiles = stable_data[n]['twiss_profile']
		else:
			zlog.error("Call get_twiss_params() before plot_twiss_params()")
			raise ValueError

		pylab.clf()
		ax1 = pylab.axes()
		ax2 = ax1.twinx()
		l2 = ax1.plot(twiss_profiles['s'],twiss_profiles['beta_y'],"-b", label=r"$\beta_y$")
		l4 = ax1.plot(twiss_profiles['s'],twiss_profiles['beta_z'],"-c", label=r"$\beta_z$")
		l6 = ax2.plot(twiss_profiles['s'],twiss_profiles['alpha_y'],"-r", label=r"$\alpha_y$")
		l8 = ax2.plot(twiss_profiles['s'],twiss_profiles['alpha_z'],"-m", label=r"$\alpha_z$")
		lns = l2+l4+l6+l8
		ax1.set_xlabel("Path length (m)")
		ax1.set_ylabel(r"$\beta$ (m)", color='k')
		ax2.set_ylabel(r"$\alpha$", color='k')
		pylab.legend(lns, [l.get_label() for l in lns],loc="best") # stackoverflow.com/questions/5484922
		pylab.savefig('%s%s.pdf'%(output_prefix, particle_ke))




def plot_cell_tracks(cell, data, particle, output_file="results/track.pdf", show=False, plot_unstable=False, draw_field_midplane=False, sector_width=None, aspect="equal", style=None, draw_tracks=True, min_y=None, max_y=None, y_steps=None, angle=0):
	"""Plot particle track through cell, starting from in the closed orbits stored in data.

	output_file: file to save plot
	show: display plot to screen
	plot_unstable: show tracks for initial coordinates marked as unstable in data structure
	draw_field_midplane: show the magnet field on the magnet midplane
	sector_width: width to draw sector magnets
	aspect: "equal", use same scale for x and y, "auto" stretch to fit
	style: style for plotting, see lab_plot.LabPlot
	draw_tracks: draw particle tracks
	min_y, max_y, y_steps, angle: starting coordinates for test particles for sampling the midplane field

	"""
	from zgoubi.lab_plot import LabPlot
	import pylab
	pylab.close()
	n_vars = len(data.dtype.names)

	if plot_unstable==False:
		# drop rows that dont start with zero
		stable_data =  numpy.extract(data['stable'], data)
	else:
		stable_data = data


	cell = uniquify_labels(cell)

	tline = Line('test_line')
	tline.add_input_files(cell.input_files)
	ob = OBJET2()
	tline.add(ob)
	part_ob, mass, charge_sign = part_info(particle)
	tline.add(part_ob)
	tline.add(DRIFT("gcpstart", XL=0* cm_))
	#add(FAISCEAU("fco"))
	tline.add(FAISCNL("gcpstart",FNAME='zgoubi.fai',))
	tline.add(cell)
	tline.add(DRIFT("gcpend", XL=0* cm_))
	tline.add(FAISCNL("gcpend", FNAME='zgoubi.fai'))
	#add(REBELOTE(NPASS=9, K=99))
	tline.add(END())
	tline.full_tracking(True, drift_to_multi=True)

	# if line contains BENDS, then line will be drawn with first energy, but zgoubi will adjust angles for other particles
	if len(data['KE']>0):
		boro = ke_to_rigidity(data['KE'][0],mass) * charge_sign
	else:
		boro = None
	lp = LabPlot(tline, boro=boro, sector_width=sector_width, aspect=aspect)
	if style:
		lp.set_style(style)
	

	if draw_tracks:
		for n, particle_ke in enumerate(stable_data['KE']):
			print "energy = ", particle_ke
			rigidity = ke_to_rigidity(particle_ke,mass) * charge_sign
			ob.set(BORO=rigidity)

			ref_Y,ref_T,ref_Z,ref_P = stable_data[n]['Y'], stable_data[n]['T'],stable_data[n]['Z'],stable_data[n]['P']
			ob.clear()
			ob.add(Y=ref_Y, T=ref_T, Z=0, P=0, X=0, D=1)

			res = tline.run()
			try:
				ftrack = res.get_all('fai')
				ptrack = res.get_all('plt')
			except IOError:
				zlog.error("Failed to read fai or plt files")
				print res.res()
				raise
			lp.add_tracks(ftrack=ftrack, ptrack=ptrack, draw=1)
	
	if draw_field_midplane:
		if min_y is None or max_y is None or y_steps is None:
			raise ValueError("When using draw_field_midplane, you must set min_y, max_y and y_steps")

		for Y in numpy.linspace(min_y, max_y, y_steps):
			rigidity = ke_to_rigidity(1e15, mass) * charge_sign
			ob.set(BORO=rigidity)

			ref_Y,ref_T,ref_Z,ref_P = Y, angle, 0, 0
			ob.clear()
			ob.add(Y=ref_Y, T=ref_T, Z=0, P=0, X=0, D=1)

			res = tline.run()
			try:
				zlog.error("Failed to read fai or plt files")
				ftrack = res.get_all('fai')
				ptrack = res.get_all('plt')
			except IOError:
				print res.res()
				raise
			lp.add_tracks(ftrack=ftrack, ptrack=ptrack, draw=0, field=1)


	if draw_field_midplane:
		lp.draw(draw_field_midplane=draw_field_midplane, draw_tracks=draw_tracks, field_steps=y_steps)
	else:
		lp.draw()

	lp.save(output_file)
	if show:
		lp.show()






