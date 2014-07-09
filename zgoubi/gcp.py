from __future__ import division
import numpy
import scipy.interpolate
from zgoubi.core import *
from zgoubi.constants import *
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
	Y, T, Z, P : closed orbit at start of cell
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
			print "No closed orbit at:", particle_ke
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


