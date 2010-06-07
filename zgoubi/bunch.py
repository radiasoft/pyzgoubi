#!/usr/bin/env python
"""A bunch object to hold the coordinates for many particles

Note that all values are in SI units, m, rad, eV, s
"""

from __future__ import division
from math import *
import numpy
import pylab
from zgoubi import rel_conv
from zgoubi import io
import struct

#from zgoubi.utils import *


class Bunch(object):
	"""Object to store a bunch of particles efficiently using numpy.
	All values are in SI units, m, rad, eV, s
	Can be use to create a bunch of particles with all there coordinates set to zero (appart from the D coordinate which is set to 1)::
	
		my_bunch = Bunch(nparticles=100, ke=1e6, mass=PROTON_MASS, charge=1)
	
	It can also be used to create a bunch from some existing coordinates stored in a numpy array::
	
		 my_bunch = Bunch(ke=1e6, mass=PROTON_MASS, charge=1, particles=existing_coords)
	
	"""
	min_data_def = [
	('D', numpy.float64), # these coorspond to zgoubi D,Y,T,Z,P,S, but in SI units
	('Y', numpy.float64),
	('T', numpy.float64),
	('Z', numpy.float64),
	('P', numpy.float64),
	('S', numpy.float64),
	('tof', numpy.float64), # these are for accumulating 
	('X', numpy.float64),
	]

	def __init__(self, nparticles=0, ke=None, rigidity=0, mass=0, charge=1, particles=None):
		"""Bunch constructor.

		"""
		if particles != None:
			self.coords = particles
		else:
			self.coords = numpy.zeros(nparticles, self.min_data_def)
			self.coords['D'] = 1
		self.mass = mass
		self.charge = charge
		self.rigidity = rigidity 
		if ke != None:
			self.set_bunch_ke(ke)

	def split_bunch(self, max_particles, n_slices):
		"Split a bunch into n_slices smaller bunches, or more if they would have too many particles in."
		if ceil(len(self.coords) / n_slices) > max_particles:
			n_slices = ceil(len(self.coords)/max_particles)

		for pslice in numpy.array_split(self.coords, n_slices):
			yield Bunch(rigidity=self.get_bunch_rigidity(), mass=self.mass, charge=self.charge,
			            particles=pslice)

	def __str__(self):
		out = "Bunch:\n"
		out += "\t"+str(len(self)) + " paricles\n"
		out += "\t"+str(self.get_bunch_ke()) + " eV\n"

		return out

	def set_bunch_ke(self, ke):
		"Set bunch kinetic energy"
		if self.mass == 0:
			raise ValueError, "Particle mass can't be Zero"
		if self.charge == 0:
			raise ValueError, "Particle charge can't be Zero"
		self.rigidity = rel_conv.ke_to_rigidity(mass=self.mass, ke=ke, charge=self.charge)

	def get_bunch_ke(self):
		"Get bunch kinetic energy"
		if self.mass == 0:
			raise ValueError, "Particle mass can't be Zero"
		if self.charge == 0:
			raise ValueError, "Particle charge can't be Zero"
		return rel_conv.rigidity_to_ke(mass=self.mass, rigidity=self.rigidity, charge=self.charge)

	def set_bunch_rigidity(self, rigidity):
		"Set bunch rigidity"
		self.rigidity = rigidity

	def get_bunch_rigidity(self):
		"Get bunch rigidity"
		return self.rigidity

	def particles(self):
		"Returns the numpy array that holds the coordinates"
		return self.coords

	def get_min_BORO(self):
		"Returns the minimum rigidity of the bunch"
		#min_BORO = ke_to_rigidity(min(self.coords['KE']), self.mass)
		min_BORO = self.rigidity * self.coords['D'].min()
		return min_BORO

	@staticmethod
	def gen_halo_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a halo bunch, i.e. an elipse outline in x-xp (Y-T) and in y-yp (Z-P). S and D are set to 0 and 1 respectively.
		example::
			my_bunch = Bunch.gen_halo_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a halo bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print "Emittance or beta can't be negative"
			print "emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z"
			print emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z
			raise ValueError

		ry = sqrt(emit_y) 
		rz = sqrt(emit_z) 

		if seed != None:
			numpy.random.seed(seed)

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ry * numpy.cos(u1)
		coords[:, 1] = ry * numpy.sin(u1)
		coords[:, 2] = rz * numpy.cos(u2)
		coords[:, 3] = rz * numpy.sin(u2)

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch =  numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_kv_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a uniform (KV) bunch, i.e. a filled elipse in x-xp (Y-T) and in y-yp (Z-P). S and D are set to 0 and 1 respectively.
		example::
			my_bunch = Bunch.gen_kv_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a KV bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print "Emittance or beta can't be negative"
			print "emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z"
			print emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z
			raise ValueError

		if seed != None:
			numpy.random.seed(seed)

		ry = sqrt(emit_y) * numpy.random.random_sample([npart])
		rz = sqrt(emit_z) * numpy.random.random_sample([npart]) 

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ry * numpy.cos(u1)
		coords[:, 1] = ry * numpy.sin(u1)
		coords[:, 2] = rz * numpy.cos(u2)
		coords[:, 3] = rz * numpy.sin(u2)

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch =  numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_gauss_x_xp_y_yp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None,
			               ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a Gaussian bunch in x-xp (Y-T) and in y-yp (Z-P). S and D are set to 0 and 1 respectively.
		example::
			my_bunch = Bunch.gen_kv_x_xp_y_yp(1000, 1e-3, 1e-3, 4, 5, 1e-3, 2e-2, ke=50e6, mass=PROTON_MASS, charge=1)
		
		creates a Gaussian bunch called my_bunch with 1000 particles of the given parameters.

		"""

		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print "Emittance or beta can't be negative"
			print "emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z"
			print emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z
			raise ValueError

		if seed != None:
			numpy.random.seed(seed)

		ry = sqrt(emit_y) * numpy.random.normal(0, 0.5, [npart])
		rz = sqrt(emit_z) * numpy.random.normal(0, 0.5, [npart]) 

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ry * numpy.cos(u1)
		coords[:, 1] = ry * numpy.sin(u1)
		coords[:, 2] = rz * numpy.cos(u2)
		coords[:, 3] = rz * numpy.sin(u2)

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)
		
		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)
		
		bunch =  numpy.zeros([npart], Bunch.min_data_def)

		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		bunch['D'] = 1
	
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def gen_gauss_x_xp_y_yp_s_dp(npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, mom_spread=0, bunch_length=0, disp=0, disp_prime=0, seed=None, ke=None, rigidity=0, mass=0, charge=1):
		"""Generate a Gaussian bunch in transverse and longitudinal phase space
		emit_y, emit_z : horizontal and vertical plane geometric emittance (1 sigma)
		beta_y, beta_z : horizontal and vertical betatron function
		alpha_y, alpha_z : horizontal and vertical alpha
		mom_spead : sigma of momentum spread (percentage) 
		bunch_length : sigma of bunch length (metres)
		disp, disp_prime : dispersion and dispersion prime (D') in horizontal plane """


		if emit_y < 0 or emit_z < 0 or beta_y < 0 or beta_z < 0:
			print "Emittance or beta can't be negative"
			print "emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z"
			print emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z
			raise ValueError


		if seed != None:
			numpy.random.seed(seed)

		#generate momentum distribution
                if mom_spread > 0.0:
                        mom_dist = numpy.random.normal(1.0, mom_spread/100, npart)
		#generate longitudinal coordinate distribution (bunch length)
		if bunch_length > 0.0:
			s_dist = numpy.random.normal(0.0, bunch_length, npart)


		#use numpy.random.multivariate_normal. 
		#set off-diagonal terms in covariance matrix to zero so that Y is uncorrelated with T (and Z with P)
		cov_yt = [[emit_y,0],[0,emit_y]]
		ryrt = numpy.random.multivariate_normal((0,0),cov_yt,(1,npart))[0]
		cov_zp = [[emit_z,0],[0,emit_z]]
		rzrp = numpy.random.multivariate_normal((0,0),cov_zp,(1,npart))[0]

		coords = numpy.zeros([npart, 6], numpy.float64)
		#coords2 = numpy.zeros([npart, 6], numpy.float64)

		coords[:, 0] = ryrt[:,0]
		coords[:, 1] = ryrt[:,1]
		coords[:, 2] = rzrp[:,0]
		coords[:, 3] = rzrp[:,1]

		matrix = Bunch._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z)

		for n, coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			coords[n] = numpy.dot(matrix, coord)

		#Add dispersion term to horizontal coord
		if disp != 0.0:
			coords[:, 0] = [y+disp*(dp-1) for y,dp in zip(coords[:, 0],mom_dist)] 

		#Add dispersion prime term to horizontal angle
		if disp_prime != 0.0:
			coords[:, 1] = [t+disp_prime*(dp-1) for t,dp in zip(coords[:, 1],mom_dist)] 
		
		bunch =  numpy.zeros([npart], Bunch.min_data_def)
		
		bunch['Y'] = coords[:, 0]
		bunch['T'] = coords[:, 1]
		bunch['Z'] = coords[:, 2]
		bunch['P'] = coords[:, 3]
		if mom_spread > 0.0:
			bunch['D'] = mom_dist
		else:
			bunch['D'] = 1

		if bunch_length > 0.0:
			bunch['S'] = s_dist
		else:
			bunch['S'] = 0.0
		
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=bunch)

	@staticmethod
	def _twiss_matrix(beta_y, beta_z, alpha_y, alpha_z):
		"Create a matrix that will convert a spherical distribution in to one with the correct twiss values"
		B = numpy.eye(6)
		B[0, 0] = sqrt(beta_y)
		B[1, 1] = 1 / sqrt(beta_y)
		B[2, 2] = sqrt(beta_z)
		B[3, 3] = 1 / sqrt(beta_z)

		A = numpy.eye(6)
		A[1, 0] = -alpha_y
		A[3, 2] = -alpha_z

		#M = numpy.dot(A, B)
		M = numpy.dot(B, A)
		return M


	@staticmethod
	def read_YTZPSD(fname, ke=None, rigidity=0, mass=0, charge=1):
		"""Read in a bunch from a file. assumes columns are Y, T, Z, P, X, D, separated by white space::
			my_bunch = Bunch.read_YTZPSD("mybunch.dat", ke=1e9, mass=ELECTRON_MASS, charge=-1)
		
		"""
		dist = numpy.loadtxt(fname)
		nparts = dist.size / 6
		dist = dist.reshape(nparts, 6)
		coords = numpy.zeros(nparts, Bunch.min_data_def)
		#self.coords['KE'] = ke
		coords['Y'] = dist[:, 0]
		coords['T'] = dist[:, 1]
		coords['Z'] = dist[:, 2]
		coords['P'] = dist[:, 3]
		coords['X'] = dist[:, 4]
		coords['D'] = dist[:, 5]
		return Bunch(ke=ke, rigidity=rigidity, mass=mass, charge=charge, particles=coords)

	def write_YTZPSD(self, fname, binary=False):
		"Output a bunch, compatible with read_YTZPSD"
		
		# it ought to be possible to do this:
		#numpy.savetxt(fh, self.coords[['Y','T','Z','P','S','D']])
		# but the field end up in the wrong order, see http://thread.gmane.org/gmane.comp.python.numeric.general/36933

		fh = open(fname, "w")
		if binary:
			#header
			for dummy in xrange(4):
				io.write_fortran_record(fh, "a"*80)
			#data
			for p in self.coords:
				record = struct.pack("6d", p['Y'], p['T'], p['Z'], p['P'], p['S'], p['D'])
				io.write_fortran_record(fh, record)

		else:
			nparts = len(self.coords)
			dist = numpy.zeros([nparts, 6])
			dist[:, 0] = self.coords['Y']
			dist[:, 1] = self.coords['T']
			dist[:, 2] = self.coords['Z']
			dist[:, 3] = self.coords['P']
			dist[:, 4] = self.coords['S']
			dist[:, 5] = self.coords['D']
			#dist = dist.reshape(nparts * 2, 3)
			fh.write("# bunch\n\n\n\n")
			numpy.savetxt(fh, dist)
		fh.close()

	
	def get_widths(self):
		"Returns the width of the bunch in each dimension Y,T,Z,P,S,D (D not calculated yet)"
		y_width = numpy.max(self.coords['Y']) - numpy.min(self.coords['Y'])
		t_width = numpy.max(self.coords['T']) - numpy.min(self.coords['T'])
		z_width = numpy.max(self.coords['Z']) - numpy.min(self.coords['Z'])
		p_width = numpy.max(self.coords['P']) - numpy.min(self.coords['P'])
		s_width = numpy.max(self.coords['S']) - numpy.min(self.coords['S'])
		d_width = numpy.max(self.coords['D']) - numpy.min(self.coords['D'])
		return (y_width, t_width, z_width, p_width, s_width, d_width)

	def get_centers(self):
		"Returns the center of the bunch in each dimension Y,T,Z,P,S,D"
		y_mean = numpy.mean(self.coords['Y'])
		t_mean = numpy.mean(self.coords['T'])
		z_mean = numpy.mean(self.coords['Z'])
		p_mean = numpy.mean(self.coords['P'])
		s_mean = numpy.mean(self.coords['S'])
		d_mean = numpy.mean(self.coords['D'])
		return (y_mean, t_mean, z_mean, p_mean, s_mean, d_mean)

	def __len__(self):
		"Returns length of bunch. Use len(my_bunch)"
		return len(self.coords)

	def plot(self, fname=None, lims=None, add_bunch=None):
		"Plot a bunch, if no file name give plot is displayed on screen. lims can be used to force axis limits eg [lY,lT,lZ,lP,lX,lD] would plot limit plot from -lY to +lY in Y, etc. Additional bunches can be passed, as add_bunch, to overlay onto the same plot."
		
		bunches = [self]
		if add_bunch != None:
			try:
				# if add_bunch is iterable, append all is memebers
				for a_bunch in add_bunch:
					bunches.append(a_bunch)
			except TypeError:
				# otherwise just append it
				bunches.append(add_bunch)

		pylab.subplot(2, 2, 1)
		pylab.grid()
		for abunch in bunches:
			pylab.plot(abunch.coords['Y'], abunch.coords['Z'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[2], lims[2])
		plotname = "X-Y (Y-T)"
		pylab.title(plotname)

		
		pylab.subplot(2, 2, 2)
		pylab.grid()
		for abunch in bunches:
			pylab.plot(abunch.coords['Y'], abunch.coords['T'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[1], lims[1])
		plotname = "X-XP (Y-T"
		pylab.title(plotname)

		pylab.subplot(2, 2, 3)
		pylab.grid()
		for abunch in bunches:
			pylab.plot(abunch.coords['Z'], abunch.coords['P'], ',')
		if lims != None:
			pylab.xlim(-lims[2], lims[2])
			pylab.ylim(-lims[3], lims[3])
		plotname = "Y-YP (z-P)"
		pylab.title(plotname)

		pylab.subplot(2, 2, 4)
		pylab.grid()
		for abunch in bunches:
			pylab.plot(abunch.coords['X'], abunch.coords['D'], ',')
		if False:# lims != None:
			pylab.xlim(-lims[4], lims[4])
			pylab.ylim(-lims[5], lims[5])
		plotname = "X-D"
		pylab.title(plotname)
		if fname == None:
			pylab.show()
		else:
			pylab.savefig(fname, dpi=300)
			pylab.clf()
	
	def get_emmitance(self):
		"return emittance h and v in m rad. Uses the bunch full width, so should only be used for a hard edge distribution"
		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]

		r_yt = numpy.sqrt(Ys**2 + Ts**2)
		theta_yt = numpy.arctan2(Ys, Ts)
		major_angle = theta_yt[r_yt.argmax()]
		theta_yt -= major_angle
		rot_y = r_yt * numpy.sin(theta_yt)
		rot_t = r_yt * numpy.cos(theta_yt)
		emmitance_h = rot_y.max()*rot_t.max()

		
		r_zp = numpy.sqrt(Zs**2 + Ps**2)
		theta_zp = numpy.arctan2(Zs, Ps)
		major_angle = theta_zp[r_zp.argmax()]
		theta_zp -= major_angle
		rot_z = r_zp * numpy.sin(theta_zp)
		rot_p = r_zp * numpy.cos(theta_zp)
		emmitance_v = rot_z.max() * rot_p.max()
		#print "Emmitance (h, v):", emmitance_h, emmitance_v
		return (emmitance_h, emmitance_v)


	def get_twiss(self, emittance):
		"Returns the twiss valuse Beta_h, Alpha_h, Beta_v, Alpha_v, calculated from bunch width extent"
		# emittance may be a single number, or tuple (emittance_h, emittance_v)
		try:
			emittance_h, emittance_v = emittance
		except TypeError:
			emittance_h = emittance_v = emittance
		widths = self.get_widths()
		centers = self.get_centers()
		Ys = self.coords['Y'] - centers[0] # work relative to center
		Ts = self.coords['T'] - centers[1]
		Zs = self.coords['Z'] - centers[2]
		Ps = self.coords['P'] - centers[3]
		beta_h = (widths[0] / 2)**2 / emittance_h
		beta_v = (widths[2] / 2)**2 / emittance_v
		#print "beta", beta_h, beta_v
		#gamma_h = (widths[1] / 2)**2 / emittance_h
		#gamma_v = (widths[3] / 2)**2 / emittance_v
		# get T of particle with bigest Y
		y_p_max = Ts[Ys.argmax()]
		y_p_min = Ts[Ys.argmin()]
		alpha_h = (y_p_min - y_p_max) / 2 / sqrt(emittance_h / beta_h)

		z_p_max = Ps[Zs.argmax()]
		z_p_min = Ps[Zs.argmin()]
		alpha_v = (z_p_min - z_p_max) / 2 / sqrt(emittance_v / beta_v)
		#print "gamma", gamma_h, gamma_v	

		# following is dangerous, due do stat fluctuations causing roots of negative numbers
		#	alpha_h = sqrt((gamma_h*beta_h)-1 )
		#	alpha_v = sqrt((gamma_v*beta_v)-1 )

		#print "twiss (bh,ah,bv,av):", beta_h, alpha_h, beta_v, alpha_v
		return beta_h, alpha_h, beta_v, alpha_v


