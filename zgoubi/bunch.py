#!/usr/bin/env python
from math import *
import numpy
import pylab
from zgoubi import rel_conv
import os

from zgoubi.utils import *


class Bunch(object):
	"""Object to store a bunch of particles efficently using numpy.
	All values are in SI units, m, rad, eV, s

	"""
	data_def = [
	('KE', numpy.float64),
	('Y', numpy.float64),
	('T', numpy.float64),
	('Z', numpy.float64),
	('P', numpy.float64),
	('tof', numpy.float64),
	('S', numpy.float64),
	('X', numpy.float64),
	#('D', numpy.float64),
	('keep', numpy.bool)
	]
	def __init__(self, nparticles=0, ke=0, mass=0):
		self.coords = numpy.zeros(nparticles, self.data_def)
		self.coords['KE'] = ke
		self.coords['keep'] = True
		self.mass = mass

	def particles(self):
		return self.coords

	def get_min_BORO(self): #FIXME proton mass hardwired 
		min_BORO = ke_to_rigidity(min(self.coords['KE']), mass)
		# FIXME could be faster with ufunc ?
	#	self.coords['D'] = [ke_to_rigidity(p['KE'], PROTON_MASS)/min_BORO for p in self.coords]
		return min_BORO

	def drop_unkept(self):
		new_count = self.coords['keep'].sum()
		new_coords = numpy.zeros(new_count, self.data_def)
		i = 0
		for p in self.coords:
			if p['keep']:
				new_coords[i] = p
				i += 1
		self.coords = new_coords


	def gen_halo_x_xp_y_yp(self, npart, emit_y, emit_z, beta_y, beta_z, alpha_y, alpha_z, seed=None):

		#r = numpy.random.random_sample([npart])
		ry = sqrt(emit_y) 
		rz = sqrt(emit_z) 

		if seed!=None:
			numpy.random.seed(seed)

		u1 = numpy.random.random_sample([npart]) * pi * 2
		u2 = numpy.random.random_sample([npart]) * pi * 2
		#u3 = numpy.random.random_sample([npart]) * pi * 2
		#u4 = numpy.random.random_sample([npart]) * pi * 2

		coords = numpy.zeros([npart, 6], numpy.float64)
		coords2 = numpy.zeros([npart, 6], numpy.float64)

		coords[:,0] = ry * numpy.cos(u1)
		coords[:,1] = ry * numpy.sin(u1)
		coords[:,2] = rz * numpy.cos(u2)
		coords[:,3] = rz * numpy.sin(u2)

		matrix = self._twiss_matrix(beta_y, beta_z, alpha_y, alpha_z )
		
		for n,coord in enumerate(coords):
			#	new_coord = numpy.dot(coord, matrix)
			new_coord = numpy.dot(matrix, coord)
			coords2[n] = new_coord 
		
		bunch =  numpy.zeros([npart],self.data_def)
		
		bunch['Y'] = coords2[:,0]
		bunch['T'] = coords2[:,1]
		bunch['Z'] = coords2[:,2]
		bunch['P'] = coords2[:,3]
		
	#	bunch['Y'] = coords[:,0]
	#	bunch['T'] = coords[:,1]
	#	bunch['Z'] = coords[:,2]
	#	bunch['P'] = coords[:,3]
		self.coords = bunch

	def _twiss_matrix(self, beta_y, beta_z, alpha_y, alpha_z ):
		B = numpy.eye(6)
		B[0,0] = sqrt(beta_y)
		B[1,1] = 1/sqrt(beta_y)
		B[2,2] = sqrt(beta_z)
		B[3,3] = 1/sqrt(beta_z)

		A = numpy.eye(6)
		A[1,0] = -alpha_y
		A[3,2] = -alpha_z

		#M = numpy.dot(A, B)
		M = numpy.dot(B, A)
		return M


	def read_YTZPSD(self, fname, ke):
		dist = numpy.loadtxt(fname)
		nparts = dist.size / 6
		dist = dist.reshape(nparts,6)
		self.coords = numpy.zeros(nparts, self.data_def)
		self.coords['keep'] = True
		self.coords['KE'] = ke
		self.coords['Y'] = dist[:,0]
		self.coords['T'] = dist[:,1]
		self.coords['Z'] = dist[:,2]
		self.coords['P'] = dist[:,3]
		self.coords['X'] = dist[:,4]
		self.coords['D'] = dist[:,5]

	def write_YTZPSD(self, fname):
		nparts = len(self.coords)
		dist = numpy.zeros([nparts, 6])
		dist[:,0] = self.coords['Y']
		dist[:,1] = self.coords['T']
		dist[:,2] = self.coords['Z']
		dist[:,3] = self.coords['P']
		dist[:,4] = self.coords['X']
		dist[:,5] = self.coords['D']
		dist = dist.reshape(nparts*2,3)
		numpy.savetxt(fname, dist)

	
	def get_widths(self):
		y_width = max(self.coords['Y']) - min(self.coords['Y'])
		t_width = max(self.coords['T']) - min(self.coords['T'])
		z_width = max(self.coords['Z']) - min(self.coords['Z'])
		p_width = max(self.coords['P']) - min(self.coords['P'])
		x_width = max(self.coords['X']) - min(self.coords['X'])
		d_width = max(self.coords['D']) - min(self.coords['D'])
		return (y_width, t_width, z_width, p_width, x_width, d_width)

	def __len__(self):
		return len(self.coords)

	def plot(self, fname=None, lims=None):
		pylab.subplot(2,2,1)
		pylab.grid()
		pylab.plot(self.coords['Y'], self.coords['Z'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[2], lims[2])
		plotname = "X-Y (Y-T)"
		pylab.title(plotname)

		
		pylab.subplot(2,2,2)
		pylab.grid()
		pylab.plot(self.coords['Y'], self.coords['T'], ',')
		if lims != None:
			pylab.xlim(-lims[0], lims[0])
			pylab.ylim(-lims[1], lims[1])
		plotname = "X-XP (Y-T"
		pylab.title(plotname)

		pylab.subplot(2,2,3)
		pylab.grid()
		pylab.plot(self.coords['Z'], self.coords['P'], ',')
		if lims != None:
			pylab.xlim(-lims[2], lims[2])
			pylab.ylim(-lims[3], lims[3])
		plotname = "Y-YP (z-P)"
		pylab.title(plotname)

		pylab.subplot(2,2,4)
		pylab.grid()
		pylab.plot(self.coords['X'], self.coords['D'], ',')
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
		"return emittance h and v in m rad"
		r_yt = numpy.sqrt(self.coords['Y']**2 + self.coords['T']**2)
		theta_yt = numpy.arctan2(self.coords['Y'], self.coords['T'])
		major_angle = theta_yt[r_yt.argmax()]
		theta_yt = theta_yt - major_angle
		rot_y = r_yt * numpy.sin(theta_yt)
		rot_t = r_yt * numpy.cos(theta_yt)
		emmitance_h = rot_y.max()*rot_t.max()

		r_zp = numpy.sqrt(self.coords['Z']**2 + self.coords['P']**2)
		theta_zp = numpy.arctan2(self.coords['Z'], self.coords['P'])
		major_angle = theta_zp[r_zp.argmax()]
		theta_zp = theta_zp - major_angle
		rot_z = r_zp * numpy.sin(theta_zp)
		rot_p = r_zp * numpy.cos(theta_zp)
		emmitance_v = rot_z.max()*rot_p.max()
		#print "Emmitance (h, v):", emmitance_h, emmitance_v
		return (emmitance_h, emmitance_v)


	def get_twiss(self, emittance):
		# emittance may be a single number, or tuple (emittance_h, emittance_v)
		try:
			emittance_h, emittance_v = emittance
		except TypeError:
			emittance_h = emittance_v = emittance
		widths = self.get_widths()
		beta_h = (widths[0]/2)**2 / emittance_h
		beta_v = (widths[2]/2)**2 / emittance_v
		#print "beta", beta_h, beta_v
		gamma_h = (widths[1]/2)**2 / emittance_h
		gamma_v = (widths[3]/2)**2 / emittance_v
		# get T of particle with bigest Y
		y_p = self.coords['T'][self.coords['Y'].argmax()]
		alpha_h = - y_p / sqrt(emittance_h/beta_h)
		z_p = self.coords['P'][self.coords['Z'].argmax()]
		alpha_v = - z_p / sqrt(emittance_v/beta_v)
		#print "gamma", gamma_h, gamma_v	
		# following is dangerous, due do stat fluctuations causing roots of negative numbers
		#alpha_h = sqrt((gamma_h*beta_h)-1 )
		#alpha_v = sqrt((gamma_v*beta_v)-1 )
		#print "twiss (bh,ah,bv,av):", beta_h, alpha_h, beta_v, alpha_v
		return beta_h, alpha_h, beta_v, alpha_v


