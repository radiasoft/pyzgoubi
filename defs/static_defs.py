

		
# beam line elements that are not generated from defs
class MARKER(zgoubi_element):
	def __init__(self, name="" ,**settings):
		self._name = name
		self._params = {}
		self._params['plt'] = True
		self._params['name'] = name
		self._zgoubi_name = "MARKER"
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		assert(len(str(self._name)) > 0)
		plt = ""
		if self.plt:
			plt = ".plt"
		return "'MARKER' " + self.name + " " + plt + nl

class OBJET1(zgoubi_element):
	"""Beam made of grid of particles, default params give just 1 reference particle.
	Equivilent to OBJET with a KOBJ=1
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._params = {}
		self._params['BORO'] = 0
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('I'+p)] = 1
			self._params[str('P'+p)] = 0
			self._params[str(p+'R')] = 0 # ref should be 0
		self._params[str('DR')] = 1 # appart from dp/p
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		assert(self.IY * self.IT * self.IZ * self.IP * self.IX * self.ID < 10000)
		assert(max(self.IY, self.IT, self.IZ, self.IP, self.IX,self.ID) < 41 )
		
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += "1" + nl
		out += i(self.IY) +' '+ i(self.IT) +' '+ i(self.IZ) +' '+ i(self.IP) +' '+ i(self.IX) +' '+ i(self.ID) +nl
		out += f(self.PY) +' '+ f(self.PT) +' '+ f(self.PZ) +' '+ f(self.PP) +' '+ f(self.PX) +' '+ f(self.PD) +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR) +nl
		return out


class OBJET2(zgoubi_element):
	"""Beam made of particles, with coords give explicitly
	Equivilent to OBJET with a KOBJ=2
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._params = {}
		self._params['BORO'] = 0
		self._params['particles'] = []
		self._params['IMAX'] = 0
		self._params['IDMAX'] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)
		self.sorted = False
		
	def add(self, **settings):
		"add a particle"
		particle_coords = dict(Y=0, T=0, Z=0, P=0, X=0, D=1, LET=' ')
		for k, v in settings.items():
			particle_coords[k] = v
			
		self._params['particles'].append(particle_coords)
		# keep list sorted by D, needed for output()
		if self.sorted:
			self._params['particles'].sort(key=itemgetter('D'))
		
		
	def clear(self):
		"remove all particles"
		self._params['particles'] = []
	
	def output(self):
		out=''
		self.IMAX = len(self.particles)
		assert(self.IMAX <= 10000)
		assert(self.IMAX >= 1)
		
		#count unique 'D' values
		if self.sorted:
			self.IDMAX = len(set([x['D'] for x in self.particles]))
		else:
			self.IDMAX = self.IMAX
		
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out += "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += "2" + nl
		out += i(self.IMAX) +' '+ i(self.IDMAX) + nl
		
		for part in self.particles:
			out += f(part['Y']) +' '+ f(part['T']) +' '+ f(part['Z']) +' '
			out += f(part['P']) +' '+ f(part['X']) +' '+ f(part['D'])
			#if (part['LET'] != ''):
			out += " '"+ part['LET'] + "'"
			out += nl
		
		# assume that we want to track all particles
		for x in xrange(len(self.particles)):
			out += i(1) + ' '
			if ((x+1)%10 == 0):
				out += nl # add a new line after 10 values
		
		#print out
		return out

		
class OBJET3(zgoubi_element):
	"""Read beam coordinates from file. KOBJ=3 read from formatted file.
		KOBJ=3.01 read from simple text file
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._params = {}
		self._params['BORO'] = 0
		self._params['IT1'] = 1
		self._params['IT2'] = 1
		self._params['ITStep'] = 1
		self._params['IP1'] = 1
		self._params['IP2'] = 1
		self._params['IPStep'] = 1
		self._params['TAG'] = '*'
		self._params['FNAME'] = 'nofilename'
		self._params['FTYPE'] = 'unspecified'
		for p in ['Y','T','Z','P','X','D','Ti']:
			self._params[str(p+'F')] = 1.0
			self._params[str(p+'R')] = 0.0
		self._params['InitC'] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		if self.FTYPE == 'unformatted':
			out += "3.01" + nl
		elif self.FTYPE == 'formatted':
			out += "3" + nl
		else:
			print "Error - specify FTYPE, formatted or unformatted"
			sys.exit()

		out += i(self.IT1) +' '+ i(self.IT2) +' '+ i(self.ITStep) +nl
		out += i(self.IP1) +' '+ i(self.IP2) +' '+ i(self.IPStep) +nl
		out += f(self.YF)+' '+f(self.TF)+' '+f(self.ZF)+' '+ f(self.PF) +' '+ f(self.XF)+' '+f(self.DF)+' '+ f(self.TiF)+' '+ self.TAG +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR)+' '+f(self.TiR)+nl
		out += i(self.InitC) +nl
		out += self.FNAME +nl
		return out


class OBJET5(zgoubi_element):
	"""Beam made of particles, for use with matrix
	Equivilent to OBJET with a KOBJ=5
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "OBJET"
		self._params = {}
		self._params['BORO'] = 0
		self._params['particles'] = []
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('P'+p)] = 1
			self._params[str(p+'R')] = 0
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'OBJET'" +nl
		out += f(self.BORO) +nl
		out += "5" + nl
		out += f(self.PY) +' '+ f(self.PT) +' '+ f(self.PZ) +' '+ f(self.PP) +' '+ f(self.PX) +' '+ f(self.PD) +nl
		out += f(self.YR) +' '+ f(self.TR) +' '+ f(self.ZR) +' '+ f(self.PR) +' '+ f(self.XR) +' '+ f(self.DR) +nl
		return out


class MCOBJET3(zgoubi_element):
	"""Monte Carlo generation of a 6-D object on a phase space ellipse
	Equivilent to MCOBJET with a KOBJ=3
	"""
	def __init__(self ,**settings):
		self._zgoubi_name = "MCOBJET"
		self._params = {}
		self._params['BORO'] = 0
		self._params['IMAX'] = 0
		self._params['particles'] = []
		for p in ['Y','T','Z','P','X','D']:
			self._params[str('K'+p)] = 1
			self._params[str(p+'0')] = 0
		for p in ['y','z','x']:
			self._params[str('alpha_'+p)] = 0
			self._params[str('beta_'+p)] = 0
			self._params[str('emit_'+p)] = 0
			self._params[str('n_cutoff_'+p)] = 0
			self._params[str('n_cutoff2_'+p)] = 0
		self._params[str('I1')] = 139339
		self._params[str('I2')] = 139397
		self._params[str('I3')] = 178393
		object.__setattr__(self, "ready", True)
		self.set(settings)

	def output(self):
		# local short cuts for long function names
		f = self.f2s
		i = self.i2s
		
		out = "'MCOBJET'" +nl
		out += f(self.BORO) +nl
		out += "3" + nl
		out += i(self.IMAX) + nl
		out += i(self.KY) +' '+ i(self.KT) +' '+ i(self.KZ) +' '+ i(self.KP) +' '+ i(self.KX) +' '+ i(self.KD) +nl
		out += f(self.Y0) +' '+ f(self.T0) +' '+ f(self.Z0) +' '+ f(self.P0) +' '+ f(self.X0) +' '+ f(self.D0) +nl

		if self.n_cutoff_y < 0:
			out += f(self.alpha_y) +' '+ f(self.beta_y) +' '+ f(self.emit_y) +' '+ i(self.n_cutoff_y) +' '+i(self.n_cutoff2_y) +nl
		else:
			out += f(self.alpha_y) +' '+ f(self.beta_y) +' '+ f(self.emit_y) +' '+ i(self.n_cutoff_y) +nl

		if self.n_cutoff_z < 0:
			out += f(self.alpha_z) +' '+ f(self.beta_z) +' '+ f(self.emit_z) +' '+ i(self.n_cutoff_z) +' '+i(self.n_cutoff2_z) +nl
		else:
			out += f(self.alpha_z) +' '+ f(self.beta_z) +' '+ f(self.emit_z) +' '+ i(self.n_cutoff_z) +nl

		if self.n_cutoff_x < 0:
			out += f(self.alpha_x) +' '+ f(self.beta_x) +' '+ f(self.emit_x) +' '+ i(self.n_cutoff_x) +' '+i(self.n_cutoff2_x) +nl
		else:
			out += f(self.alpha_x) +' '+ f(self.beta_x) +' '+ f(self.emit_x) +' '+ i(self.n_cutoff_x) +nl

		out += i(self.I1) +' '+ i(self.I2) +' '+ i(self.I3) +nl
		return out
		
# generate classes defined in the defs file
#exec(generate_classes())



# define some useful particles
# constants defined in zgoubi_constants.py
class ELECTRON(zgoubi_element):
	def __init__(self):
		self._params = {}
		self._zgoubi_name = "PARTICUL"
		
	def output(self):
		f = self.f2s
		i = self.i2s
		
		out = "'PARTICUL'" +nl
		out += f(ELECTRON_MASS/1e6) +' '+ f(ELECTRON_CHARGE) +' '+ f(ELECTRON_ANOM_MAG_MOM) +' '+ f(ELECTRON_HALF_LIFE) +' 0' +nl
		return out

class PROTON(zgoubi_element):
	def __init__(self):
		self._params = {}
		self._zgoubi_name = "PARTICUL"
		
	def output(self):
		f = self.f2s
		i = self.i2s
		
		out = "'PARTICUL'" +nl
		out += f(PROTON_MASS/1e6) +' '+ f(PROTON_CHARGE) +' '+ f(PROTON_ANOM_MAG_MOM) +' '+ f(PROTON_HALF_LIFE) +' 0' +nl
		return out

class MUON(zgoubi_element):
	def __init__(self):
		self._params = {}
		self._zgoubi_name = "PARTICUL"
               
	def output(self):
		f = self.f2s
		i = self.i2s

		out = "'PARTICUL'" +nl
		out += f(MUON_MASS/1e6) +' '+ f(MUON_CHARGE) +' '+ f(MUON_ANOM_MAG_MOM) +' '+ f(MUON_HALF_LIFE) +' 0' +nl
		return out

class IMMORTAL_MUON(zgoubi_element):
	"non decaying muon"
	def __init__(self):
		self._params = {}
		self._zgoubi_name = "PARTICUL"
               
	def output(self):
		f = self.f2s
		i = self.i2s

		out = "'PARTICUL'" +nl
		out += f(MUON_MASS/1e6) +' '+ f(MUON_CHARGE) +' '+ f(MUON_ANOM_MAG_MOM) +' '+ f(0) +' 0' +nl
		return out
		
class FAKE_ELEM(zgoubi_element):
	def __init__(self, data=""):
		self.data = data
		self._params = {}
		self._zgoubi_name = "####"

	def output(self):
		return self.data+'\n'
