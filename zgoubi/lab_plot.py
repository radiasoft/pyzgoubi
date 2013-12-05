from __future__ import division
from math import *
import numpy as np
from zgoubi.core import zlog
import scipy.interpolate
import scipy.spatial

null_elements = "END FAISCEAU FAISCNL FAISTORE MCOBJET OBJET PARTICUL ".split()
rect_elements = "DRIFT MULTIPOL QUADRUPO BEND".split()



class LabPlotElement(object):
	def __init__(self, z_element, prev_coord, prev_angle, boro, sector_width):
		self.z_element = z_element
		self.prev_coord = list(prev_coord)
		self.prev_angle = prev_angle

		self.entry_coord = list(prev_coord) # coord of the entry of the ref line (exit of previous element)
		self.entry_angle = prev_angle # angle of the entry of the ref line
		self.element_type =  z_element._zgoubi_name

		zlog.debug("LabPlotElement(%s %s %s)"%(z_element._zgoubi_name, prev_coord, prev_angle))

		self.exit_coord = list(self.entry_coord) # coord of the exit of the ref line
		self.exit_angle = self.entry_angle #  angle of the exit of the ref line

		if self.element_type in null_elements:
			# does not effect locations
			pass

		elif self.element_type in rect_elements:
			# step by length, no angle change
			self.length = self.z_element.XL
			self.exit_coord = self.transform(self.length,0)
			self.entrance_wedge_angle = 0
			self.exit_wedge_angle = 0

			if self.element_type in "MULTIPOL QUADRUPO".split():
				self.width =  self.z_element.R_0
			else:
				self.width = 20

			if self.element_type in "MULTIPOL QUADRUPO".split():
				if self.z_element.KPOS != 1:
					raise ValueError("Only %s with KPOS=1 is currently implemented"%self.element_type)
			elif self.element_type in "BEND":
				if boro == None : raise ValueError("Must set boro for %s"%self.element_type)
				angle = asin(self.z_element.B1 * self.z_element.XL / 2 / boro)
				self.entrance_wedge_angle = self.z_element.W_E - angle
				self.exit_wedge_angle = self.z_element.W_S - angle 
				if self.z_element.KPOS  == 1:
					pass
				elif self.z_element.KPOS  == 3:
					self.entry_angle -= angle
					self.exit_angle -= 2*angle
					self.exit_coord = self.transform(self.length,0)
				else:
					raise ValueError("Only %s with KPOS=1 or 3 is currently implemented"%self.element_type)




		elif self.element_type == "CHANGREF":
			# change reference element
			# angle
			self.exit_angle += radians(self.z_element.ALE)
			# X
			self.exit_coord[0] += self.z_element.XCE * cos(self.exit_angle)
			self.exit_coord[1] += self.z_element.XCE * sin(self.exit_angle)
			# Y
			self.exit_coord[0] += self.z_element.YCE * -sin(self.exit_angle)
			self.exit_coord[1] += self.z_element.YCE * cos(self.exit_angle)


		elif self.element_type in ["DIPOLE", "DIPOLES"]:
			self.dip_at = radians(self.z_element.AT) # sector angle of magnet region
			self.dip_re = self.z_element.RE # radius at entry
			self.dip_rs = self.z_element.RS # radius at exit
			self.dip_te = self.z_element.TE # radius at entry
			self.dip_ts = self.z_element.TS # radius at exit
			self.width = min(self.dip_re, 500) # FIXME need proper method for setting widths
			if sector_width:
				self.width = float(sector_width)

			if self.dip_te != 0 or self.dip_ts != 0:
				raise ValueError("Non zero TE or TS not implemented for %s"%self.element_type)

			# assume local coords have origin at entry, with Y in middle of magnet pointing and machine center
			self.entry_angle += self.dip_at / 2
			self.exit_angle -= self.dip_at
			self.sector_center_coord = list(self.entry_coord)
			self.sector_center_coord[0] += self.dip_re * sin(self.entry_angle - self.dip_at/2)
			self.sector_center_coord[1] -= self.dip_re * cos(self.entry_angle - self.dip_at/2)

			self.exit_coord[0] = self.sector_center_coord[0] + self.dip_rs * sin(-self.exit_angle)
			self.exit_coord[1] = self.sector_center_coord[1] + self.dip_rs * cos(-self.exit_angle)
			

		elif self.element_type in ["POLARMES"]:
			# to get opening angle read the list of angles in the field map
			self.fmap_file_path = self.z_element.FNAME
			fmap_file = open(self.fmap_file_path)
			for n in xrange(4): line = fmap_file.readline()
			phi_line = fmap_file.readline()
			phi_line = phi_line.split()
			self.dip_at = float(phi_line[-1])

			self.dip_re = self.dip_rs = 0
			self.width =  500
			if sector_width:
				self.width = float(sector_width)

			self.entry_angle -= self.dip_at / 2
			self.exit_angle -= self.dip_at

			self.sector_center_coord = list(self.entry_coord)
			self.exit_coord[0] = self.entry_coord[0]
			self.exit_coord[1] = self.entry_coord[1]
			

		else:
			raise ValueError("Can't handle element "+ self.element_type)
	
	def transform(self, x, y):
		if self.element_type not in ["DIPOLE", "DIPOLES", "POLARMES"]:
			x0, y0 = self.entry_coord # FIXME how to handle transform in changref
			a0 = self.entry_angle

			x1 = x0 + x * cos(a0) - y * sin(a0)
			y1 = y0 + y * cos(a0) + x * sin(a0)
		else:
			# coords are in polar
			x0, y0 = self.sector_center_coord
			a0 = self.prev_angle

			x1 = x0 + y * sin(-a0 + x)
			y1 = y0 + y * cos(-a0 + x)
		return [x1, y1]

	def draw_ref_line(self, lpd):
		t = self.transform
		if self.element_type in rect_elements:
			points = [t(0,0), t(self.length,0)]
			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "k-")
		if self.element_type in ["DIPOLE", "DIPOLES"]:
			re = self.dip_re
			a = self.dip_at
			arcsteps = 20
			points = []
			for a1 in np.linspace(0,a,arcsteps):
				points.append(t(a1, re))
			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "k-")

	
	def draw_outline(self, lpd):
		t = self.transform
		if (self.element_type in rect_elements
			and self.element_type != 'DRIFT'):
			if self.entrance_wedge_angle == 0 and self.exit_wedge_angle == 0:
				points = [t(0,self.width/2),
						  t(0,-self.width/2),
						  t(self.length,-self.width/2),
						  t(self.length,self.width/2),
						  t(0,self.width/2)]
			elif abs(self.entrance_wedge_angle) > pi/2:
				raise ValueError("entrance_wedge_angle of %s greater than pi/2"%self.element_type)
			elif abs(self.exit_wedge_angle) > pi/2:
				raise ValueError("exit_wedge_angle of %s greater than pi/2"%self.element_type)
			else:
				entry_offset = self.width/2 * sin(self.entrance_wedge_angle)
				exit_offset = self.width/2 * sin(self.exit_wedge_angle)
				points = [t(0+entry_offset,self.width/2),
						  t(0-entry_offset,-self.width/2),
						  t(self.length+exit_offset,-self.width/2),
						  t(self.length-exit_offset,self.width/2),
						  t(0+entry_offset,self.width/2)]

			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "b-")

		if self.element_type in ["DIPOLE", "DIPOLES","POLARMES"]:
			# in polar
			re = self.dip_re
			w = self.width
			if self.element_type in ["POLARMES"]:
				wp = self.width
				wm = 0
			else:
				wp = self.width/2
				wm = -self.width/2
			a = self.dip_at
			arcsteps = 20
			points = [ t(0, re + wm),
			           t(0, re + wp)]
			for a1 in np.linspace(0,a,arcsteps):
				points.append(t(a1, re + wp))
			points.append(t(a, re + wp))
			points.append(t(a, re + wm))
			for a1 in np.linspace(a,0,arcsteps):
				points.append(t(a1, re + wm))

			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "b-")

		
class LabPlotDrawer(object):
	def __init__(self, mode="matplotlib"):
		self.mode = mode
		
		if self.mode == "matplotlib":
			global matplotlib,plt, Line2D
			import matplotlib
			import matplotlib.pyplot as plt
			from matplotlib.lines import Line2D

			self.fig = plt.figure()
			self.fig.clf()
			self.ax = self.fig.add_subplot(111)
			self.ax.set_aspect('equal', adjustable='datalim')
		else:
			raise ValueError("Can't handle mode "+ self.mode)

	
	def draw_line(self, xs,ys, style, linewidth=1):
		xs = np.array(xs)
		ys = np.array(ys)
		#if np.any(xs < -10): raise ValueError
		if self.mode == "matplotlib":
			self.ax.plot(xs, ys, style, linewidth=linewidth)
			
		else: ValueError("Can't handle mode "+ self.mode)
	
	def draw_label(self, x, y, l, marker=""):
		if marker != "":
			self.ax.plot([x],[y], marker)
		self.ax.annotate(str(l), (x,y))
	
	def draw_im(self, im, extent, cm, vmin, vmax, colorbar=True, colorbar_label=""):
		i = self.ax.imshow(im, extent=extent, origin='lower', cmap=plt.get_cmap(cm),vmin=vmin, vmax=vmax)
		if colorbar:
			cbar = self.fig.colorbar(i)
			if colorbar_label:
				cbar.ax.set_ylabel(colorbar_label)

	def finish(self):
		if self.mode == "matplotlib":
			plt.xlabel("Lab x (cm)")
			plt.ylabel("Lab y (cm)")
			version = matplotlib.__version__.split(".")
			if int(version[0]) >= 1 and int(version[1]) >= 1:
				plt.tight_layout()

	def show(self):
		if self.mode == "matplotlib":
			plt.show()
		else: ValueError("Can't handle mode "+ self.mode)

	def save(self,fname):
		if self.mode == "matplotlib":
			plt.savefig(fname)
		else: ValueError("Can't handle mode "+ self.mode)


class LabPlot(object):
	"""A plotter for beam lines and tracks.
	
	"""
	def __init__(self, line, boro=None, sector_width=None):
		"""Creates a new plot from the line.
		If using an element that adjusts it shape based on BORO, then it must be passed in
		if sector_width is a number it used for the width of sector elements

		"""
		
		self.line = line
		self.elements = []
		self.element_label1 = []
		self.tracks = []
		self.field_map_data = []
		self.boro = boro
		self.duped_labels = []
		self.sector_width = sector_width
		
		self._scan_line()
		

	def _scan_line(self):
		"""Scan through the line, and make a note of where all the elements are.
		
		"""
		angle = 0
		position = [0, 0]

		for elem in self.line.elements():
			if hasattr(elem, 'label1'):
				label = elem.label1
				if label in self.element_label1:
					self.duped_labels.append(label)
					#zlog.warn("Repeated label '%s'"%label)
			else:
				label = ""
			lpelem = LabPlotElement(elem, position, angle, boro=self.boro, sector_width=self.sector_width)
			self.elements.append(lpelem)
			self.element_label1.append(label)
			angle = lpelem.exit_angle
			position = lpelem.exit_coord


	def draw(self, draw_tracks=True, draw_field_points=False, draw_field_midplane=False, field_component='z', field_steps=100, field_int_mode="kd"):
		self.lpd = LabPlotDrawer()

		if field_component not in ['x','y','z']:
			raise ValueError("field_component should be 'y', 'z' or 'x'")

		if draw_field_midplane:
			if field_int_mode=="griddata":
				label = r"$B_%s$ (kG)"%field_component
				field_map_data = np.array(self.field_map_data)
				points = field_map_data[:,0:3:2]
				if field_component == 'y':
					values = field_map_data[:,3].reshape([-1])
				elif field_component == 'z':
					values = field_map_data[:,4].reshape([-1])
				elif field_component == 'x':
					values = field_map_data[:,5].reshape([-1])

				# http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
				# interested in y-x plane
				ymin,ymax,xmin,xmax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()
				nsteps = field_steps*1j # j because of how mgrid works
				grid_y, grid_x = np.mgrid[ymin:ymax:nsteps, xmin:xmax:nsteps]
				
				# check that paths don't deviate from z=0 plane
				if np.abs(field_map_data[:,1]).max() > 1e-6:
					zlog.warn("Some field points are not at z=0. Plot will be of projection onto z=0")

				int_field = scipy.interpolate.griddata(points, values, (grid_y, grid_x))
				#print int_field.nanmax(), int_field.nanmin(),  np.abs(int_field).nanmax()
				vmax = np.nanmax(np.abs(int_field))
				self.lpd.draw_im(int_field.T, (ymin,ymax,xmin,xmax), cm='bwr', vmin=-vmax, vmax=vmax, colorbar_label=label)
			elif field_int_mode=="kd":			
				label = r"$B_%s$ (kG)"%field_component
				field_map_data = np.array(self.field_map_data)
				points = field_map_data[:,0:3:2]
				if field_component == 'y':
					values = field_map_data[:,3].reshape([-1])
				elif field_component == 'z':
					values = field_map_data[:,4].reshape([-1])
				elif field_component == 'x':
					values = field_map_data[:,5].reshape([-1])

				xmin,xmax,ymin,ymax = points[:,0].min(), points[:,0].max(), points[:,1].min(), points[:,1].max()

				print points.shape
				print points
				print values.shape
				print xmin,xmax,ymin,ymax

				nxsteps = field_steps
				xstep_size = (xmax-xmin) / nxsteps
				nysteps = (ymax-ymin) / xstep_size

				kd = scipy.spatial.cKDTree(points)
				
				field_map = np.zeros([nxsteps, nysteps])
				for nx,x in enumerate(np.linspace(xmin,xmax,nxsteps)):
					for ny,y in enumerate(np.linspace(ymin,ymax,nysteps)):
						# get nearby points
						kd_f, kd_i = kd.query([x,y],k=5, distance_upper_bound=xstep_size*1.2)
						# remove index outside range, then mean not points found in distance_upper_bound
						kd_i = kd_i[kd_i < values.size]
						# remove zero fields
						near_values = values[kd_i][ values[kd_i] != 0 ]
						if near_values.size:
							field_map[nx,ny] = near_values.mean()
				vmax = np.nanmax(np.abs(field_map))
				self.lpd.draw_im(field_map.T, (xmin,xmax,ymin,ymax), cm='bwr', vmin=-vmax, vmax=vmax, colorbar_label=label)
			else:
				raise ValueError("field_int_mode must be griddata or kd")



		for elem in self.elements:
			elem.draw_ref_line(self.lpd)
			elem.draw_outline(self.lpd)

		if draw_tracks:
			for track in self.tracks:
				xs, ys, dummy, dummy, dummy = zip(*track)
				self.lpd.draw_line(xs, ys, "r-", linewidth=0.1)

		if draw_field_points:
			for track in self.tracks:
				for p in track:
					xs, ys, by, bz, bx = p
					if by is None: continue
					if field_component == 'y':
						self.lpd.draw_label(xs, ys, "%.3g"%by, 'rx')
					elif field_component == 'z':
						self.lpd.draw_label(xs, ys, "%.3g"%bz, 'rx')
					elif field_component == 'x':
						self.lpd.draw_label(xs, ys, "%.3g"%bx, 'rx')

		self.lpd.finish()
	
	def show(self):
		self.lpd.show()

	def save(self,fname):
		self.lpd.save(fname)


	def add_tracks(self, ftrack=None, ptrack=None):
		#tracks = []
 		# find the list of particles and laps/passes
		pids = set()
		passes = set()
		noels = set()
		if ftrack is not None:
			pids  |= set(np.unique(ftrack['ID']))
			passes  |= set(np.unique(ftrack['PASS']))
			noels  |= set(np.unique(ftrack['NOEL']))
		if ptrack is not None:
			pids |= set(np.unique(ptrack['ID']))
			passes  |= set(np.unique(ptrack['PASS']))
			noels  |= set(np.unique(ptrack['NOEL']))
		if ftrack == ptrack == None:
			raise ValueError("Must pass a fai track or plt track (or both)")

		pids = sorted(list(pids))
		passes = sorted(list(passes))
		noels = sorted(list(noels))

		# dummy tracks if not passed
		dummy = np.zeros([0], dtype=[('ID',int),('PASS',int),('NOEL',int),('element_label1','S10')])
		if ftrack is None: ftrack = dummy
		if ptrack is None: ptrack = dummy


		# per particle, per lap, per element
		for pid in pids:
			ftrack_p = ftrack[ftrack['ID'] == pid]
			ptrack_p = ptrack[ptrack['ID'] == pid]
			for ppass in passes:
				ftrack_pp = ftrack_p[ftrack_p['PASS'] == ppass]
				ptrack_pp = ptrack_p[ptrack_p['PASS'] == ppass]
				this_track = []
				for noel in noels:
					ftrack_ppn = ftrack_pp[ftrack_pp['NOEL'] == noel]
					ptrack_ppn = ptrack_pp[ptrack_pp['NOEL'] == noel]
					if len(ftrack_ppn) == 0 and len(ptrack_ppn) == 0:
						continue

					label = ""
					if len(ftrack_ppn) > 0: label=ftrack_ppn[0]['element_label1']
					elif len(ptrack_ppn) > 0: label=ptrack_ppn[0]['element_label1']
					else:
						ValueError("No label at element number %s" % noel)
					label = label.strip()

					try:
						el_ind = self.element_label1.index(label)
					except ValueError:
						raise ValueError("Track contains label '%s' not found in line. NOEL=%s"%(label, noel))

					if label in self.duped_labels:
						zlog.warn("Track point at element with duplicated label '%s'. Points may be drawn in wrong element"%label)


					for t in ptrack_ppn:
						if t['IEX'] != 1: break
						y = t['Y']
						x = t['X']
						z = t['Z']
						if isnan(x) or isnan(y) or isnan(z):
							continue
						by, bz, bx = t['BY'], t['BZ'], t['BX']
						if abs(by) < 1e-50 and abs(bz) < 1e-50 and abs(bx) < 1e-50:
							by, bz, bx = 0, 0, 0 # ignore tiny fields
						xt, yt = self.elements[el_ind].transform(x,y)
						this_track.append([xt,yt, by, bz, bx])
						self.field_map_data.append([xt,z,yt,by, bz, bx])
					for t in ftrack_ppn:
						if t['IEX'] != 1: break
						# fai has no x coord, and takes label from element before it
						y = t['Y']
						x = 0

						# index error here may mean plt track passed as fai track, maybe there should be some way to check
						xt, yt = self.elements[el_ind+1].transform(x,y)
						this_track.append([xt,yt,None,None,None])

				#print this_track
				if len(this_track) > 0 :
					self.tracks.append(this_track)

					



