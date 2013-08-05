from __future__ import division
from math import *
import numpy as np
from zgoubi.core import zlog

null_elements = "END FAISCEAU FAISCNL FAISTORE MCOBJET OBJET PARTICUL ".split()
rect_elements = "DRIFT MULTIPOL QUADRUPO".split()



class LabPlotElement(object):
	def __init__(self, z_element, entry_coord, entry_angle):
		self.z_element = z_element
		self.entry_coord = entry_coord
		self.entry_angle = entry_angle
		self.element_type =  z_element._zgoubi_name

		print "LabPlotElement(",z_element._zgoubi_name, entry_coord, entry_angle, ")"

		self.exit_coord = list(self.entry_coord)
		self.exit_angle = self.entry_angle

		if self.element_type in null_elements:
			# does not effect locations
			pass

		elif self.element_type in rect_elements:
			# step by length, no angle change
			self.length = self.z_element.XL
			self.exit_coord = self.transform(self.length,0)
			if self.element_type != 'DRIFT':
				self.width =  self.z_element.B_0 * 2

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


		else:
			raise ValueError("Can't handle element "+ self.element_type)
	
	def transform(self, x, y):
		x0, y0 = self.entry_coord # FIXME how to handle transform in changref
		a0 = self.entry_angle

		x1 = x0 + x * cos(a0) - y * sin(a0)
		y1 = y0 + y * cos(a0) + x * sin(a0)
		return [x1, y1]

	def draw_ref_line(self, lpd):
		t = self.transform
		if self.element_type in rect_elements:
			points = [t(0,0), t(self.length,0)]
			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "k-")

	
	def draw_outline(self, lpd):
		t = self.transform
		if (self.element_type in rect_elements
			and self.element_type != 'DRIFT'):

			points = [t(0,self.width/2),
			          t(0,-self.width/2),
			          t(self.length,-self.width/2),
			          t(self.length,self.width/2),
			          t(0,self.width/2)]

			xs, ys = zip(*points)
			lpd.draw_line(xs, ys, "b-")



		
class LabPlotDrawer(object):
	def __init__(self, mode="matplotlib"):
		self.mode = mode
		
		if self.mode == "matplotlib":
			global plt, Line2D
			import matplotlib.pyplot as plt
			from matplotlib.lines import Line2D

			self.fig = plt.figure()
			self.ax = self.fig.add_subplot(111)
			self.ax.set_aspect('equal', adjustable='datalim')
		else:
			raise ValueError("Can't handle mode "+ self.mode)

	
	def draw_line(self, xs,ys, style):
		if self.mode == "matplotlib":
			self.ax.plot(xs, ys, style)
			
		else: ValueError("Can't handle mode "+ self.mode)


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
	def __init__(self, line):
		"""Creates a new plot from the line.

		"""
		
		self.line = line
		self.elements = []
		self.element_label1 = []
		self.tracks = []

		
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
					zlog.warn("Repeated label '%s'"%label)
			else:
				label = ""
			lpelem = LabPlotElement(elem, position, angle)
			self.elements.append(lpelem)
			self.element_label1.append(label)
			angle = lpelem.exit_angle
			position = lpelem.exit_coord


	def draw(self):
		self.lpd = LabPlotDrawer()
		for elem in self.elements:
			elem.draw_ref_line(self.lpd)
			elem.draw_outline(self.lpd)

		for track in self.tracks:
			xs, ys = zip(*track)
			self.lpd.draw_line(xs, ys, "r-")


		#self.lpd.show()
	
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
			print ftrack.dtype.names
			pids  |= set(np.unique(ftrack['ID']))
			passes  |= set(np.unique(ftrack['PASS']))
			noels  |= set(np.unique(ftrack['NOEL']))
		if ptrack is not None:
			print ptrack.dtype.names
			pids |= set(np.unique(ptrack['ID']))
			passes  |= set(np.unique(ptrack['PASS']))
			noels  |= set(np.unique(ptrack['NOEL']))
		if ftrack == ptrack == None:
			raise ValueError("Must pass a fai track or plt track (or both)")

		pids = sorted(list(pids))
		passes = sorted(list(passes))
		noels = sorted(list(noels))

		# dummy tracks if not passed
		dummy = np.zeros([0], dtype=[('ID',int),('PASS',int),('NOEL',int)])
		if ftrack is None: ftrack = dummy
		if ptrack is None: ptrack = dummy

		print pids, passes, noels
		print np.unique(ftrack['element_label1'])
		print np.unique(ptrack['element_label1'])

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

					label = ""
					if len(ftrack_ppn) > 0: label=ftrack_ppn[0]['element_label1']
					elif len(ptrack_ppn) > 0: label=ptrack_ppn[0]['element_label1']
					else:
						ValueError("No label at element number %s" % noel)
					label = label.strip()
					print label

					try:
						el_ind = self.element_label1.index(label)
					except ValueError:
						raise ValueError("Track contains label '%s' not found in line"%label)

					for t in ptrack_ppn:
						y = t['Y']
						x = t['X']
						xt, yt = self.elements[el_ind].transform(x,y)
						this_track.append([xt,yt])
					for t in ftrack_ppn:
						# fai has no x coord, and takes label from element before it
						y = t['Y']
						x = 0
						xt, yt = self.elements[el_ind+1].transform(x,y)
						this_track.append([xt,yt])

				#print this_track
				self.tracks.append(this_track)

					



