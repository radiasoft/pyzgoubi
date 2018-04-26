#!/usr/bin/env python
# -*- coding: utf-8 -*-

#		pyzgoubi - python interface to zgoubi
#		Copyright 2008-2015 Sam Tygier <Sam.Tygier@hep.manchester.ac.uk>
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

"""Contains the Line and Results objects. Also does all the Pyzgoubi setup needed, eg reading settings, generating definitions, and importing definitions.

"""


from __future__ import division
from string import ascii_letters, digits
from math import *
import tempfile
import shutil
import os
import sys
import time
import re
import struct
from glob import glob
import copy
import threading
import Queue
import subprocess
import weakref
import warnings
dep_warn = " is deprecated, and will likely be removed in a future version. Please see the upgrade notes in the docs for more info."
import logging
logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(filename)s:%(lineno)d - %(funcName)s(): %(message)s")
zlog = logging.getLogger('PyZgoubi')


try:
	import numpy
except ImportError:
	zlog.error("Numpy is required for PyZGoubi")
	sys.exit(1)
try:
	from operator import itemgetter
except ImportError:
	zlog.error("PyZgoubi requires Python 2.5 or newer (2.7 recommended)")
	sys.exit(1)

from zgoubi.constants import *
from zgoubi.common import *
from zgoubi.exceptions import *
import zgoubi.io as io
import zgoubi.bunch


from zgoubi.settings import zgoubi_settings

# in python3 we will be able to use zlog.setLevel(zgoubi_settings['log_level']), can be changed in pyzgoubi too
zlog.setLevel(logging._levelNames[zgoubi_settings['log_level']])

sys.setcheckinterval(10000)

zgoubi_module_path = os.path.dirname(os.path.realpath(__file__))
# something like
# $PREFIX/lib/python2.6/site-packages/zgoubi
# $PREFIX\Lib\site-packages\zgoubi

zgoubi_path = zgoubi_settings['zgoubi_path']

pyzgoubi_egg_path = glob(os.path.dirname(zgoubi_module_path) + "/pyzgoubi*egg-info")

#definitions of elements
home = os.path.expanduser('~')
config_dir = os.path.join(home, ".pyzgoubi")
if not os.path.exists(config_dir):
	os.mkdir(config_dir)

def_cache_dir = os.path.join(config_dir, "def_cache")
if not os.path.exists(def_cache_dir):
	os.mkdir(def_cache_dir)

compiled_defs_path = os.path.join(def_cache_dir, "user_defs.py")

definitions_paths = [os.path.join(config_dir, x) for x in os.listdir(config_dir) if x.endswith('.defs')]
definitions_paths += zgoubi_settings['extra_defs_files']
nl = '\n'

#if needed recompile defs
if not os.path.exists(compiled_defs_path):
	zlog.debug("no compiled defs. compiling")
	need_def_compile = True
else:
	need_def_compile = False
	for f in definitions_paths:
		if not os.path.exists(f):
			zlog.error("Definitions file: "+f+" does not exist")
		if os.path.exists(f) and os.path.getmtime(f) >= os.path.getmtime(compiled_defs_path):
			zlog.debug("need to recompile"+ f)
			need_def_compile = True
	if not need_def_compile: # also need to recompile after a install or update
		for f in pyzgoubi_egg_path:
			if os.path.exists(f) and os.path.getmtime(f) >= os.path.getmtime(compiled_defs_path):
				zlog.debug("pyzgoubi first run, compiling defs")
				need_def_compile = True
				break
if need_def_compile:
	from zgoubi import makedefs
	zlog.debug("Compiling definitions")
	makedefs.make_element_classes(definitions_paths, compiled_defs_path)

max_label_size = zgoubi_settings["max_label_size"]


def yield_n_lines(fh, n):
	"yields n lines at a time"

	lines = []
	for line in fh:
		lines.append(line)
		if len(lines) == n:
			yield lines
			lines = []

	yield lines


def read_n_lines(fh, n):
	"reads n lines at a time"

	lines = []
	for dummy in xrange(n):
		lines.append(fh.readline())
	return lines

def trans_to_regex(fmt):
	"Transform a printf style format in to a regular expression"
	fmt = fmt.replace('%c', '(.)')
	fmt = re.sub(r'%(\d+)c', r'(.{\1})', fmt)
	fmt = fmt.replace('%d', r'([-+]?\d+)')
	fmt = fmt.replace('%e', r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%E', r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%f', r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%g', r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%i', r'([-+]?(?:0[xX][\dA-Fa-f]+|0[0-7]*|\d+))')
	fmt = fmt.replace('%o', '(0[0-7]*)')
	fmt = fmt.replace('%s', r'(\S+)')
	fmt = fmt.replace('%u', r'(\d+)')
	fmt = fmt.replace('%x', r'(0[xX][\dA-Fa-f]+")')
	fmt = fmt.replace('%X', r'(0[xX][\dA-Fa-f]+")')
	#print "fmt=", fmt
	return re.compile(fmt)

def scanf(in_str, fmt):
	"something like scanf, used to read the fai and plt files in Results.get_all()"
	if type(fmt) == type(re.compile('')):
		res = fmt.search(in_str)
	else:
		res = trans_to_regex(fmt).search(in_str)
	if res is None:
		return None
	else:
		return res.groups()

# a base class for all the beam line objects
class zgoubi_element(object):
	"A base class for zgoubi elements"
	def __init__(self):
		pass
	
	def set(self, *dsettings, **settings):
		"""Set a parameter value::
			my_element.set(XL=5)

		can also use a dictionary to set values::
			s = {'XL':5, B_0=0.2}
			my_element.set(s)
		
		"""
		#try to merge the two dicts
		try:
			settings.update(dsettings[0])
		except IndexError:
			pass
			
		for key, val in settings.items():
			#self.__setattr__(key,val)
			if key in self._params.keys():
				self._params[key] = val
			else:
				raise ValueError("no such param: '" + str(key) + "' In element " + self._zgoubi_name)
	def get(self, key):
		"Get a parameter"
		return self._params[key]

	#FIXME investigae if these would be faster as staticmethods
	def f2s(self, f):
		"format float for printing"
		#out = "%e" % float(f)
		#out = "%s" % float(f)
		out = "%.12e" % float(f)
		return out
		
	def i2s(self, i):
		"format integer for printing"
		out = str(int(i))
		return out

	def l2s(self, l):
		"format label for printing"
		out = l[:max_label_size]
		return out

	def x2s(self, i):
		"format xpas for printing"
		try:
			out = self.f2s(i)
		except TypeError:
			out = '#'+self.i2s(i[0])+ "|"+self.i2s(i[1])+ "|"+self.i2s(i[2])
		return out

	def __getattr__(self, name):
		"allow dot access to parameters"
		if name != '_params':
			try:
				return self._params[name]
			except KeyError:
				pass
		return object.__getattribute__(self, name)

	def list_params(self):
		"Return a list of parameters"
		return list(self._params)
	
	def reverse(self):
		"Flip the element along the beam line direction, i.e. the entrance and exit properties are swapped"
		if self._zgoubi_name in  ["DIPOLES", "FFAG"]:
			sub_swap_pairs = "G0_E,G0_S KAPPA_E,KAPPA_S NCE,NCS CE_0,CS_0 CE_1,CS_1 CE_2,CS_2 CE_3,CS_3 CE_4,CS_4 CE_5,CS_5 SHIFT_E,SHIFT_S OMEGA_E,OMEGA_S THETA_E,THETA_S R1_E,R1_S U1_E,U1_S U2_E,U2_S R2_E,R2_S"
			for sub_element in self._looped_data:
				for swap_pair in sub_swap_pairs.split():
					p1, p2 = swap_pair.split(",")
					sub_element[p1], sub_element[p2] = sub_element[p2], sub_element[p1]
				sub_element["ACN"] = self._params["AT"] - sub_element["ACN"]
				sub_element["OMEGA_E"] *= -1
				sub_element["OMEGA_S"] *= -1
				sub_element["THETA_E"] *= -1
				sub_element["THETA_S"] *= -1
		elif self._zgoubi_name == "CHANGREF":
			self._params["YCE"] *= -1

	def __neg__(self):
		new_e = copy.deepcopy(self)
		new_e.reverse()
		return new_e

	def set_plot_hint(self, **hints):
		"Add hints to help lab_plot"
		if not hasattr(self, "plot_hints"): self.plot_hints = {}
		self.plot_hints.update(hints)



from zgoubi.static_defs import *
from zgoubi.simple_defs import *
sys.path.append(os.path.join(def_cache_dir))
from user_defs import *


try:
	test_element = BEND()
except NameError:
	zlog.error("Elements did not load correctly")
	path_info = " zgoubi_module_path: %s\nzgoubi_path: %s\npyzgoubi_egg_path: %s" % (zgoubi_module_path, zgoubi_path, pyzgoubi_egg_path)
	zlog.error(path_info)
	zlog.error("Try deleting: "+static_defs+ " " +compiled_defs_path)
	exit(1)


class Line(object):
	"The Line object holds a series of elements, to represent an accelerator lattice."

	def __init__(self, name):
		"Create a Line. The name parameter is passed to zgoubi."
		self.element_list = []
		self.name = name
		self.tmp_folders = [] # keep track of old tmp folders
		self.results = [] # results generated by this line
		self.last_result = None
		self.no_more_xterm = False
		self.input_files = []
		
		self.shutil = shutil # need to keep a reference to shutil
							# otherwise the intepreter may have thrown it away
							# by the time the __del__() is run
		self.has_run = False
		self.full_line = False # has an OBJET, dont allow full lines to be added to each other
								# only a full line outputs its name into zgoubi.dat

	def __copy__(self):
		"A shallow copy, contains the same elements"
		new_line = type(self)(self.name)
		new_line.element_list = copy.copy(self.element_list)
		new_line.no_more_xterm = self.no_more_xterm
		new_line.input_files = self.input_files
		new_line.full_line = self.full_line
		return new_line
	
	def __deepcopy__(self, memo):
		"A deep copy, contains the copies of elements"
		new_line = type(self)(self.name)
		new_line.element_list = copy.deepcopy(self.element_list, memo)
		new_line.no_more_xterm = self.no_more_xterm
		new_line.input_files = self.input_files
		new_line.full_line = self.full_line
		return new_line

	def __neg__(self):
		"return a reversed line"
		new_line = copy.copy(self)
		new_line.element_list = copy.deepcopy(new_line.element_list)
		new_line.element_list.reverse()
		for e in new_line.element_list:
			e.reverse()
		new_line.name = "-"+self.name
		return new_line

	def __add__(self, rhs):
		new_line = Line(self.name)
		for element in self.element_list:
			new_line.add(element)
		for element in rhs.element_list:
			new_line.add(element)
		new_line.add_input_files(self.input_files)
		new_line.add_input_files(rhs.input_files)
		return new_line

	def __rmul__(self, lhs):
		new_line = Line(self.name)
		for x in xrange(lhs):
			for element in self.element_list:
				new_line.add(element)
		new_line.add_input_files(self.input_files)
		return new_line

	def __mul__(self, rhs):
		new_line = Line(self.name)
		for x in xrange(rhs):
			for element in self.element_list:
				new_line.add(element)
		new_line.add_input_files(self.input_files)
		return new_line


	def __str__(self, prefix=""):
		"For printing a line"
		out = self.name + "\n"
		for element in self.element_list:
			if isinstance(element, Line):
				out += element.__str__(prefix+" ")
			else:
				out += "%s %s %s %s\n" % (prefix, element._zgoubi_name, element.label1, element.label2)
		return out

	def elements(self):
		"Iterator for elements in line, including sub lines."
		for element in self.element_list:
			try:
				# decend into sub lines
				for sub_element in element.elements():
					yield sub_element
			except AttributeError:
				yield element

	def add(self, *elements):
		"Add an elements to the line. Can also be used to add one line into another."
		for element in elements:
			self.element_list.append(element)
			try:
				if 'OBJET' in element._zgoubi_name:
					self.full_line = True
			except AttributeError:
				pass
			if hasattr(element, "input_files"):
				self.add_input_files(element.input_files)
	
	def check_line(self):
		"Check that line has OBJET or MCOBJET at the start, and an END at the end. Gives warnings otherwise. Called by run() if in debug mode."
		has_end = False
		line_good = True
		for n, element in enumerate(self.elements()):
			if has_end:
				line_good = False
				try:
					zlog.warn("Element (%s) after END", element._zgoubi_name)
				except AttributeError:
					zlog.warn("Element after END")

			isobjet = False
			try:
				if 'OBJET' in element._zgoubi_name:
					isobjet = True
				if element._zgoubi_name == "END":
					has_end = True
			except AttributeError:
				pass

			if n == 0 and not isobjet:
				zlog.warn("First element in line no OBJET/MCOBJET")
				line_good = False
			if n != 0 and isobjet:
				zlog.warn("OBJET/MCOBJET appears as element number %d. (Should only be first)", n)
				line_good = False

			if hasattr(element, 'XPAS'):
				if element.XPAS == 0:
					zlog.warn("Element %s type %s has XPAS=0 (integration step)", n, element._zgoubi_name)
		if not has_end:
			zlog.warn("No END element found")
			line_good = False

		return line_good
				
	
	def full_tracking(self, enable=True, drift_to_multi=False):
		"""Enable full tracking on magnetic elements.
		This works by setting IL=2 for any element with an IL parameter.
		use line.full_tracking(False) to disable tracking
		drift_to_multi=True replaces drifts with weak (B_1 = 1e-99) multipoles

		"""
		if drift_to_multi:
			for element in self.elements():
				if element._zgoubi_name == "DRIFT" and element.XL != 0:
					fake_drift = MULTIPOL(element.label1, element.label2,
					                      XL=element.XL, XPAS=(10, 10, 10),
					                      B_1=1e-99, IL=2, KPOS=1)
					self.replace(element, fake_drift)

		found_trackable_element = False

		for element in self.elements():
			#t = str(type(element)).split("'")[1] #get type, eg zgoubi22.QUADRUPO
			#print t

			# if this gives "RuntimeError: maximum recursion depth exceeded"
			#errors, then the object may not have a parameters dict.
			try:
				#print element.output()
				if enable == True:
					element.set(IL=2)
				elif enable == False:
					element.set(IL=0)
				found_trackable_element = True
			#	print element.output()
			except ValueError:
				pass

		if not found_trackable_element:
			zlog.warn("Did not find any magnets in line")

			
	def remove_looping(self):
		"removes any REBELOTE elements from the line"
		self.element_list = [element for element in self.element_list if ("REBELOTE" not in str(type(element)).split("'")[1])]
				
		
	def output(self):
		"Generate the zgoubi.dat file, and return it as a string"
		out = ""
		if self.full_line:
			out = self.name + nl
		
		for element in self.element_list:
			out += element.output() + nl
		
		return out

		
	def run(self, xterm=False, tmp_prefix=zgoubi_settings['tmp_dir'], silence=False, timer=False):
		"""Run zgoubi on line.
		If xterm is true, stop after running zgoubi, and open an xterm for the user in the tmp dir. From here zpop can be run.
		Returns a :py:class:`Results` object
		"""
		if timer: t0 = time.time()
		if zlog.isEnabledFor(logging.DEBUG):
			self.check_line()
		orig_cwd = os.getcwd()

		tmpdir = tempfile.mkdtemp(prefix="zgoubi_", dir=tmp_prefix)
		zlog.debug("running zgoubi in"+tmpdir)
	
		for input_file in self.input_files:
			src = os.path.join(orig_cwd, input_file)
			dst = os.path.join(tmpdir, os.path.basename(input_file))
			#print src, dst
			# if possible make a symlink instead of copying (should work on linux/unix)
			if os.path.exists(dst):
				os.remove(dst) # can't over write an existing symlink
			try:
				os.symlink(src, dst)
			except AttributeError: # should catch windows systems which don't have symlink
				shutil.copyfile(src, dst)
	
		self.tmp_folders.append(tmpdir)
		#os.chdir(tmpdir)
		
		for element in self.elements():
			# some elements may have a setup function, to be run before zgoubi
			if hasattr(element, "setup"):
				element.setup(tmpdir)
		
		infile = open(tmpdir+"/zgoubi.dat", 'w')
		infile.write(self.output())
		infile.close()

		command = zgoubi_settings['zgoubi_path']
		if timer: t1 = time.time()
		if silence:
			command += " > zgoubi.stdout"
			z_proc = subprocess.Popen(command, shell=True, cwd=tmpdir)
		else:
			z_proc = subprocess.Popen(command, shell=False, cwd=tmpdir)
		exe_result = z_proc.wait()
		if timer: t2 = time.time()

		if exe_result != 0:
			zlog.error("zgoubi failed to run\nIt returned:%s", exe_result)

		if xterm and not self.no_more_xterm:
			print "Do you want an xterm? (y=yes/n=no/s=stop asking)"
			ans = raw_input()
			if ans.startswith('y'):
				subprocess.Popen("xterm", shell=True, cwd=tmpdir).wait()
			elif ans.startswith('s'):
				self.no_more_xterm = True

		if exe_result != 0:
			if silence:
				print open(os.path.join(tmpdir, "zgoubi.sdterr")).read()
			if exe_result == 32512:
				raise ZgoubiRunError("Check that fortran runtime libraries are installed")
			if exe_result == -9:
				raise ZgoubiRunError("Check that you have sufficient RAM or enable memory overcommit")
			if exe_result == 2:
				raise ZgoubiRunError("Fortran runtime error")

		res_file = tmpdir+"/zgoubi.res"
		#output = outfile.read()
		
		for n, line in enumerate(open(res_file)):
			lline = line.lower()
			if "error" in lline or "warning" in lline or "sbr" in lline:
				print "zgoubi.res:", n, ":", line
				if "SBR OBJ3 -> error in  reading  file" in line:
					raise ZgoubiRunError(line)

		#os.chdir(orig_cwd)
		
		element_types = [str(type(element)).split("'")[1].rpartition(".")[2] for element in self.elements()]
		self.has_run = True	
		result = Results(line=self, rundir=tmpdir, element_types=element_types)
		self.results.append(weakref.ref(result))
		self.last_result = result
		if timer:
			result.timer_setup = t1-t0
			result.timer_run = t2-t1

		return result
	
	def track_bunch(self, bunch, binary=False, **kwargs):
		"Track a bunch through a Line, and return the bunch. This function will uses the OBJET_bunch object, and so need needs a Line that does not already have a OBJET. If binary is true then particles are sent to zgoubi in binary (needs a version of zgoubi that supports this)"
		if self.full_line:
			raise BadLineError("If line already has an OBJET use run()")

		bunch_len = len(bunch)
		if bunch_len == 0:
			zlog.error("Bunch has zero particles")
			raise ValueError
		#build a line with the bunch OBJET and segment we were passed
		new_line = Line("bunch_line")
		new_line.add(OBJET_bunch(bunch, binary=binary))
		new_line.add(self)
		#mark the faiscnl that we are interested in
		new_line.add(MARKER("trackbun"))
		new_line.add(FAISCNL(FNAME='b_zgoubi.fai'))
		new_line.add(END())

		# run the line
		result = new_line.run(**kwargs)
		del new_line
		# return the track bunch
		done_bunch = result.get_bunch('bfai', end_label="trackbun", old_bunch=bunch)
		done_bunch_len = len(done_bunch)
		if bunch_len != done_bunch_len:
			zlog.warn("Started with %s particles, finished with %s", bunch_len, done_bunch_len)
		result.clean()
		return done_bunch
		
	def track_bunch_mt(self, bunch, n_threads=4, max_particles=None, binary=False, **kwargs):
		"This function should be used identically to the track_bunch function, apart from the addition of the n_threads argument. This will split the bunch into several slices and run them simultaneously. Set n_threads to the number of CPU cores that you have. max_particle can be set to limit how many particles are sent at a time."
		in_q = Queue.Queue()
		out_q = Queue.Queue()
		if max_particles is None:
			max_particles = 1e3

		def worker(in_q, out_q, work_line, name, stop_flag):
			"A worker function. Gets run in threads"
			while True:
				try:
					start_index, work_bunch = in_q.get(block=True, timeout=0.1)
				except Queue.Empty:
					#print "get() timed out in thread", name
					if stop_flag.is_set():
						#print "exiting thread", name
						return
					continue
				#print "Thread", name, "working"
				try:
					done_bunch = work_line.track_bunch(work_bunch, binary=binary, **kwargs)
				except:
					zlog.error("Exception in track_bunch() thread")
					out_q.put((sys.exc_info()))
				else:
					out_q.put((start_index, done_bunch.particles()))
				in_q.task_done()
				#print "Thread", name, "task done"

		# pre process line output, so it does not have to be done in each thread
		line_output = self.output()
		new_line = Line(self.name)
		new_line.add(FAKE_ELEM(line_output))

		stop_flag = threading.Event()
		for thread_n in xrange(n_threads):
			t = threading.Thread(target=worker,
			                     kwargs={'in_q':in_q, 'out_q':out_q,
			                     'work_line':new_line, 'name':thread_n, 'stop_flag':stop_flag})
			t.setDaemon(True)
			t.start()
			#print "Created thread", x
		
		start_index = 0
		n_tasks = 0
		#print "Queuing work"
		for item in bunch.split_bunch(max_particles=max_particles, n_slices=n_threads):
			#print "Queue task", n_tasks
			in_q.put((start_index, item))
			start_index += len(item)
			n_tasks += 1
		#print "Work queued"
		#in_q.join()
		#print "Work done"

		final_bunch = zgoubi.bunch.Bunch(nparticles=len(bunch), rigidity=bunch.get_bunch_rigidity(), mass=bunch.mass, charge=bunch.charge)
		survive_particles = numpy.zeros(len(bunch), dtype=numpy.bool) # bit map, set true when filling with particles

		# workers may return out of order, so use start_index to put the coords in the correct place
		for x in xrange(n_tasks):
			#print "collecting task", x
			result = out_q.get()
			try:
				start_index, done_bunch = result
			except ValueError:
				import traceback, time
				stop_flag.set()
				time.sleep(2) # give the other threads a couple of seconds, to make output prettier
				zlog.error("Exception retrieved by main thread")
				print
				traceback.print_exception(*result)
				print
				#reraise error message
				raise result[0](result[1])
			out_q.task_done()
			final_bunch.particles()[start_index:start_index+len(done_bunch)] = done_bunch
			survive_particles[start_index:start_index+len(done_bunch)] = True

		if not numpy.all(survive_particles):
			final_bunch.coords = final_bunch.particles()[survive_particles]
			zlog.warn("Started with %s particles, finished with %s", len(bunch), len(final_bunch))


		#print "all done"
		stop_flag.set()
		return final_bunch

	def clean(self):
		"clean up temp directories"
		for result in self.results:
			obj = result()
			#print "in Line.clean", obj
			if obj is not None:
				obj.clean()

		self.results = []
		#for dir in self.tmp_folders:
		#	print "removing", dir
		#	print shutil
		#	self.shutil.rmtree(dir)
		
		#self.tmp_folders = [] # and blank list
		

	def add_input_files(self, file_paths=None, pattern=None):
		"""Add some extra input files to the directory where zgoubi is run.
		This is useful for field map files.
		file_paths can be a string::

			l.add_input_files('map')
		
		an iterable, for example a list::

			l.add_input_files(['map1', 'map2', 'map3'])
		
		To add many files use a pattern eg::

			l.add_input_files(pattern="maps/*")
		
		Will use symlinks when avaliable (Linux/UNIX), falls back to copying otherwise.
		"""

		if file_paths is None:
			file_paths = []

		if hasattr(file_paths, "lower"):
			file_paths = [file_paths]

		if pattern is not None:
			globbed_files = glob(pattern)
			file_paths += globbed_files
		self.input_files += file_paths
	

	def _find_by_index(self, index):
		"""Iterate through sub lines to find element as if indexed in a flat list
		returns [line, index_in_line]"""
		stack = [[self, -1]]
		for n in xrange(index+1):
			if stack[-1][1]+1 >= len(stack[-1][0].element_list):
				# step back up
				stack.pop()
				if len(stack) == 0:
					raise ValueError("Index %s out of range"%index)

			# step along
			stack[-1][1] += 1

			if isinstance(stack[-1][0].element_list[stack[-1][1]], Line):
				# step in
				stack.append([stack[-1][0].element_list[stack[-1][1]], 0])
		return  stack[-1]
		
	def replace(self, elementold, elementnew, select_index=0):
		"""Replace an element in the line. setting select_index to n will replace the nth occurence of that item. If select index is not set, the first occurence is replaced.
		
		Note: elementold must be the same instance as an element in the Line"""

		indices = self.find_elements(elementold)
		if len(indices) == 0:
			raise ValueError("elementold not found in line")
		if select_index >= len(indices):
			raise ValueError("only %s instances of elementold in line"%len(indices))
		index = indices[select_index]
		subline, pos = self._find_by_index(index)
		subline.element_list[pos] = elementnew

	def insert(self, index, *elements):
		"Insert elements into the line before position given by index, treats sub lines as linear list"
		subline, pos = self._find_by_index(index)
		for element in elements:
			subline.element_list.insert(pos, element)

	def prepend(self, *elements):
		"Add a elements to the start of the line"
		for element in elements:
			self.element_list.insert(0, element)

			try:
				if 'OBJET' in element._zgoubi_name:
					self.full_line = True
			except AttributeError:
				pass
			
	def remove(self, index):
		"Remove element at index, treats sub lines as linear list"
		subline, pos = self._find_by_index(index)
		subline.element_list.pop(pos)


	def find_elements(self, element):
		"Returns all the positions of element in line, treats sub lines as linear list"
		indices = []
		for n, e in enumerate(self.elements()):
			if e is element:
				indices.append(n)

		return indices
	
	def get_objet(self):
		"Find the OBJET element"
		for e in self.elements():
			if e._zgoubi_name == "OBJET":
				return e
		raise ValueError("Line has no OBJET element")
			
	
class Results(object):
	"""This class lets you analyse the results after running a line.

	It is created automatically and returned by :py:meth:`Line.run()`

	"""
	def __init__(self, line=None, rundir=None, element_types=None):
		#self.line = line
		self.rundir = rundir
		self.element_types = element_types
		self.shutil = shutil

	def clean(self):
		"clean up temp directory"
		try:
			self.shutil.rmtree(self.rundir)
		except OSError:
			# dir must already be deleted
			pass
	
	def __del__(self):
		self.clean()

	def _get_fh(self, f):
		"return f file handle"
		path = os.path.join(self.rundir, f)
		if not os.path.exists(path):
			raise IOError("No file: %s in %s" % (f, self.rundir))
		fh = open(path)
		return fh
		
	def _get_str(self, f):
		"return f as a string"
		return self._get_fh(f).read()
	
	def _save_file(self, f, path):
		"save f to path"
		spath = os.path.join(self.rundir, f)
		if not os.path.exists(spath):
			raise IOError("No file: %s in %s" % (f, self.rundir))
		shutil.copyfile(spath, path)
	
	#generate specific functions
	def res_fh(self):
		"return file handle for res file"
		return self._get_fh("zgoubi.res")
	def res(self):
		"return res file as string"
		return self._get_str("zgoubi.res")
	def save_res(self, path):
		"save res file to path"
		return self._save_file("zgoubi.res", path)
		
	def plt_fh(self):
		"return file handle for plt file"
		return self._get_fh("zgoubi.plt")
	def plt(self):
		"return plt file as string"
		return self._get_str("zgoubi.plt")
	def save_plt(self, path):
		"save plt file to path"
		return self._save_file("zgoubi.plt", path)
		
	def dat_fh(self):
		"return file handle for dat file"
		return self._get_fh("zgoubi.dat")
	def dat(self):
		"return dat file as string"
		return self._get_str("zgoubi.dat")
	def save_dat(self, path):
		"save dat file to path"
		return self._save_file("zgoubi.dat", path)
		
	def fai_fh(self):
		"return file handle for fai file"
		return self._get_fh("zgoubi.fai")
	def fai(self):
		"return fai file as string"
		return self._get_str("zgoubi.fai")
	def save_fai(self, path):
		"save fai file to path"
		return self._save_file("zgoubi.fai", path)

	def spn_fh(self):
		"return file handle for spn file"
		return self._get_fh("zgoubi.spn")
	def spn(self):
		"return spn file as string"
		return self._get_str("zgoubi.spn")
	def save_spn(self, path):
		"save spn file to path"
		return self._save_file("zgoubi.spn", path)
		
	def b_fai_fh(self):
		"return file handle for binary fai file"
		return self._get_fh("b_zgoubi.fai")
	def b_fai(self):
		"return binary fai file as string"
		return self._get_str("b_zgoubi.fai")
	def save_b_fai(self, path):
		"save binary fai file to path"
		return self._save_file("b_zgoubi.fai", path)
		
	def b_plt_fh(self):
		"return file handle for binary plt file"
		return self._get_fh("b_zgoubi.plt")
	def b_plt(self):
		"return binary plt file as string"
		return self._get_str("b_zgoubi.plt")
	def save_b_plt(self, path):
		"save binary plt file to path"
		return self._save_file("b_zgoubi.plt", path)

	def impdev_fh(self):
		"return file handle for impdev file"
		return self._get_fh("zgoubi.impdev.out")
	def impdev(self):
		"return impdev file as string"
		return self._get_str("zgoubi.impdev.out")
	def save_impdev(self, path):
		"save impdev file to path"
		return self._save_file("zgoubi.impdev.out", path)
	
	def opticsout_fh(self):
		"return file handle for optics out file"
		return self._get_fh("zgoubi.OPTICS.out")
	def opticsout(self):
		"return optics out file as string"
		return self._get_str("zgoubi.OPTICS.out")
	def save_opticsout(self, path):
		"save optics out file to path"
		return self._save_file("zgoubi.OPTICS.out", path)
	


	def _bad_float(self, text):
		"""A wrapper around float to deal with zgoubi output numbers like
		2.67116100-102 when it means 2.67116100e-102
		

		"""
		try:
			return float(text)
		except ValueError:
			assert('-' in text[1:])
			error = 'cant make float from "' + text + '". '
			error += 'assuming it to be zero'
			zlog.warn(error)
			return 0
		
		
	def get_all_bin(self, file='bplt'):
		
		if file == 'bplt':
			return io.read_file(os.path.join(self.rundir, 'b_zgoubi.plt'))
		elif file == 'bfai':
			try:
				return io.read_file(os.path.join(self.rundir, 'b_zgoubi.fai'))
			except io.OldFormatError:
				pass

			fh = self.b_fai_fh()
			file_len = os.path.getsize(os.path.join(self.rundir, 'b_zgoubi.fai'))
			head_len = 352
			chunk_len = 251
			assert((file_len-head_len)%chunk_len == 0), "File size does not seem right for a binary fai file. File length is %s"%file_len
		else:
			raise ValueError("get_all_bin() expects name to be 'bplt' or 'bfai'")

		warnings.warn("Support for reading zgoubi <6.0 formats will be removed in the future")
		plt_data = []
		for pos in xrange(head_len, file_len, chunk_len):
			fh.seek(pos)
			chunk = fh.read(chunk_len)
			particle = {}
			# in vim use to create this next bit of code 
			# :r!devtools/format_to_code.py devtools/fai.format

			fmt = '=icidddddddddddddddiiixxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxdi10s10s10sii'
			chunk_part = chunk[0: struct.calcsize(fmt)]
			bits = struct.unpack(fmt, chunk_part)
			particle['LET'] = bits[1]
			particle['IEX'] = bits[2]
			particle['D0'] = bits[3]
			particle['Y0'] = bits[4]
			particle['T0'] = bits[5]
			particle['Z0'] = bits[6]
			particle['P0'] = bits[7]
			particle['X0'] = bits[8]
			particle['D'] = bits[9]
			particle['Y'] = bits[11]
			particle['T'] = bits[12]
			particle['Z'] = bits[13]
			particle['P'] = bits[14]
			particle['S'] = bits[15]
			particle['tof'] = bits[16]
			particle['KE'] = bits[17]
			particle['ID'] = bits[18]
			particle['IREP'] = bits[19]
			particle['SORT'] = bits[20]
			particle['BORO'] = bits[21]
			particle['PASS'] = bits[22]
			particle['element_type'] = bits[23].strip()
			particle['element_label1'] = bits[24].strip()
			particle['element_label2'] = bits[25].strip()
			particle['NOEL'] = bits[26]

			plt_data.append(particle)
		return plt_data
		
	def get_all(self, file='plt', id=None):
		"""Read all the data out of the file.
		Set file to plt or fai
		if using the new headered data formats will return a numpy array with named columns, otherwise returns a list of dicts. each dict has the particle coordinates at a point
		"""
		#if (isinstance(file_nh, str)):
		#	fh = open(file_nh)
		#elif (isinstance(file_nh, file)):
		#	fh = file_nh
		#else:
		#	error = "get_all(file_nh), file_nh should be string or file handle. "
		#	error += "It actually was"+ str(type(file_nh))
		#	raise TypeError, error

		if file == 'plt':
			fh = self.plt_fh()
		elif file == 'fai':
			fh = self.fai_fh()
		elif file == 'spn':
			fh = self.spn_fh()
		elif file == 'bfai':
			return self.get_all_bin(file=file)
		else:
			#open previously saved file
			fh = open(file)
			#decide what type of file we have, look at extension
			fileext = os.path.splitext(file)[1]
			if fileext == '.plt':
				file = 'plt'
			elif fileext == '.fai':
				file = 'fai'
			elif fileext == '.spn':
				file = 'spn'
			else:
				raise ValueError("get_all() expects '.plt', '.fai' or '.spn' extension")
			

		header = list(read_n_lines(fh, 4))
		#print "header", header
		#Versions of Zgoubi up to 236 had just '...' in lines 3 and 4 of the fai file 
		test_version = header[2][0:3]
		old_format = False
		if test_version == '...':
			old_format = True
			warnings.warn("Support for reading zgoubi <6.0 formats will be removed in the future")
		else:
			return io.read_file(fh)

		#for n,x in enumerate([9,3,5,10,8]):
		#   for i in xrange(x):
		#       print '"['+str(n)+','+str(i)+']",',


		l0_re =     trans_to_regex('%c +%d +%f +%f +%f +%f +%f +%f +%f')
		l1_re =     trans_to_regex('%f +%f +%f')
		l2_re =     trans_to_regex('%f +%f +%f +%f +%f')
		l3_re =     trans_to_regex('%d +%d +%f +%f +%f +%f +%f +%f +%f +%f')
		l3_re_plt = trans_to_regex('%d +%d +%d +%f +%f +%f +%f +%f +%f +%f')
		l4_re =     trans_to_regex('%f +%d %8c %8c%8c +%d')
		l4_re_plt = trans_to_regex('%f +%f +%f +%f +%d %8c %8c%8c +%d')

		l0_spn =    trans_to_regex('%d +%f +%f +%f +%f +%f +%f +%f +%f')
		l1_spn =    trans_to_regex('%f +%d +%d +%d %8c +%d')

		#no of lines per fai, plt or spn data point
		if file != 'spn':
			if old_format:
				n_lines = 5
			else:
				n_lines = 1
		else:
			n_lines = 2

		plt_data = []
		for lines in yield_n_lines(fh, n_lines):
			
			if id is not None:
				try:
					if lines[3].split()[1] != id:
						continue # escape quickly if this is not a point we are interested in
				except:
					print lines
					continue

			if len(lines) == n_lines: # dont go through this if we have just few dangling lines on the end of the file
				# fill dictionary with the bits that have been identified
				#
				p = {}
				# for compatability with new routines that expect numpy arrays
				if file=="fai" and old_format:
					p = numpy.zeros([1], dtype=[('element_label2', 'S8'), ('element_label1', 'S8'), ('BORO', '<f8'), ('X0', '<f8'), ('tof', '<f8'), ('SORT', '<f8'), ('IEX', '<i8'), ('P0', '<f8'), ('D', '<f8'), ('D-1', '<f8'), ('IREP', '<i8'), ('P', '<f8'), ('S', '<f8'), ('T0', '<f8'), ('PASS', '<i8'), ('Y', '<f8'), ('X', '<f8'), ('Z', '<f8'), ('ID', '<f8'), ('KE', '<f8'), ('Z0', '<f8'), ('element_type', 'S8'), ('LET', 'S1'), ('Y0', '<f8'), ('D0', '<f8'), ('D0-1', '<f8'), ('T', '<f8')])
				if file=="plt" and old_format:
					p = numpy.zeros([1], dtype=[('element_label2', 'S8'), ('element_label1', 'S8'), ('KART', '<f8'), ('BORO', '<f8'), ('X0', '<f8'), ('tof', '<f8'), ('By', '<f8'), ('Bz', '<f8'), ('SORT', '<f8'), ('IEX', '<i8'), ('P0', '<f8'), ('D', '<f8'), ('IREP', '<i8'), ('P', '<f8'), ('S', '<f8'), ('T0', '<f8'), ('PASS', '<i8'), ('Y', '<f8'), ('X', '<f8'), ('Z', '<f8'), ('ID', '<f8'), ('D-1', '<f8'), ('Z0', '<f8'), ('element_type', 'S8'), ('LET', 'S1'), ('D0-1', '<f8'), ('Y0', '<f8'), ('D0', '<f8'), ('T', '<f8')])

				if old_format:

					# line 0
					if file == 'spn':
						l0 = scanf(lines[0], l0_spn)
						p['SPIN_X'] = float(l0[5])
						p['SPIN_Y'] = float(l0[6])
						p['SPIN_Z'] = float(l0[7])
					else:
						l0 = scanf(lines[0], l0_re)
						p['LET'] = l0[0] 
						p['IEX'] = int(l0[1]) #flag
						p['D0'] = float(l0[2]) #initial D-1
						p['D0-1'] = float(l0[2]) #initial D-1
						p['Y0'] = self._bad_float(l0[3]) #initial Y
						p['T0'] = self._bad_float(l0[4]) #initial T (remember that T is dY/dX not time)
						p['Z0'] = float(l0[5]) #initial Z
						p['P0'] = float(l0[6]) #initial P
						p['X0'] = self._bad_float(float(l0[7])) # path length at origin of structure ?

					# line 1
					if file == 'spn':
						l1 = scanf(lines[1], l1_spn)
						p['PASS'] = int(l1[3])
					else:
						l1 = scanf(lines[1], l1_re)
						p['D'] = float(l1[0]) #D-1
						p['D-1'] = float(l1[0]) #D-1
						p['Y'] = self._bad_float(l1[1]) #current corrods
						p['T'] = self._bad_float(l1[2])

					#line 2
					if file != 'spn':
						l2 = scanf(lines[2], l2_re)
						p['Z'] = self._bad_float(l2[0])
						p['P'] = self._bad_float(l2[1])
						p['S'] = self._bad_float(l2[2])
						p['tof'] = self._bad_float(l2[3]) # time of flight
						if file == 'plt': # a strange bug prehaps. anyways this makes it so that you get microseconds
							p['tof'] = p['tof'] / 1e5
						if file == 'fai':
							p['KE'] = self._bad_float(l2[4]) # Kinetic energy
					
					#line 3
					# see 20080501 in lab book
					if file == 'plt':
						l3 = scanf(lines[3], l3_re_plt)
						p['KART'] = int(l3[0])
						p['ID'] = l3[1]
						p['IREP'] = l3[2]
						p['SORT'] = l3[3]
						p['X'] = float(l3[4])
						p['By'] = float(l3[6])
						p['Bz'] = float(l3[7])
					elif file == 'fai':
						l3 = scanf(lines[3], l3_re)
						p['ID'] = l3[0]
						p['IREP'] = l3[1]
						p['SORT'] = l3[2]
						p['X'] = -1 # X is not defined in FAI file, particle is always at end of the element

					#if  particle['X'] > 1e18: # this comes out as a silly giant number in some cases, maybe due to FAICNL or something
					#	particle['X'] = 0 

					#line4
					if file == 'plt':
						l4 = scanf(lines[4], l4_re_plt)
						p['BORO'] = float(l4[3])
						p['PASS'] = int(l4[4])
						p['element_type'] = l4[5].strip()
						p['element_label1'] = l4[6].strip()
						p['element_label2'] = l4[7].strip()
					elif file == 'fai':
						l4 = scanf(lines[4], l4_re)
						p['BORO'] = float(l4[0])
						p['PASS'] = int(l4[1])
						p['element_type'] = l4[2].strip()
						p['element_label1'] = l4[3].strip()
						p['element_label2'] = l4[4].strip()
				else:
					#new format - first separate the strings at end of l0 from the rest
					l0 = lines[0].split("\'")
					#remove spurious ' ' strings
					for elem in l0:
						if elem == ' ':
							l0.remove(elem)
					l0_stringpart = l0[1:-1]
					#split the numerical part of l0
					l0 = l0[0].split()

					p['IEX'] = int(l0[0])
					p['D0'] = float(l0[1]) #initial D-1
					p['Y0'] = self._bad_float(l0[2]) #initial Y
					p['T0'] = self._bad_float(l0[3]) #initial T (remember that T is dY/dX not time)
					p['Z0'] = float(l0[4]) #initial Z
					p['P0'] = float(l0[5]) #initial P
					p['X0'] = self._bad_float(float(l0[6])) # path length at origin of structure ?
					p['tof0'] = self._bad_float(l0[7]) # time of flight
					p['D'] = float(l0[8])
					p['Y'] = self._bad_float(l0[9])
					p['T'] = self._bad_float(l0[10])
					p['Z'] = self._bad_float(l0[11])
					p['P'] = self._bad_float(l0[12])
					p['S'] = self._bad_float(l0[13])
					p['tof'] = self._bad_float(l0[14]) # time of flight in microseconds
					if file == 'fai':
						p['KE'] = self._bad_float(l0[15]) # Kinetic energy
						p['E'] = self._bad_float(l0[16]) # Total energy
						p['ID'] = int(l0[17])
						p['IREP'] = int(l0[18])
						p['SORT'] = float(l0[19])
						p['AMQ'] = l0[20:25]
						p['RET'] = float(l0[25])
						p['DPR'] = float(l0[26])
						p['PS'] = float(l0[27])
						p['BORO'] = float(l0[28])
						p['PASS'] = int(l0[29])
						p['NOEL'] = int(l0[30])
						p['element_type'] = l0_stringpart[0].strip('\'')
						p['element_label1'] = l0_stringpart[1].strip('\'')
						p['element_label2'] = l0_stringpart[2].strip('\'')
						p['LET'] = l0_stringpart[3].strip('\'')
					elif file == 'plt':
						p['beta'] = self._bad_float(l0[15]) # Relativistic beta
						p['DS'] = self._bad_float(l0[16]) # Relativistic beta
						p['ID'] = int(l0[18])
						p['IREP'] = int(l0[19])
						p['SORT'] = float(l0[20])
						p['X'] = float(l0[21])
						p['Bx'] = float(l0[22])
						p['By'] = float(l0[23])
						p['Bz'] = float(l0[24])
						p['RET'] = float(l0[25])
						p['DPR'] = float(l0[26])
						p['PS'] = float(l0[27])
						p['BORO'] = float(l0[39])
						p['PASS'] = int(l0[40])
						p['NOEL'] = int(l0[41])
						p['element_type'] = l0_stringpart[0].strip()
						p['element_label1'] = l0_stringpart[1].strip()
						p['element_label2'] = l0_stringpart[2].strip()
						p['LET'] = l0_stringpart[3].strip()

				
				plt_data.append(p)
		if (file=="fai" or file=="plt") and old_format:
			plt_data = numpy.array(plt_data, dtype=plt_data[0].dtype)
			plt_data =  plt_data.reshape([-1])
		return plt_data


	def get_track(self, file, coord_list, multi_list=None):
		"""
		returns a list of coordinates specified, and multiply them by the multiplier if not none
		eg get_trac('plt', ['X','Y','element_label'],[0.01,0.01,None])
		
		If all the columns requested are numerical, and new headered data formats are being used then this function will return a numpy array
		"""
		#FIXME can probably give a rec array for mixed case. numpy is a requirement these days
		alldata = self.get_all(file)
		#check if we are using the new zgoubi.io version
		if type(alldata) == type(numpy.zeros(0)):
			#coords = numpy.zeros([alldata.size, len(coord_list)])
			
			# is set logic to see what best array type is
			# mixed arrays might break some old code, so fall back to dicts
			dtypes = set([alldata[c].dtype for c in coord_list])
			if not dtypes - set([numpy.dtype('i')]):
				# if there is nothing but ints, use int
				best_type = 'i'
			elif not dtypes - set([numpy.dtype('i'), numpy.dtype('f'), numpy.dtype('d')]):
				# if there is nothing but ints and float/double, use double
				best_type = 'd'
			else:
				# otherwise fall back to dicts
				best_type = "mixed"


			if best_type != "mixed":
				# all cols numeric, so give a fast numpy array
				coords = numpy.zeros([alldata.size, len(coord_list)], dtype=(best_type))

				for n, c in enumerate(coord_list):
					coords[:, n] = alldata[c]
					if multi_list:
						if multi_list[n]:
							coords[:, n] *= multi_list[n]
				return coords


		coords = []
		for p in alldata:
			this_coord = []
			for n, c in enumerate(coord_list):
				if multi_list is None:
					this_coord.append(p[c])
				elif multi_list[n] is None:
					this_coord.append(p[c])
				else:
					this_coord.append(p[c] * multi_list[n])
			coords.append(this_coord)
		return coords
	
	def loss_summary(self, coords=None, file='plt'):
		"""Returns False if no losses, otherwise returns a summery of losses
		::
		
			loss = res.loss_summary(file='plt')
			#or
			all = res.get_all('plt')
			loss = res.loss_summary(all) # if you already have got the coordinates
			
		"""
		if coords is None:
			coords = self.get_all(file)
	
		if all(coords["IEX"] == 1):
			return False
		
		loss_types = {-1:"the trajectory happened to wander outside the limits of a field map",
		              -2:"too many integration steps in an optical element",
		              -3:"deviation happened to exceed pi/2in an optical element",
		              -4:"stopped by walls (procedures CHAMBR, COLLIMA)",
		              -5:"too many iterations in subroutine DEPLA",
		              -6:"energy loss exceeds particle energy",
		              -7:"field discontinuities larger than 50% wthin a field map",
		              -8:"reached field limit in an optical element",
		              }
		loss_res = {}
		for iexval in loss_types.keys():
			lossnum = (coords["IEX"] == iexval).sum()
			if lossnum >= 1:
				loss_res[loss_types[iexval]] = lossnum
		return loss_res




	def get_bunch(self, file, end_label=None, old_bunch=None, drop_lost=True):
		""""Get back a bunch object from the fai file. It is recommended that you put a MARKER before the last FAISCNL, and pass its label as end_label, so that only the bunch at the final position will be returned. All but the final lap is ignored automatically.
		Optionally the an old_bunch can be passed to the function, its mass and charge will be copyed to the new bunch.
		"""
		try:
			all_c = self.get_all(file)
		except IOError:
			zlog.warn("Could not read %s. returning empty bunch", file)
			empty_bunch = zgoubi.bunch.Bunch(nparticles=0, rigidity=0)
			if old_bunch is not None:
				empty_bunch.mass = old_bunch.mass
				empty_bunch.charge = old_bunch.charge
			return empty_bunch
		except EmptyFileError:
			zlog.warn("%s empty. returning empty bunch", file)
			empty_bunch = zgoubi.bunch.Bunch(nparticles=0, rigidity=0)
			if old_bunch is not None:
				empty_bunch.mass = old_bunch.mass
				empty_bunch.charge = old_bunch.charge
			return empty_bunch

		loss_sum = self.loss_summary(all_c)
		if loss_sum:
			for k, v in loss_sum.items():
				zlog.warn("%s particles lost: %s" % (v, k))

			if drop_lost:
				all_c = all_c[all_c['IEX'] == 1]



		if not type(all_c) == type(numpy.zeros(0)):
			raise OldFormatError("get_bunch() only works with the new fai format")
		
		# select only the particles that made it to the last lap
		last_lap = all_c[all_c['PASS'] == all_c['PASS'].max()]
		# also select only particles at FAISTORE with matching end_label
		if end_label:
			end_label = end_label.ljust(last_lap.dtype['element_label1'].itemsize) # pad to match zgoubi, as of Zgoubi SVN r290 this has changed from 8 to 10
			last_lap = last_lap[last_lap['element_label1'] == end_label]

		#print last_lap[:10]['BORO']
		#print last_lap[:10]['D-1']

		# Create a new bunch and copy the values arcoss (with conversion to SI)
		last_bunch = zgoubi.bunch.Bunch(nparticles=last_lap.size, rigidity=last_lap[0]['BORO']/1000)
		if old_bunch is not None:
			last_bunch.mass = old_bunch.mass
			last_bunch.charge = old_bunch.charge
		particles = last_bunch.particles()
#		print last_lap.dtype
		particles['Y'] = last_lap['Y'] /100
		particles['T'] = last_lap['T'] /1000
		particles['Z'] = last_lap['Z'] /100
		particles['P'] = last_lap['P'] /1000
		particles['S'] = last_lap['S'] /100
		particles['D'] = last_lap['D-1'] +1

		return last_bunch


	def get_extremes(self, file, element_label=None, coord='Y', id=None):
		"""
		Get the max and min positions of a certain coordinate, in a certain element
		"""
		# FIXME could be done better with numpy
		all = self.get_all(file, id=id)
		if element_label is None:
			el_points = all
		else:
			el_points = [p for p in all if (p['element_label1'] == element_label or p['element_label2'] == element_label)]

		
		if len(el_points) == 0:
			print "no points for particle id:", id, "element_label", element_label
			raise NoTrackError
			
		
		mini = min([p[coord] for p in el_points])
		maxi = max([p[coord] for p in el_points])

		return mini, maxi

	def list_particles(self, file):
		if file == 'plt':
			fh = self.plt_fh()
		elif file == 'fai':
			fh = self.fai_fh()
		else:
			raise ValueError("get_all() expects name to be 'plt' or 'fai'")

		#ignore header lines
		dummy = list(read_n_lines(fh, 4))
		particle_ids = set()
		for lines in yield_n_lines(fh, 5):
			if len(lines) == 5: # dont go through this if we have just few dangling lines on the end of the file
				#bits = [] # not computerscience bits. just the bits of data on each line
				#bits.append(line.split()) # split them by whitespace
				#particle_ids.add(int(bits[3][1]))
				particle_ids.add(int(lines[3].split()[1]))
		return particle_ids

	def in_bounds(self, file, element_label, min_bound, max_bound, coord='Y', verbose=False, id=None):
		"""
		check if particle exceeded bounds in this element. pass bounds in zgoubi's default unit for that coordinate
		"""
		if id is not None:
			part_list = [id]
		else:
			part_list = self.list_particles(file)

		
		for id in part_list:
			try:
				track_min, track_max = self.get_extremes(file, element_label, coord, id=id)
			except NoTrackError:
				print "No Track"
				return True
			in_b = True
			if track_min <= min_bound:
				if verbose: print track_min, "<=", min_bound, "hit lower bound, in element", element_label
				in_b = False
				
			if track_max >= max_bound:
				if verbose: print track_max, ">=", max_bound, "hit upper bound, in element", element_label
				in_b = False

		return in_b

	def check_bounds(self, file, min_bounds, max_bounds, coord='Y', part_ids=None, assume_in_order=False):
		"""read whole file. for each chunk, check if it is an element with bounds, if yes check if it is inside
		record crashes, by particle id
		min_bounds and max_bounds should be of the form {'label': bound, 'label': bound}
		returns numpy arrays
		note: lost particles dont crash
		"""
		if part_ids is None:
			particle_ids = self.list_particles(file)
		else:
			particle_ids = part_ids
		if file == 'plt':
			fh = self.plt_fh()
		elif file == 'fai':
			fh = self.fai_fh()
		else:
			raise ValueError("get_all() expects name to be 'plt' or 'fai'")

		header = list(read_n_lines(fh, 4))
		
		crashes = numpy.zeros(len(particle_ids), dtype=int)
		not_crashes = numpy.zeros(len(particle_ids), dtype=int)
		laps = numpy.zeros(len(particle_ids), dtype=int)
		crash_lap = numpy.zeros(len(particle_ids), dtype=int)
		crash_lap += 100000
		
		for lines in yield_n_lines(fh, 5):
			if len(lines) != 5:
				continue
			bits = {} # not computerscience bits. just the bits of data on each line
			particle = {}
			#for line in lines:
			#	bits.append(line.split()) # split them by whitespace
			bits[1] = lines[1].split()
			bits[3] = lines[3].split()
			bits[4] = lines[4].split()
			
			particle['ID'] = bits[3][1]
			particle['Y'] = float(bits[1][1])
			particle['T'] = float(bits[1][2])
			
			if assume_in_order:
				if crashes[int(particle['ID'])-1] > 0:
					continue
			
			#the last line of each chunk is a pain.
			#lets find the first none numeric element and work from that
			for n, test_string in enumerate(bits[4]):
				try:
					float(test_string)
				except ValueError:
					first_non_numeric = n
					break

			#particle['BORO'] = bits[4][first_non_numeric -2]
			particle['PASS'] = bits[4][first_non_numeric -1]
			#print bits[4]
			laps[int(particle['ID'])-1] = max(laps[int(particle['ID'])-1], int(particle['PASS'])) # store the biggest lap seen
			#particle['element_type'] = bits[4][first_non_numeric]
			#set these to empty, and fill them in later if possible
			particle['element_label1'] = ""
			particle['element_label2'] = ""
			#particle['element_number'] = bits[4][len(bits[4]) -1]
			if (len(bits[4]) - first_non_numeric) == 4: # 2 labels
				particle['element_label1'] = bits[4][first_non_numeric + 1]
				particle['element_label2'] = bits[4][first_non_numeric + 2]
			if (len(bits[4]) - first_non_numeric) == 3: # 1 labels
				particle['element_label1'] = bits[4][first_non_numeric + 1]
			
			has_crashed = False
			for k, v in min_bounds.items():
				if particle['element_label1'] == k or particle['element_label2'] == k:
					if particle['Y'] <= v:
						#crash
						has_crashed = True
			for k, v in max_bounds.items():
				if particle['element_label1'] == k or particle['element_label2'] == k:
					if particle['Y'] >= v:
						#crash
						has_crashed = True
			if has_crashed:
				# zgoubi counts from 1, we count from zero
				crashes[int(particle['ID'])-1] += 1
				crash_lap[int(particle['ID'])-1] = min(crash_lap[int(particle['ID'])-1], int(particle['PASS'])) # record the first crash
			else:
				not_crashes[int(particle['ID'])-1] += 1

		return crashes, not_crashes, laps, crash_lap
	
	def parse_matrix(self):
		"""Parse data from output of MATRIX element

		Used by get_tune(), get_transfer_matrix() and get_twiss_parameters()
		"""

		has_object5 = 'OBJET5' in self.element_types
		has_matrix = 'MATRIX' in self.element_types
				
		if not (has_object5 and has_matrix):
			raise BadLineError("beamline need to have an OBJET with kobj=5 (OBJET5), and a MATRIX element with IORD=1 and IFOC>10 to get tune")

		found_output = False
		in_matrix = False

		parsed_info = dict(matrix1=None, twiss=None, tune=(-1, -1))
		tp = twiss_param_array()
		matrix_lines = []
		
		fh = self.res_fh()
		while True:
			try: line = fh.next()
			except StopIteration: break
			if line.startswith("******"): # find start of an output block
				found_output = True
				try: element_line = fh.next()
				except StopIteration: break
				if element_line.strip() == "": continue

				# find element type
				if "Keyword, label(s)" in element_line: # changed somewhere between svn261 and svn312
					element_type = element_line.split()[4] #new
				else:
					element_type = element_line.split()[1] #old

				line = element_line
				in_matrix = element_type == "MATRIX"

			if not found_output: continue
			if in_matrix:
				matrix_lines.append(line.strip())

		if len(matrix_lines) == 0:
			raise BadLineError("Could not find MATRIX output in res file. Maybe beam lost.")

		#print "\n".join(matrix_lines)

		for n, line in enumerate(matrix_lines):
			if line == "TRANSFER  MATRIX  ORDRE  1  (MKSA units)":
				matrix1 = numpy.zeros([6, 6])
				for x in xrange(6):
					matrix1[x] = matrix_lines[x+n+2].split()
				parsed_info['matrix1'] = matrix1
			if line == "Beam  matrix  (beta/-alpha/-alpha/gamma) and  periodic  dispersion  (MKSA units)":
				twissmat = numpy.zeros([6, 6])
				for x in xrange(6):
					row = []
					for rowe in matrix_lines[x+n+2].split():
						try:
							row.append(float(rowe))
						except ValueError:
							zlog.warn("Could not convert %s to float"%rowe)
							row.append(numpy.nan)
					if len(row) != 6:
						zlog.warn("Could not convert get 6 values from row: %s"%matrix_lines[x+n+2] )
						row = numpy.nan
					twissmat[x] = row
				tp['beta_y'] = twissmat[0, 0]
				tp['alpha_y'] = -twissmat[1, 0]
				tp['gamma_y'] = twissmat[1, 1]
				tp['disp_y'] = twissmat[0, 5]
				tp['disp_py'] = twissmat[1, 5]

				tp['beta_z'] = twissmat[2, 2]
				tp['alpha_z'] = -twissmat[3, 2]
				tp['gamma_z'] = twissmat[3, 3]
				tp['disp_z'] = twissmat[2, 5]
				tp['disp_pz'] = twissmat[3, 5]
				parsed_info['twiss'] = tp
			if line.startswith("NU_Y"):
				#print "\n# ".join(matrix_lines)
				nu_y = nu_z = -1
				bits = line.split()
				if bits[2] == "undefined":
					zlog.error("Horizontal tune is undefined: %s"%line)
				if bits[5] == "undefined":
					zlog.error("Vertical tune is undefined: %s"%line)
				try:
					nu_y = float(bits[2])
				except ValueError:
					pass
				try:
					nu_z = float(bits[5])
				except ValueError:
					pass

				parsed_info['tune'] = (nu_y, nu_z)

		return parsed_info


	def get_tune(self):
		"""Returns a tuple (NU_Y, NU_Z) of the tunes.
		Needs a beam line is an OBJET type 5, and a MATRIX element.


		"""
		tune = self.parse_matrix()["tune"]

		if tune[0] == -1 or tune[1] == -1:
			zlog.error("could not get tune res file setting to -1")

		return tune


	def get_transfer_matrix(self):
		"""Returns a transfer matrix of the line in (MKSA units).
		Needs a beam line is an OBJET type 5, and a MATRIX element.

		"""
		tm = self.parse_matrix()["matrix1"]

		if tm is None:
			zlog.error("could not get matrix from res file")

		return tm


	def get_twiss_parameters(self):
		"""Returns a tuple (beta_y, alpha_y, gamma_y, disp_y, disp_py, beta_z, alpha_z, gamma_z, disp_z, disp_pz) from 
			the twiss parameter matrix.

		The results are stored a in structured numpy array containing the following elements 

		#Outputs
		beta_y, alpha_y, gamma_y: horizontal Twiss parameters
		beta_y, alpha_y, gamma_y: vertical Twiss parameters
		disp_y, disp_py : horizontal dispersion and dispersion-prime
		disp_z, disp_pz : vertical dispersion and dispersion-prime

		Needs a beam line is an OBJET type 5, and a MATRIX element.

		"""
		tp = self.parse_matrix()["twiss"]

		if tp is None:
			zlog.error("could not get twiss parameters from res file")

		return tp

	def show_particle_info(self):
		"show the particle info, a good check of energies, mass etc"
		in_particle = False
		past_input = False
		fh = self.res_fh()
		while True:
			line = fh.readline()
			if line.startswith("***"):
				if in_particle: # going into a new section, so can break
					break
				past_input = True
				section = fh.readline().strip().split()[-1]
				if section == "PARTICUL": #entering particle
					in_particle = True
				continue

			if not past_input:
				continue
			if in_particle:
				if line.strip().startswith("I, AMQ(1,I)"):
					break
				print line,


	def test_rebelote(self):
		"""Return true if end of REBELOTE procedure reported
		Needs a REBELOTE element.

		"""
		#has_reb = False
		#for e in self.line.elements():
		#	t = str(type(e)).split("'")[1].rpartition(".")[2]
		#	if t == 'REBELOTE':
		#		has_reb = True
		has_reb = 'REBELOTE' in self.element_types
		if not has_reb:
			raise BadLineError, "beamline need to have a REBELOTE for this function"

		for line in self.res_fh():
			if "End  of  'REBELOTE'  procedure" in line:
				print "REBELOTE completed"
				return True
		return False

	def run_success(self):
		"""Checks that zgoubi completed

		"""
		for line in self.res_fh():
			if "MAIN PROGRAM : Execution ended upon key  END" in line or "ZGOUBI RUN COMPLETED" in line:
				return True
			elif "Execution ended normally, upon keyword END or FIN" in line:
				return True
		return False
		





