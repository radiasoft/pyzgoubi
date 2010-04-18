#!/usr/bin/env python
# -*- coding: utf-8 -*-

#		pyzgoubi - python interface to zgoubi
#		Copyright 2008 Sam Tygier <Sam.Tygier@hep.manchester.ac.uk>
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


from __future__ import division
from string import ascii_letters, digits
from math import *
import tempfile
import shutil
import os
import sys
import re
import struct
from glob import glob
import copy
import threading
import Queue
import subprocess

try:
	import numpy
except ImportError:
	print "could not import numpy, some functions will give errors about numpy not being defined"
try:
	from operator import itemgetter
except ImportError:
	print "please use python 2.5 or newer"
	sys.exit(1)
try:
	import cairo
except ImportError:
	pass
#	print "cairo not found, no plotting available"

from zgoubi.constants import *
from zgoubi.exceptions import *
import zgoubi.io as io
import zgoubi.bunch


from zgoubi.settings import zgoubi_settings

sys.setcheckinterval(10000)

zgoubi_module_path = os.path.dirname( os.path.realpath( __file__ ) )
zgoubi_path = zgoubi_settings['zgoubi_path']

pyzgoubi_egg_path = glob(os.path.dirname(zgoubi_module_path) + "/pyzgoubi*egg-info")

if (not zgoubi_module_path.startswith("/")) and os.name == 'posix' :
	# windows style paths, but os.name is 'posix', so need to mess with paths by hand
	bits = zgoubi_module_path.split('\\')[:-3]+['share','pyzgoubi','definitions']
	static_defs = '\\'.join(bits+['static_defs.py'])
	simple_defs = '\\'.join(bits+['simple_elements.defs'])
else:
	static_defs = os.path.join(zgoubi_module_path,'..','..','..','..' ,"share","pyzgoubi","definitions", "static_defs.py")
	simple_defs = os.path.join(zgoubi_module_path,'..','..','..','..' ,"share","pyzgoubi","definitions", "simple_elements.defs")


#definitions of elements
home = os.path.expanduser('~')
config_dir = os.path.join(home, ".pyzgoubi")
if not os.path.exists(config_dir):
	os.mkdir(config_dir)

def_cache_dir = os.path.join(config_dir, "def_cache")
if not os.path.exists(def_cache_dir):
	os.mkdir(def_cache_dir)

compiled_defs_path = os.path.join(config_dir, "def_cache","defs.py")

definitions_paths = [os.path.join(config_dir, x) for x in os.listdir(config_dir) if x.endswith('.defs')]
definitions_paths.append(simple_defs)
definitions_paths += zgoubi_settings['extra_defs_files']
nl = '\n'

#print "definitions at", definitions_paths
#if needed recompile defs
if not os.path.exists(compiled_defs_path):
	print "no compiled defs. compiling"
	need_def_compile = True
else:
	need_def_compile = False
	for f in definitions_paths:
#		print "checking ", f
		if os.path.exists(f) and os.path.getmtime(f) >= os.path.getmtime(compiled_defs_path):
			print "need to recompile", f
			need_def_compile = True
	if not need_def_compile: # also need to recompile after a install or update
		for f in pyzgoubi_egg_path:
			if os.path.exists(f) and os.path.getmtime(f) >= os.path.getmtime(compiled_defs_path):
				print "pyzgoubi first run, compiling defs"
				need_def_compile = True
				break
if need_def_compile:
	from zgoubi import makedefs
	print "Compiling definitions"
	makedefs.make_element_classes(definitions_paths, compiled_defs_path)




#force rebuild everytime
#import zgoubi_makedefs
#zgoubi_makedefs.make_element_classes(definitions_paths, compiled_defs)




def yield_n_lines(fh, n):
	"yields n lines at a time"

	lines = []
	for line in fh:
		lines.append(line)
		if (len(lines) == n):
			yield lines
			lines = []

	yield lines


def read_n_lines(fh, n):
	"reads n lines at a time"

	lines = []
	for x in xrange(n):
		lines.append(fh.readline())
	return lines

def trans_to_regex(fmt):
	"Transform a printf style format in to a regular expression"
	fmt = fmt.replace('%c', '(.)')
	fmt = re.sub('%(\d+)c', r'(.{\1})', fmt)
	fmt = fmt.replace('%d', '([-+]?\d+)')
	fmt = fmt.replace('%e', '([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%E', '([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%f', '([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%g', '([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
	fmt = fmt.replace('%i', '([-+]?(?:0[xX][\dA-Fa-f]+|0[0-7]*|\d+))')
	fmt = fmt.replace('%o', '(0[0-7]*)')
	fmt = fmt.replace('%s', '(\S+)')
	fmt = fmt.replace('%u', '(\d+)')
	fmt = fmt.replace('%x', '(0[xX][\dA-Fa-f]+")')
	fmt = fmt.replace('%X', '(0[xX][\dA-Fa-f]+")')
	#print "fmt=", fmt
	return re.compile(fmt)

def scanf(str, fmt):
	"something like scanf, used to read the fai and plt files in Results.get_all()"
	if (type(fmt) == type(re.compile(''))):
		res = fmt.search(str)
	else:
		res = trans_to_regex(fmt).search(str)
	if res == None:
		return None
	else:
		return res.groups()

# a base class for all the beam line objects


class zgoubi_element(object):
	def __init__(self):
		pass
	
	def set(self, *dsettings, **settings):
		
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
				raise ValueError,  "no such param: '" + str(key) + "' In element " + self._zgoubi_name
	def get(self,key):
		return self._params[key]

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
		return list(self._params)
			
try:
	execfile(static_defs)
except IOError:
	print "Could not load static element definitions, maybe you are running pyzgoubi from the source directory"
	sys.exit(1)
execfile(compiled_defs_path)

class Line(object):
	def __init__(self, name):
		self.element_list = []
		self.name = name
		self.tmp_folders = [] # keep track of old tmp folders
		self.no_more_xterm = False
		self.input_files = []
		
		self.shutil = shutil # need to keep a reference to shutil
							# otherwise the intepreter may have thrown it away
							# by the time the __del__() is run
		self.has_run = False
		self.full_line = False # has an OBJET, dont allow full lines to be added to each other
								# only a full line outputs its name into zgoubi.dat

	def __del__(self):
		self.clean()
	
	def __neg__(self):
		"return a reversed line"
		new_line = copy.copy(self)
		new_line.element_list = copy.copy(new_line.element_list)
		new_line.element_list.reverse()
		new_line.name = "-"+self.name
		return new_line

	def __str__(self, prefix=""):
		"For printing a line"
		out = self.name + "\n"
		for element in self.element_list:
			if isinstance(element, Line):
				out += element.__str__(prefix+" ")
			else:
				out +=  "%s %s %s %s\n" % (prefix, element._zgoubi_name, element.label1, element.label2)
		return out

	def elements(self):
		for element in self.element_list:
			try:
				# decend into sub lines
				for sub_element in element.elements():
					yield sub_element
			except AttributeError:
				yield element

	def add(self, *elements):
		"Add an elements to the line"
		for element in elements:
			self.element_list.append(element)
			try:
				if 'OBJET' in element._zgoubi_name:
					self.full_line = True
			except AttributeError:
				pass
	
	def full_tracking(self, enable=True):
		"""Enable full tracking on magnetic elements.
		This works by setting IL=2 for any element with an IL parameter.
		use line.full_tracking(False) to disable tracking

		"""
		for e in self.elements():
			#t = str(type(e)).split("'")[1] #get type, eg zgoubi22.QUADRUPO
			#print t

			# if this gives "RuntimeError: maximum recursion depth exceeded"
			#errors, then the object may not have a parameters dict.
			try:
				#print e.output()
				if enable == True:
					e.set(IL=2)
				elif enable == False:
					e.set(IL=0)
			#	print e.output()
			except ValueError:
				pass
			
	def remove_looping(self):
		"removes any REBELOTE elements from the line"
		self.element_list = [e for e in self.element_list if ("REBELOTE" not in str(type(e)).split("'")[1])]
				
		
	def output(self):
		"Generate the zgoubi.dat file, and return it as a string"
		out = ""
		if self.full_line:
			out = self.name + nl
		
		for element in self.element_list:
			out += element.output() + nl
		
		return out

		
	def run(self, xterm=False, tmp_prefix=zgoubi_settings['tmp_dir'], silence=False):
		"Run zgoubi on line. If break is true, stop after running zgoubi, and open an xterm for the user in the tmp dir. From here zpop can be run."
		orig_cwd = os.getcwd()
		tmpdir = tempfile.mkdtemp("zgoubi", prefix=tmp_prefix)
		self.tmpdir = tmpdir
	
		for file in self.input_files:
			src = os.path.join(orig_cwd, file)
			dst = os.path.join(tmpdir, os.path.basename(file))
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
		if silence:
			command += " > zgoubi.stdout 2> zgoubi.sdterr"
		#exe_result = os.system(command)
		z_proc = subprocess.Popen(command, shell=True, cwd=tmpdir)
		exe_result = z_proc.wait()

		if exe_result != 0:
			print "zgoubi failed to run"
			print "It returned:", exe_result
			if exe_result == 32512: print "check that fortran runtime libraries are installed"

		if (xterm and not self.no_more_xterm):
			print "Do you want an xterm? (y=yes/n=no/s=stop asking)"
			ans = raw_input()
			if ans.startswith('y'):
				subprocess.Popen("xterm", shell=True, cwd=tmpdir).wait()
			elif ans.startswith('s'):
				self.no_more_xterm = True
		
		#os.system('ls')
		
		self.res_file = tmpdir+"/zgoubi.res"
		self.dat_file = tmpdir+"/zgoubi.dat"
		self.plt_file = tmpdir+"/zgoubi.plt"
		self.fai_file = tmpdir+"/zgoubi.fai"
		self.spn_file = tmpdir+"/zgoubi.spn"
		#output = outfile.read()
		
		#os.chdir(orig_cwd)
		
		self.has_run = True	
		return Results(line=self, rundir=tmpdir)
	
	def track_bunch(self, bunch, **kwargs):
		"Track a bunch through a Line, and return the bunch. This function will uses the OBJET_bunch object, and so need needs a Line that does not already have a OBJET"
		if self.full_line:
			raise BadLineError("If line already has an OBJET use run()")

		#build a line with the bunch OBJET and segment we were passed
		new_line = Line("bunch_line")
		new_line.add(OBJET_bunch(bunch, binary=False))
		new_line.add(self)
		#mark the faiscnl that we are interested in
		new_line.add(MARKER("trackbun"))
		new_line.add(FAISCNL(FNAME='b_zgoubi.fai'))
		new_line.add(END())

		# run the line
		result = new_line.run(**kwargs)
		
		# return the track bunch
		return  result.get_bunch('bfai', end_label="trackbun", old_bunch=bunch)
		
	def track_bunch_mt(self, bunch, n_threads=4, **kwargs):
		in_q = Queue.Queue()
		out_q = Queue.Queue()

		def worker(in_q, out_q, work_line, name):
			while True:
				start_index, work_bunch = in_q.get()
				#print "Thread", name, "working"
				try:
					done_bunch = work_line.track_bunch(work_bunch, **kwargs)
				except:
					print "Exception in track_bunch() thread"
					out_q.put((sys.exc_info()))
				else:
					out_q.put((start_index, done_bunch.particles()))
				in_q.task_done()
				#print "Thread", name, "task done"

		for x in xrange(n_threads):
			t = threading.Thread(target=worker,
			                     kwargs={'in_q':in_q, 'out_q':out_q, 'work_line':self, 'name':x})
			t.setDaemon(True)
			t.start()
			#print "Created thread", x
		
		start_index = 0
		n_tasks = 0
		#print "Queuing work"
		for item in bunch.split_bunch(max_particles=10000, n_slices=n_threads):
			#print "Queue task", n_tasks
			in_q.put((start_index, item))
			start_index += len(item)
			n_tasks +=1
		#print "Work queued"
		#in_q.join()
		#print "Work done"

		final_bunch = zgoubi.bunch.Bunch(nparticles = len(bunch), rigidity=bunch.get_bunch_rigidity(), mass=bunch.mass, charge=bunch.charge)
		for x in xrange(n_tasks):
			#print "collecting task", x
			result = out_q.get()
			try:
				start_index, done_bunch = result
			except ValueError:
				import traceback, time
				time.sleep(2) # give the other threads a couple of seconds, to make output prettier
				print "Exception retrieved by main thread"
				print
				traceback.print_exception(*result)
				print
				#reraise error message
				raise result[0](result[1])
			out_q.task_done()
			final_bunch.particles()[start_index:start_index+len(done_bunch)] = done_bunch

		#print "all done"
		return final_bunch





	def clean(self):
		"clean up temp directories"
		for dir in self.tmp_folders:
		#	print "removing", dir
		#	print shutil
			self.shutil.rmtree(dir)
		
		self.tmp_folders = [] # and blank list
		
	def res(self):
		"return zgoubi.res as a string"
		if (not self.has_run): print "Line has not been run"
		fh = open(self.res_file)
		return fh.read()
		
	def dat(self):
		"return zgoubi.dat as a string"
		if (not self.has_run): print "Line has not been run"
		fh = open(self.dat_file)
		return fh.read()
		
	def plt(self):
		"return zgoubi.plt as a string"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.plt_file)
		return fh.read()
	
	def fai(self):
		"return zgoubi.fai as a string"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.fai_file)
		return fh.read()

	def spn(self):
		"return zgoubi.spn as a string"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.spn_file)
		return fh.read()

	def res_fh(self):
		"return zgoubi.res file handle"
		if (not self.has_run): print "Line has not been run"
		fh = open(self.res_file)
		return fh
		
	def dat_fh(self):
		"return zgoubi.dat file handle"
		if (not self.has_run): print "Line has not been run"
		fh = open(self.dat_file)
		return fh
		
	def plt_fh(self):
		"return zgoubi.plt file handle"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.plt_file)
		return fh
		
	def fai_fh(self):
		"return zgoubi.fai file handle"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.fai_file)
		return fh

	def spn_fh(self):
		"return zgoubi.spn file handle"
		if (not self.has_run): print "Line has not been run"
		fh= open(self.spn_file)
		return fh

	def save_res(self, path):
		"save zgoubi.res to path"
		if (not self.has_run): print "Line has not been run"
		shutil.copyfile(self.res_file, path)

	def save_dat(self, path):
		"save zgoubi.dat to path"
		if (not self.has_run): print "Line has not been run"
		shutil.copyfile(self.dat_file, path)
		
	def save_plt(self, path):
		"save zgoubi.plt to path"
		if (not self.has_run): print "Line has not been run"
		shutil.copyfile(self.plt_file, path)
		
	def save_fai(self, path):
		"save zgoubi.fai to path"
		if (not self.has_run): print "Line has not been run"
		shutil.copyfile(self.fai_file, path)

	def save_spn(self, path):
		"save zgoubi.spn to path"
		if (not self.has_run): print "Line has not been run"
		shutil.copyfile(self.spn_file, path)

	def add_input_files(self, file_paths=[], pattern=None):
		"""Add some extra input files to the directory where zgoubi is run.
		This is useful for field map files.
		file_paths must be an iterable, for example a list
		To add many files use a pattern eg:
		l.add_input_files(pattern="maps/*")

		"""

		if pattern != None:
			globbed_files = glob.glob(pattern)
			file_paths += globbed_files
		self.input_files += file_paths
		
	def replace(self, elementold, elementnew, select_index=0):
		"Replace an element in the line. setting select_index to n will replace the nth occurence of that item. If select index is not set, the first occurence is replaced"

		indices = []
		i = -1
		try:
			while 1:
				i = self.element_list.index(elementold, i+1)
				indices.append(i)
		except ValueError:
			pass
		index = indices[select_index]

		self.element_list.pop(index)
		self.element_list.insert(index,elementnew)

	def insert(self, index, *elements):
		"Insert elements into the line before position given by index"
		for element in elements:
			self.element_list.insert(index,element)

	def remove(self,index):
		"Remove element at index"
		self.element_list.pop(index)

	def find_elements(self, element):
		indices = []
		i = -1
		try:
			while 1:
				i = self.element_list.index(element, i+1)
				indices.append(i)
		except ValueError:
			pass	

		return indices
			
	
class Results(object):
	"""This class lets you analyse the results after running a line.

	"""
	def __init__(self, line=None, rundir=None):
		self.line = line
		self.rundir = rundir

	# generic functions for accessing files in results folder
	def _get_fh(self, f):
		"return f file handle"
		path = os.path.join(self.rundir, f)
		if not os.path.exists(path):
			raise IOError, "No file: %s in %s"%(f, self.rundir)
		fh = open(path)
		return fh
		
	def _get_str(self, f):
		"return f as a string"
		return self._get_fh(f).read()
	
	def _save_file(self,f, path):
		"save f to path"
		spath = os.path.join(self.rundir, f)
		if not os.path.exists(spath):
			raise IOError, "No file: %s in %s"%(f, self.rundir)
		shutil.copyfile(spath, path)
	
	#generate specific functions
	def res_fh(self): return self._get_fh("zgoubi.res")
	def res(self): return self._get_str("zgoubi.res")
	def save_res(self, path): return self._save_file("zgoubi.res", path)
		
	def plt_fh(self): return self._get_fh("zgoubi.plt")
	def plt(self): return self._get_str("zgoubi.plt")
	def save_plt(self, path): return self._save_file("zgoubi.plt", path)
		
	def dat_fh(self): return self._get_fh("zgoubi.dat")
	def dat(self): return self._get_str("zgoubi.dat")
	def save_dat(self, path): return self._save_file("zgoubi.dat", path)
		
	def fai_fh(self): return self._get_fh("zgoubi.fai")
	def fai(self): return self._get_str("zgoubi.fai")
	def save_fai(self, path): return self._save_file("zgoubi.fai", path)

	def spn_fh(self): return self._get_fh("zgoubi.spn")
	def spn(self): return self._get_str("zgoubi.spn")
	def save_spn(self, path): return self._save_file("zgoubi.spn", path)
		
	def b_fai_fh(self): return self._get_fh("b_zgoubi.fai")
	def b_fai(self): return self._get_str("b_zgoubi.fai")
	def save_b_fai(self, path): return self._save_file("b_zgoubi.fai", path)
		
	def b_plt_fh(self): return self._get_fh("b_zgoubi.plt")
	def b_plt(self): return self._get_str("b_zgoubi.plt")
	def save_b_plt(self, path): return self._save_file("b_zgoubi.plt", path)

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
			print error
			return 0
		
		
	def get_all_bin(self, file='bplt'):
		
		if(file == 'bplt'):
			return io.read_file(os.path.join(self.rundir,'b_zgoubi.plt'))
		elif(file == 'bfai'):
			try:
				return io.read_file(os.path.join(self.rundir,'b_zgoubi.fai'))
			except io.OldFormatError:
				pass

			fh = self.b_fai_fh()
			file_len = os.path.getsize(os.path.join(self.rundir,'b_zgoubi.fai'))
			head_len = 352
			chunk_len = 251
			assert((file_len-head_len)%chunk_len == 0), "File size does not seem right for a binary fai file. File length is %s"%file_len
		else:
			raise ValueError, "get_all_bin() expects name to be 'bplt' or 'bfai'"

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
			particle['LET']=bits[1]
			particle['IEX']=bits[2]
			particle['D0']=bits[3]
			particle['Y0']=bits[4]
			particle['T0']=bits[5]
			particle['Z0']=bits[6]
			particle['P0']=bits[7]
			particle['X0']=bits[8]
			particle['D']=bits[9]
			particle['Y']=bits[11]
			particle['T']=bits[12]
			particle['Z']=bits[13]
			particle['P']=bits[14]
			particle['S']=bits[15]
			particle['tof']=bits[16]
			particle['KE']=bits[17]
			particle['ID']=bits[18]
			particle['IREP']=bits[19]
			particle['SORT']=bits[20]
			particle['BORO']=bits[21]
			particle['PASS']=bits[22]
			particle['element_type']=bits[23].strip()
			particle['element_label1']=bits[24].strip()
			particle['element_label2']=bits[25].strip()
			particle['NOEL']=bits[26]

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

		if(file == 'plt'):
			fh = self.plt_fh()
		elif(file == 'fai'):
			fh = self.fai_fh()
		elif(file == 'spn'):
			fh = self.spn_fh()
		elif(file == 'bfai'):
			return self.get_all_bin(file=file)
		else:
			#open previously saved file
			fh = open(file)
			#decide what type of file we have, look at extension
			fileext = os.path.splitext(file)[1]
			if(fileext == '.plt'):
				file = 'plt'
			elif(fileext == '.fai'):
				file = 'fai'
			elif(fileext == '.spn'):
				file = 'spn'
			else:
				raise ValueError, "get_all() expects '.plt', '.fai' or '.spn' extension"
			

		header = list(read_n_lines(fh,4))
		#print "header", header
		#Versions of Zgoubi up to 236 had just '...' in lines 3 and 4 of the fai file 
		test_version = header[2][0:3]
		old_format = False
		if test_version == '...':
			old_format = True
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
			
			if (id != None):
				try:
					if(lines[3].split()[1] != id):
						continue # escape quickly if this is not a point we are interested in
				except:
					print lines
					continue

			if (len(lines) == n_lines): # dont go through this if we have just few dangling lines on the end of the file
				# fill dictionary with the bits that have been identified
				#
				p = {}

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
						p['Y'] = self._bad_float(l1[1]) #current corrods
						p['T'] = self._bad_float(l1[2])

					#line 2
					if file != 'spn':
						l2 = scanf(lines[2], l2_re)
						p['Z'] = self._bad_float(l2[0])
						p['P'] = self._bad_float(l2[1])
						p['S'] = self._bad_float(l2[2])
						p['tof'] = self._bad_float(l2[3]) # time of flight
						if file=='plt': # a strange bug prehaps. anyways this makes it so that you get microseconds
							p['tof'] = p['tof'] / 1e5
						if file=='fai':
							p['KE'] = self._bad_float(l2[4]) # Kinetic energy
					
					#line 3
					# see 20080501 in lab book
					if (file=='plt'):
						l3 = scanf(lines[3], l3_re_plt)
						p['KART'] = int(l3[0])
						p['ID'] = l3[1]
						p['IREP'] = l3[2]
						p['SORT'] = l3[3]
						p['X'] = float(l3[4])
						p['By'] = float(l3[6])
						p['Bz'] = float(l3[7])
					elif (file=='fai'):
						l3 = scanf(lines[3], l3_re)
						p['ID'] = l3[0]
						p['IREP'] = l3[1]
						p['SORT'] = l3[2]
						p['X'] = -1 # X is not defined in FAI file, particle is always at end of the element

					#if  particle['X'] > 1e18: # this comes out as a silly giant number in some cases, maybe due to FAICNL or something
					#	particle['X'] = 0 

					#line4
					if (file=='plt'):
						l4 = scanf(lines[4], l4_re_plt)
						p['BORO'] = float(l4[3])
						p['PASS'] = int(l4[4])
						p['element_type'] = l4[5].strip()
						p['element_label1'] = l4[6].strip()
						p['element_label2'] = l4[7].strip()
					elif (file=='fai'):
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
					if file=='fai':
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
						p['element_label1']  = l0_stringpart[1].strip('\'')
						p['element_label2']  = l0_stringpart[2].strip('\'')
						p['LET'] = l0_stringpart[3].strip('\'')
					elif file=='plt':
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
						p['element_label1']  = l0_stringpart[1].strip()
						p['element_label2']  = l0_stringpart[2].strip()
						p['LET'] = l0_stringpart[3].strip()

				
				plt_data.append(p)
		return plt_data


	def get_track(self, file, coord_list, multi_list=None):
		"""
		returns a list of coordinates specified, and multiply them by the multiplier if not none
		eg get_trac('plt', ['X','Y','element_label'],[0.01,0.01,None])
		
		If all the columns requested are numerical, and new headered data formats are being used then this function will return a numpy array
		"""
		all = self.get_all(file)
		#check if we are using the new zgoubi.io version
		if (type(all) == type(numpy.zeros(0))):
			#coords = numpy.zeros([all.size, len(coord_list)])
			
			# is set logic to see what best array type is
			# mixed arrays might break some old code, so fall back to dicts
			dtypes =  set([all[c].dtype for c in coord_list])
			if not( dtypes - set([numpy.dtype('i')]) ):
				# if there is nothing but ints, use int
				best_type = 'i'
			elif not( dtypes - set([numpy.dtype('i'), numpy.dtype('f'), numpy.dtype('d')]) ):
				# if there is nothing but ints and float/double, use double
				best_type = 'd'
			else:
				# otherwise fall back to dicts
				best_type = "mixed"


			if best_type != "mixed":
				# all cols numeric, so give a fast numpy array
				coords = numpy.zeros([all.size, len(coord_list)], dtype=(best_type))

				for n,c in enumerate(coord_list):
					print c
					coords[:,n] = all[c]
					if multi_list:
						if multi_list[n]:
							coords[:,n] *= multi_list[n]
				return coords


		coords = []
		for p in all:
			this_coord = []
			for n,c in enumerate(coord_list):
				if (multi_list == None):
					this_coord.append(p[c])
				elif (multi_list[n] == None):
					this_coord.append(p[c])
				else:
					this_coord.append(p[c] * multi_list[n])
			coords.append(this_coord)
		return coords
	
	def get_bunch(self, file, end_label=None, old_bunch=None):
		""""Get back a bunch object from the fai file. It is recommended that you put a MARKER before the last FAISCNL, and pass its label as end_label, so that only the bunch at the final position will be returned. All but the final lap is ignored automatically.
		Optionally the an old_bunch can be passed to the function, its mass and charge will be copyed to the new bunch.
		"""
		all = self.get_all(file)

		if not (type(all) == type(numpy.zeros(0))):
			raise OldFormatError("get_bunch() only works with the new fai format")
		
		# select only the particles that made it to the last lap
		last_lap = all[ all['PASS'] == all['PASS'].max() ]
		# also select only particles at FAISTORE with matching end_label
		if end_label:
			end_label = end_label.ljust(8) # pad to match zgoubi
			last_lap = last_lap[ last_lap['element_label1'] == end_label ]

		#print last_lap[:10]['BORO']
		#print last_lap[:10]['D-1']

		# Create a new bunch and copy the values arcoss (with conversion to SI)
		last_bunch = zgoubi.bunch.Bunch(nparticles=last_lap.size, rigidity=last_lap[0]['BORO']/1000)
		if old_bunch != None:
			last_bunch.mass = old_bunch.mass
			last_bunch.charge = old_bunch.charge
		particles = last_bunch.particles()
#		print last_lap.dtype
		particles['Y'] = last_lap['Y'] /100
		particles['T'] = last_lap['T'] /1000
		particles['Z'] = last_lap['Z'] /100
		particles['P'] = last_lap['P'] /1000
		particles['D'] = last_lap['D-1'] +1

		return last_bunch


	def get_extremes(self, file, element_label=None, coord='Y', id=None):
		"""
		Get the max and min positions of a certain coordinate, in a certain element
		"""
		all = self.get_all(file, id=id)
		if (element_label == None):
			el_points = all
		else:
			el_points = [p for p in all if (p['element_label1'] == element_label or p['element_label2'] == element_label)]

		
		if (len(el_points) == 0):
			print "no points for particle id:", id, "element_label", element_label
			raise NoTrackError
			
		
		mini = min( [p[coord] for p in el_points] )
		maxi = max( [p[coord] for p in el_points] )

		return mini, maxi

	def list_particles(self, file):
		if(file == 'plt'):
			fh = self.line.plt_fh()
		elif(file == 'fai'):
			fh = self.line.fai_fh()
		else:
			raise ValueError, "get_all() expects name to be 'plt' or 'fai'"

		header = list(read_n_lines(fh,4))
		particle_ids = set()
		for lines in yield_n_lines(fh, 5):
			if (len(lines) == 5): # dont go through this if we have just few dangling lines on the end of the file
				#bits = [] # not computerscience bits. just the bits of data on each line
				#bits.append(line.split()) # split them by whitespace
				#particle_ids.add(int(bits[3][1]))
				particle_ids.add(int(lines[3].split()[1]))
		return particle_ids

	def in_bounds(self, file, element_label, min_bound, max_bound, coord='Y', verbose=False, id=None):
		"""
		check if particle exceeded bounds in this element. pass bounds in zgoubi's default unit for that coordinate
		"""
		if (id != None):
			part_list = [id]
		else:
			part_list =  self.list_particles(file)

		
		for id in part_list:
			try:
				track_min, track_max = self.get_extremes(file, element_label, coord, id=id)
			except NoTrackError:
				print "No Track"
				return True
			v=verbose # to save typing
			in_b = True
			if (track_min <= min_bound):
				if(v): print track_min, "<=", min_bound, "hit lower bound, in element", element_label
				in_b = False
				
			if (track_max >= max_bound):
				if(v): print track_max, ">=", max_bound, "hit upper bound, in element", element_label
				in_b = False

		return in_b

	def check_bounds(self, file, min_bounds, max_bounds, coord='Y', part_ids=None, assume_in_order=False):
		"""read whole file. for each chunk, check if it is an element with bounds, if yes check if it is inside
		record crashes, by particle id
		min_bounds and max_bounds should be of the form {'label': bound, 'label': bound}
		returns numpy arrays
		note: lost particles dont crash
		"""
		if (part_ids == None):
			particle_ids = self.list_particles(file)
		else:
			particle_ids = part_ids
		if(file == 'plt'):
			fh = self.line.plt_fh()
		elif(file == 'fai'):
			fh = self.line.fai_fh()
		else:
			raise ValueError, "get_all() expects name to be 'plt' or 'fai'"

		header = list(read_n_lines(fh,4))
		
		crashes = numpy.zeros(len(particle_ids), dtype=int)
		not_crashes = numpy.zeros(len(particle_ids), dtype=int)
		laps = numpy.zeros(len(particle_ids), dtype=int)
		crash_lap = numpy.zeros(len(particle_ids), dtype=int)
		crash_lap += 100000
		
		for lines in yield_n_lines(fh, 5):
			if (len(lines) != 5):
				continue
			bits = {} # not computerscience bits. just the bits of data on each line
			particle = {}
			#for line in lines:
			#	bits.append(line.split()) # split them by whitespace
			bits[1]= lines[1].split()
			bits[3]= lines[3].split()
			bits[4]= lines[4].split()
			
			particle['ID'] = bits[3][1]
			particle['Y'] = float(bits[1][1])
			particle['T'] = float(bits[1][2])
			
			if (assume_in_order):
				if (crashes[int(particle['ID'])-1]>0):
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
			if((len(bits[4]) - first_non_numeric) == 4): # 2 labels
				particle['element_label1'] = bits[4][first_non_numeric + 1]
				particle['element_label2'] = bits[4][first_non_numeric + 2]
			if((len(bits[4]) - first_non_numeric) == 3): # 1 labels
				particle['element_label1'] = bits[4][first_non_numeric + 1]
			
			has_crashed = False
			for k,v in min_bounds.items():
				if (particle['element_label1'] == k or particle['element_label2'] == k):
					if (particle['Y'] <= v):
						#crash
						has_crashed = True
			for k,v in max_bounds.items():
				if (particle['element_label1'] == k or particle['element_label2'] == k):
					if (particle['Y'] >= v):
						#crash
						has_crashed = True
			if (has_crashed):
				# zgoubi counts from 1, we count from zero
				crashes[int(particle['ID'])-1] +=1
				crash_lap[int(particle['ID'])-1] = min(crash_lap[int(particle['ID'])-1], int(particle['PASS'])) # record the first crash
			else:
				not_crashes[int(particle['ID'])-1] +=1

		return crashes, not_crashes, laps, crash_lap
	
	def get_tune(self):
		"""Returns a tuple (NU_Y, NU_Z) of the tunes.
		Needs a beam line is an OBJET type 5, and a MATRIX element.


		"""
		has_object5 = False
		has_matrix = False
		for e in self.line.elements():
			t = str(type(e)).split("'")[1].rpartition(".")[2]
			if t == 'OBJET5':
				has_object5 = True
			if t == 'MATRIX':
				if e.IORD == 1 and e.IFOC>10:
					has_matrix = True
				
		if not (has_object5 and has_matrix):
			raise BadLineError, "beamline need to have an OBJET with kobj=5 (OBJET5), and a MATRIX element with IORD=1 and IFOC>10 to get tune"

		found_matrix = False
		for line in self.line.res_fh():
			if "MATRIX" in line:
				bits = line.split()
				if bits[0].isdigit() and bits[1] == "MATRIX":
					found_matrix = True
			elif found_matrix and "NU" in line:
				bits = line.split()
				if (bits[0], bits[1], bits[3], bits[4]) == ("NU_Y", "=", "NU_Z", "="):
					try:
						NU_Y = float(bits[2])
					except ValueError:
						print "could not get Y tune from:"
						print line
						print "setting NU_Y to -1"
						NU_Y = -1 
					try:
						NU_Z = float(bits[5])
					except ValueError:
						print "could not make floats from line"
						print line
						print "setting NU_Z to -1"
						NU_Z = -1 
					print "Tune: ", (NU_Y, NU_Z)
					return (NU_Y, NU_Z)
		raise NoTrackError, "Could not find MATRIX output, maybe beam lost"

	def get_twiss_parameters(self):
		"""Returns a tuple (BETA_X,GAMMA_X,BETA_Y,GAMMA_Y) from the twiss parameter matrix.
		Needs a beam line is an OBJET type 5, and a MATRIX element.


		"""
		has_object5 = False
		has_matrix = False
		for e in self.line.elements():
			t = str(type(e)).split("'")[1].rpartition(".")[2]
			if t == 'OBJET5':
				has_object5 = True
			if t == 'MATRIX':
				has_matrix = True
		if not (has_object5 and has_matrix):
			raise BadLineError, "beamline need to have an OBJET with kobj=5 (OBJET5), and a MATRIX elementi to get tune"

		found_matrix = False
		found_twiss = False
		found_row1 = False
		found_row2 = False
		found_row3 = False
		found_row4 = False
		for line in self.line.res_fh():
			if "MATRIX" in line:
				bits = line.split()
				if bits[0].isdigit() and bits[1] == "MATRIX":
					found_matrix = True
			elif found_matrix and "beta/-alpha" in line:
				found_twiss = True
			elif found_twiss and not found_row1 and len(line) > 1:
				row = line.split()
				beta_y = float(row[0])
				alpha_y = -1*float(row[1])
				found_row1 = True
			elif found_row1 and not found_row2 and len(line) > 1:
				row = line.split()
				gamma_y = float(row[1])
				found_row2 = True
			elif found_row2 and not found_row3 and len(line) > 1:
				row = line.split()
				beta_z = float(row[2])
				alpha_z = -1*float(row[3])
				found_row3 = True
			elif found_row3 and not found_row4 and len(line) > 1:
				row = line.split()
				gamma_z = float(row[3])
				found_row4 = True
				return (beta_y, alpha_y, gamma_y, beta_z, alpha_z, gamma_z)

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
		has_reb = False
		for e in self.line.elements():
			t = str(type(e)).split("'")[1].rpartition(".")[2]
			if t == 'REBELOTE':
				has_reb = True
		if not (has_reb):
			raise BadLineError, "beamline need to have a REBELOTE for this function"

		for line in self.line.res_fh():
			if "End  of  'REBELOTE'  procedure" in line:
				print "REBELOTE completed"
				return True
		return False

	def run_success(self):
		"""Checks that zgoubi completed

		"""
		for line in self.res_fh():
			if "MAIN PROGRAM : Execution ended upon key  END" in line:
				return True
		return False
		
class Plotter(object):
	"""A plotter for beam lines and tracks.
	
	"""
	def __init__(self, line):
		"""Creates a new plot from the line.

		"""
		try:
			dummy = cairo.version
		except NameError:
			raise StandardError, "Cairo python bindings not found\nFor Debian/Ubuntu systems please install python-cairo"
		
		self.line = line
		
		self.sections = {}
		self.elements = {}
		#self.pipes = {}
		self.tracks = []
		
		self.line_length = 0

		self.width = 100 # physical size of canvas in cm
		self.height = 100 
		#self.last_added_elem = None
		self.pretty_names={}
		self.magnet_width = {}
		self._scan_line()
		#self.save_pdf('test.pdf')
	
	def _calc_physial_size(self):
		self.min_x = 1e6
		self.min_y = 1e6
		self.max_x = -1e6
		self.max_y = -1e6
		
		points = [] # list of points we want to fit on canvas
		
		for e in self.elements.values():
			points += e['box']

		for t in self.tracks:
			points += t['points']

		for x,y in points:
			self.min_x = min(self.min_x, x)
			self.max_x = max(self.max_x, x)
			self.min_y = min(self.min_y, y)
			self.max_y = max(self.max_y, y)
		
		#add a 5% margin
		self.min_x -= 0.05 * (self.max_x - self.min_x)
		self.max_x += 0.05 * (self.max_x - self.min_x)
		self.min_y -= 0.05 * (self.max_y - self.min_y) + 11# and a bit more for labels
		self.max_y += 0.05 * (self.max_y - self.min_y)
		
		self.aspect = (self.max_x - self.min_x)/(self.max_y - self.min_y)
		

	def _scan_line(self):
		"""Scan through the line, and make a note of where all the elements are.
		
		"""
		count = 0
		angle = 0
		position = [0,0]
		self.elements = {}
		for elem in self.line.elements():
			classtype =  str(type(elem))
			classtype = classtype.split("'")[-2]
			classtype = classtype.split(".")[-1]
			#print classtype

			if (classtype=='DRIFT'):
				section = {}
				section['label'] = elem.label1
				if section['label'] == "":
					print classtype, "in line with no label"
					section['label'] = str(count)
					count += 1
				section['len'] = elem.XL
				section['start'] = position
				section['ang'] = angle

				self.sections[section['label']] = section
				section['end'] = self.transform(section['label'], (section['len'], 0))		
				position = section['end']
				#print section
			elif(classtype=='CHANGREF'):
				x = elem.XCE
				y = elem.YCE
				position[0] += x * cos(angle) - y * sin(angle)
				position[1] += x * sin(angle) + y * cos(angle)
				angle += radians(elem.ALE)
				#print "changref", elem.XCE,elem.YCE,elem.ALE
				#print "new coord, angle", position, angle
			elif(classtype in ['MULTIPOL','QUADRUPO']):
				section = {}
				section['label'] = elem.label1
				if section['label'] == "":
					print classtype, "in line with no label"
					section['label'] = str(count)
					count += 1
				section['len'] = elem.XL
				section['start'] = position
				section['ang'] = angle

				self.sections[section['label']] = section
				section['end'] = self.transform(section['label'], (section['len'], 0))		
				position = section['end']
				#print section
				element = {}
				width = 5* elem.R_0
				if self.magnet_width.has_key(section['label']):
					width = self.magnet_width[section['label']]
				length = elem.XL
				a = self.transform(section['label'],(0, width/2))
				b = self.transform(section['label'],(length, width/2))
				c = self.transform(section['label'],(length, -width/2))
				d = self.transform(section['label'],(0, -width/2))
				element['box'] = (a,b,c,d)
				self.elements[section['label']] = element

	def transform(self, section, point):
		"""transform a point relative to a section into global coords.

		"""

		x0, y0 = self.sections[section]['start']
		ang = self.sections[section]['ang']
		
		x, y = point
		if x == -1: # adjustment for FAI files
			x = self.sections[section]['len']
		
		#if (cum_len):
		#	x -= self.sections[section]['prev_cum_len']
		#print ang, x0, y0	
		xp = x * cos(ang) - y * sin(ang) + x0
		yp = x * sin(ang) + y * cos(ang) + y0
		
		return [xp, yp]
		
	def add_results(self, results, colour=None):
		"""Add a results object to the plot. All tracks will be extracted.
		Multiple passes made with REBELOTE will be nicely draw ontop of each other. Multiple particles will make a mess
		Colour is optional, pass a tuple of (red, blue, green) with values from 0 off, to 1 on.

		"""
		raw_track = []
		try:
			raw_track += results.get_track('plt', ['X','Y','element_label1','PASS','tof'])
		except IOError, e:
			print "No PLT file for tracks:", e
		try:
			raw_track += results.get_track('fai', ['X','Y','element_label1','PASS','tof'])
		except IOError, e:
			print "No FAI file for tracks:", e
		if len(raw_track)==0:
			raise NoTrackError, "No tracks to plot"
		#for t in raw_track:
		#	print t
		raw_track.sort(key=itemgetter(4)) # sort points by time of flight
		passes = set([x[3] for x in raw_track])
		for thispass in passes: # deal with one lap at a time
			track = []
			for raw_p in [x for x in raw_track if x[3]==thispass]:
				p = self.transform(raw_p[2], raw_p[0:2]) # transform to element coordinates
				#print raw_p[2], self.sections[raw_p[2]], raw_p[0:2], p
				track.append(p)
			self.tracks.append(dict(points=track, colour=colour))


	def save_pdf(self, fname, canv_width=None, canv_height=None):
		"""Write plot to a PDF file.
		canv_width or canv_height can be specified in points (1/72 inch)
		"""
		self._calc_physial_size()
		if (canv_width != None):
			canv_height = canv_width/self.aspect
		elif (canv_height != None):
			canv_width = canv_height*self.aspect
		else:
			canv_width = 595
			canv_height = canv_width/self.aspect
			
		pdf_surf = cairo.PDFSurface(fname, canv_width, canv_height)
		cr = cairo.Context(pdf_surf)
		self.draw(cr, canv_width, canv_height, line_width=1)
	
	def save_png(self, fname, canv_width=None, canv_height=None):
		"""Write plot to a png file.
		canv_width or canv_height can be specified in pixels
		"""
		self._calc_physial_size()
		if (canv_width != None):
			canv_height = canv_width/self.aspect
		elif (canv_height != None):
			canv_width = canv_height*self.aspect
		else:
			canv_width = 800
			canv_height = canv_width/self.aspect
			
		img_surf = cairo.ImageSurface(cairo.FORMAT_RGB24,int(canv_width), int(canv_height))
		cr = cairo.Context(img_surf)
		self.draw(cr, canv_width, canv_height, line_width=1)
		img_surf.write_to_png(fname)
		
	def save_ps(self, fname, eps=False, canv_width=None, canv_height=None):
		"""Write plot to a PDF file.
		canv_width or canv_height can be specified in points (1/72 inch)
		"""
		self._calc_physial_size()
		if (canv_width != None):
			canv_height = canv_width/self.aspect
		elif (canv_height != None):
			canv_width = canv_height*self.aspect
		else:
			canv_width = 595
			canv_height = canv_width/self.aspect
			
		ps_surf = cairo.PSSurface(fname, canv_width, canv_height)
		if eps:  # needs a hot new version of cairo, so not implemented yet
			pass #ps_surf.PSSurface_set_eps(True)
		cr = cairo.Context(ps_surf)
		self.draw(cr, canv_width, canv_height, line_width=1)

	def set_pretty_names(self, pretty_names):
		self.pretty_names = pretty_names
	
	def set_magnet_width(self, magwidth):
		self.magnet_width = magwidth
		self._scan_line() # need to go through the line again with the new widths
		
	def draw(self, cr, canv_width, canv_height, line_width=1):
		"""Draw plot to a cairo surface.

		"""
		cr.set_font_size(18)
		#self._calc_physial_size()
		#self.width = self.max_x - self.min_x
		self.width = (self.max_x - self.min_x)
		self.height = (self.max_y - self.min_y)
		scale_fac = min(canv_width/self.width, canv_height/self.height)

		cr.scale(scale_fac, -scale_fac)

		#cr.translate(0,-self.height/2)	
		cr.translate(-self.min_x,-self.max_y)	

		cr.set_source_rgb(1,1,1)
		cr.rectangle(0, -self.height/2, self.width, self.height)
		cr.fill()
		
		
		cr.set_line_width(max(cr.device_to_user_distance(line_width, line_width)))


		cr.set_source_rgb(0.5,0.5,0.5)
		#draw line sections
		for name, l in self.sections.items():
			cr.move_to(*l['start'])
			cr.line_to(*l['end'])
			
		cr.stroke()
		
		cr.set_source_rgb(0.5,0.5,0.5)
		# draw_elements
		for name, e in self.elements.items():
			cr.move_to(*e['box'][0])
			cr.line_to(*e['box'][1])
			cr.line_to(*e['box'][2])
			cr.line_to(*e['box'][3])
			cr.line_to(*e['box'][0])
			cr.stroke()
			box_min_y = box_min_x = 10000
			box_max_x = -10000
			for p in e['box']:
				box_min_y = min(box_min_y, p[1])
				box_min_x = min(box_min_x, p[0])
				box_max_x = max(box_max_x, p[0])
			
			pname = self.pretty_names.get(name, name)
			#this is ugly. need to draw in name coordinate space, because the real space coords are upside down
			x_bearing, y_bearing, extent_width, extent_height = cr.text_extents(pname)[:4]
			cr.move_to((box_min_x+box_max_x)/2-extent_width/2/scale_fac, box_min_y-extent_height/scale_fac)
			cr.identity_matrix()
			cr.show_text(pname)
			cr.scale(scale_fac, -scale_fac)
			cr.translate(-self.min_x,-self.max_y)	

		#draw tracks
		cr.set_source_rgb(0,0,0)
		for track in self.tracks:
			if track['colour'] != None:
				cr.set_source_rgb(*track['colour'])
			cr.move_to(*track['points'][0])
			for p in track['points']:
				cr.line_to(*p)
			cr.stroke()

		
		#draw scale
		scale_start = 0
		scale_end = 100
		x_bearing, y_bearing, extent_width, extent_height = cr.text_extents('1')[:4]
		scale_pos = self.min_y + 2.1 + extent_height*2/scale_fac
		cr.set_source_rgb(0,0,0)
		cr.move_to(scale_start, scale_pos)
		cr.line_to(scale_end, scale_pos)
		cr.stroke()
		for x in xrange(11):
			x_pos = (scale_end - scale_start)/10*x + scale_start
			cr.move_to(x_pos, scale_pos)
			cr.line_to(x_pos,scale_pos-1)
			cr.stroke()
			if (x%5 ==0):
				x_bearing, y_bearing, extent_width, extent_height = cr.text_extents(str((scale_end - scale_start)/10*x/100))[:4]
				cr.move_to(x_pos-extent_width/2/scale_fac, scale_pos-extent_height/scale_fac-1.5)
				cr.identity_matrix()
				cr.show_text(str((scale_end - scale_start)/10*x/100))
				cr.scale(scale_fac, -scale_fac)
				cr.translate(-self.min_x,-self.max_y)

			if (x ==0):
				x_bearing, y_bearing, extent_width, extent_height = cr.text_extents(str((scale_end - scale_start)/10*x))[:4]
				cr.move_to(x_pos-extent_width/2/scale_fac, scale_pos-extent_height/scale_fac*2-2)
				cr.identity_matrix()
				cr.show_text("m")
				cr.scale(scale_fac, -scale_fac)
				cr.translate(-self.min_x,-self.max_y)





