#!/usr/bin/env python

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
from math import *
from zgoubi.utils import *
from zgoubi.core import *
from zgoubi.constants import *

zgoubi_version = "0.3dev"

# create a set of functions that work on default_line and default_results objects
# these will make it simpler to run simple simulations

def make_line(*a, **k):
	global default_line
	default_line = Line(*a, **k)

def run(*a, **k):
	global default_line, default_result
	default_result = default_line.run(*a, **k)
	return default_result

# most of these functions can be built from a template

function_template="""def %(func_name)s(*a, **k):
	global default_%(class_name)s
	return default_%(class_name)s.%(func_name)s(*a, **k)
"""
# templatable functions for the line class
line_funcs = ['add', 'output', 'full_tracking', 'remove_looping', 'add_input_files', 'replace']

for func_name in line_funcs:
	code = function_template%dict(func_name=func_name, class_name='line')
	exec(code)

# templatable functions for the res class

# this pattern builds all the file access functions by pattern
res_funcs = [a%b for a in ['%s_fh', '%s', 'save_%s'] for b in ['res', 'fai', 'dat', 'plt', 'b_fai']]

res_funcs += ['get_all', 'get_track', 'get_tune', 'get_twiss_parameters', 'run_success']

for func_name in res_funcs:
	code = function_template%dict(func_name=func_name, class_name='result')
	exec(code)

def _show_help():
	if len(sys.argv) == 2:
		print "For docementation see http://www.hep.manchester.ac.uk/u/sam/pyzgoubi/"
		sys.exit(0)
	if sys.argv[2].lower() == "elements":
		print "avaliable elements"
		elements = []

		for k,v in globals().items():
			try:
				for b in v.__bases__:
					if "zgoubi_element" in str(b):
						elements.append(k)
			except AttributeError:
				pass
		elements.sort()
		print '\n'.join(elements)
		sys.exit(0)

	try:
		help_elem = globals()[sys.argv[2].upper()]
	except KeyError:
		print "There is no help for %s"%sys.argv[2]
		sys.exit(1)

	print sys.argv[2].upper()
	help_elem_inst = eval("%s()"%sys.argv[2].upper())
	print "zgoubi name:",help_elem_inst._zgoubi_name
	print "Parameters:"
	print '\n'.join(help_elem_inst.list_params())


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print "Usage:"
		print sys.argv[0], "input_file_name"
		print sys.argv[0], "--help"
		sys.exit(1)
	
	if sys.argv[1] in ['version', '--version']:
		print "Pyzgoubi version:%s"%zgoubi_version
		sys.exit(0)
	if sys.argv[1] in ['help', '--help']:
		_show_help()
		sys.exit(0)
	try:
		input_file_name = sys.argv[1]
	except IndexError:
		print "No input file, try:"
		print "pyzgoubi inputfile"
		sys.exit(1)
	
	if not os.path.exists(input_file_name):
		print "no such input file in current directory"
		sys.exit(1)
		
	execfile(input_file_name)
	

	for left_over in locals().values():
		if "Line"  in str(type(left_over)).split("'")[1]:
			del(left_over)




