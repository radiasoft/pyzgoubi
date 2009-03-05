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
from zgoubi_utils import *
from zgoubi import *
from zgoubi_constants import *

# create a set of functions that work on default_line and default_results objects
# these will make it simpler to run simple simulations

def make_line(*a, **k):
	global default_line
	default_line = Line(*a, **k)

def run(*a, **k):
	global default_line, default_result
	default_result = default_line.run(*a, **k)

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

res_funcs += ['get_all', 'get_track', 'get_tune', 'get_twiss_parameters']

for func_name in res_funcs:
	code = function_template%dict(func_name=func_name, class_name='result')
	exec(code)

if __name__ == '__main__':
	try:
		input_file_name = sys.argv[1]
	except IndexError:
		print "No input file, try:"
		print "pyzgoubi inputfile"
		sys.exit()
	
	if not os.path.exists(input_file_name):
		print "no such input file in current directory"
		sys.exit()
		
	execfile(input_file_name)
	

	for left_over in locals().values():
		if "Line"  in str(type(left_over)).split("'")[1]:
			del(left_over)




