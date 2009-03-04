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


if __name__ == '__main__':
	from zgoubi_utils import *
	from zgoubi import *
	from zgoubi_constants import *
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




