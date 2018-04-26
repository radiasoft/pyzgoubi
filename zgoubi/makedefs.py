#!/usr/bin/env python

"Parser for definitions formats"

from __future__ import division
from string import ascii_letters, digits


def expand_types(typestring):
	"""convert 'I,3E,2I' to ['I','E','E','E','I','I']
	This is used internally to parse the definitions file
	
	"""
	types = []
	chunks = [x.strip() for x in typestring.split(',')]
	for chunk in chunks:
		multi = "" # how many
		#print chunk
		ptype = ""  # int, float, strin
		count = "" # sting len
		chunk = list(chunk)
		# pop the digits at the front into multi
		while (chunk[0] in digits):
			multi += chunk.pop(0)
		if multi == "":
			multi = '1'
		ptype = chunk.pop(0) # first non digit is type
		#print multi, ptype  
		if (ptype in ['I', 'E', 'X']):
			for x in range(int(multi)):
				types.append(ptype)
			continue
		#print chunk		
		if (ptype == 'A'):
			count = int(''.join(chunk)) # what is left must be count
			for x in range(int(multi)):
				types.append(ptype+str(count))
			continue
		print multi, ptype, count
		print "can't parse", typestring
		raise ValueError
	return types


def make_element_classes(definitions_paths, compiled_path):
	"This reads the defs files, picks out an element and passes it to make_element_class"

	#read all the files into def_lines
	try:
		out_file = open(compiled_path, 'w')
	except IOError:
		print "Cant create file: " + compiled_path + "\n"
		print "Check that you have permission to modify files in that directory\n"
		raise
	
	out_file.write("# Generated definitions file\n")
	out_file.write("from zgoubi.core import zgoubi_element\n\n")

		
	def_lines = []
	for file_path in definitions_paths:
		try:
			fh = open(file_path)
		except IOError:
			print "can't read", file_path
			continue

		def_lines += fh.readlines()
		def_lines += '\n\n' # elements are separated by new lines

	defs = []
	for line in def_lines:
		pline = line.split('#')[0] # ignore anything after a '#' (comments)
		pline = pline.strip()

		if (line.strip() == "" and len(defs) > 0):
			# if we read a blank line, and we have gathered some lines
			# then we have a whole element
		   
			# create python class code from definition
			out_file.write(make_element_class(defs))
			out_file.write('\n')

			defs = []
			continue
		elif(len(pline) > 0):
			#otherwise gather the line into defs
			defs.append(pline)


def make_element_class(defs):
	"Convert the definition into a python class"
	# class name is the name used for the python class
	# zname is the zgoubi element name 
	class_name = defs[0]

	# must be a valid identifier name
	for letter in class_name:
		if letter not in ascii_letters + digits:
			raise ValueError, "invalid classname in definitions file: " + class_name

	if class_name[0] in digits:
		raise ValueError, "invalid classname in definitions file: " + class_name
		
	zname = defs[1]
	class_code = ""
	init_params_code = ""
	output_code = ""
	loop_output_code = ""
	add_func_code = ""
	
#	print zname
	
	class_code += "class %s(zgoubi_element):\n" % class_name
	class_code += "\t_class_name='%s'\n" % class_name
	class_code += "\t_zgoubi_name='%s'\n" % class_name
	class_code += "\tdef __init__(self, label1='', label2='', **settings):\n"
	class_code += "\t\tself._params={}\n"
	class_code += "\t\tself._types={}\n"
	class_code += "\t\tself._params['label1'] = label1\n"
	class_code += "\t\tself._params['label2'] = label2\n"
	
	output_code += "\tdef output(self):\n"
	output_code += "\t\tI=self.i2s\n"
	output_code += "\t\tE=self.f2s\n"
	output_code += "\t\tA=str\n"
	output_code += "\t\tL=self.l2s\n"
	output_code += "\t\tX=self.x2s\n"

	output_code += "\t\tout = ''\n"
	output_code += "\t\tnl = '\\n'\n"
	output_code += "\t\tsq = '\\''\n"
	output_code += "\t\tout += sq + '%s' + sq + ' ' + L(self._params['label1']) + ' ' + L(self._params['label2'])  + nl \n" % (zname)
	
	add_func_code += "\tdef add(self, **settings):\n"
	add_func_code += "\t\tparams = {}\n"
	
	#get all params  and types
	cond_code = ""
	in_loop = has_loops = False
	for l in defs[2:]:
		# conditional output
		#print '"', l, '"'
		if l.startswith('!'):
			conditions = l[1:]
			# parameter equals x type conditions
			if '==' in conditions:
				cond_param, cond_value = conditions.split('==')
				cond_code = "if self._params['%s'] == %s : " % (cond_param, cond_value)
			# blank condition, clears condition
			elif l.strip() == "!":
				cond_code = ""
			# looping for elements such as FFAG, OBJET2
			elif '*{' in conditions:
				loop_on = conditions.strip('!*{') # the param that holds number of parts
				in_loop = has_loops = True
		#		print "start loop, ", loop_on
			elif l.strip() == "!}":
				in_loop = False
				output_code += "\t\tfor part in self._looped_data:\n"
				output_code += loop_output_code
				loop_output_code = ''
			continue
		# lines that define output look like
		# NAME1, NAME2 : TYPE1, TYPE2
		# maybe one day might look like
		# NAME1, NAME2 : TYPE1, TYPE2 : new extra stuff
		namess, typess = l.split(':')[0:2]
		names = [x.strip() for x in namess.split(',')] # split out elements
		types = expand_types(typess.strip())
		if(len(names) != len(types)):
			error = "Error in defs file:\n"
			error = "In element %s:\n" % class_name
			error += l + '\n'
			error += str(len(names)) + " names, but only"+ str(len(types))+ " types\n"
			raise ValueError, error
		
		for name, ptype in zip(names, types):
			if not in_loop:
				if ptype in ['I', 'E', 'X']:
					# let these types default to 0
					init_params_code += "\t\tself._params['%s']=0\n" % name
				elif ptype.startswith('A'):
					# let these types default to an empty string
					init_params_code += "\t\tself._params['%s']=''\n" % name
				else:
					print "unknown type", ptype, "for", name, "in", class_name
					raise ValueError
			else:
				if ptype in ['I', 'E', 'X']:
					# let these types default to 0
					add_func_code += "\t\tparams['%s']=0\n" % name
				elif ptype.startswith('A'):
					# let these types default to an empty string
					add_func_code += "\t\tparams['%s']=''\n" % name
				else:
					print "unknown type", ptype, "for", name, "in", class_name
					raise ValueError

				
			init_params_code += "\t\tself._types['%s']='%s'\n" % (name, ptype)
		
		#the lines that out put the parameters to zgoubi.dat
		if not in_loop:
			out_bits = []
			for name, ptype in zip(names, types):
				out_bits.append("%s(self._params['%s'])" % (ptype[0], name)) # type[0] so that A80 calls the A function
			output_code += "\t\t%sout +=  %s +nl \n" % (cond_code, " + ' ' + ".join(out_bits))
		else:
			out_bits = []
			for name, ptype in zip(names, types):
				out_bits.append("%s(part['%s'])" % (ptype[0], name))
			loop_output_code += "\t\t\t%sout +=  %s +nl \n" % (cond_code, " + ' ' + ".join(out_bits))
		
		
	output_code += "\t\treturn out\n"

	init_params_code += "\t\tself.set(settings)\n"
	if has_loops:
		init_params_code += "\t\tself._looped_data = []\n"
	
		add_func_code += "\t\tfor k, v in settings.items():\n"
		add_func_code += "\t\t\tif not params.has_key(k):\n"
		add_func_code += "\t\t\t\traise ValueError('Sub element of %s does not have parameter %s'%(self._class_name, k))\n"
		add_func_code += "\t\t\tparams[k] = v\n"
		add_func_code += "\t\tself._looped_data.append(params)\n"
		add_func_code += "\t\tself._params['%s'] = len(self._looped_data)\n" % loop_on
		add_func_code += "\tdef clear(self):\n"
		add_func_code += "\t\tself._looped_data = []\n"
	

	code = ''
	code += class_code
	code += init_params_code
	code += output_code
	if has_loops:
		code += add_func_code
	return code

