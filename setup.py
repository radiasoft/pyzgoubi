#!/usr/bin/env python
from distutils.core import setup
from distutils.sysconfig import get_python_lib
import sys
import os
from glob import glob


setup(name='pyzgoubi',
	version='0.3',
	packages=['zgoubi'],
	scripts=['pyzgoubi'],
	#package_data={'zgoubi': ['defs/*.py', 'defs/*.defs']},
	data_files=[('share/pyzgoubi/definitions',glob('defs/*')),
	            ('share/pyzgoubi/examples',glob('examples/*')),
	            ('share/pyzgoubi/test',glob('tests/*.py')),
	            ('share/pyzgoubi/doc', glob('doc/*'))
	            ]
	)


if ("install" in sys.argv) and not ( "--help" in sys.argv):
	#find the log
	try:
		logfile = [x for x in sys.argv if x.startswith('--record')][-1]
		#strip away quotes, and expand path
		logfile = logfile.split('=')[-1].strip('"\'')
		logfile = os.path.expanduser(logfile)
	except IndexError:
		logfile = "install.log"

	try:
		for line in open(logfile):
			line = line.strip()
			if line.endswith("bin/pyzgoubi"):
				bin_path = os.path.dirname(line)
			if line.endswith("zgoubi/__init__.py"):
				lib_path = os.path.normpath(os.path.join(os.path.dirname(line),".."))
	except IOError:
		#error below will show
		pass
		
	try:
		print "\nyou may need to add the following to your .bashrc"
		print "export PYTHONPATH=$PYTHONPATH:%s"%lib_path
		print "export PATH=$PATH:%s"%bin_path
		print "or"
		print 'alias pyzgoubi="PYTHONPATH=%s python %s/pyzgoubi"'%(lib_path, bin_path)
	except NameError:
		print "Could not see install log, install may have failed."
		print "Can't give help with setting up path"
