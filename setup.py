#!/usr/bin/env python
from distutils.core import setup
from distutils.sysconfig import get_python_lib
import sys
import os
from glob import glob


setup(name='pyzgoubi',
	version='0.3',
	packages=['zgoubi'],
	scripts=['pyzgoubi.py'],
	#package_data={'zgoubi': ['defs/*.py', 'defs/*.defs']},
	data_files=[('share/pyzgoubi/definitions',glob('defs/*')),
	            ('share/pyzgoubi/examples',glob('examples/*')),
	            ('share/pyzgoubi/test',glob('tests/*')),
	            ('share/pyzgoubi/doc', glob('doc/*'))
	            ]
	)


if "install" in sys.argv:
	#find the prefix argument
	prefix = [x for x in sys.argv if x.startswith('--prefix')][-1]
	#strip away quotes, and expand path
	prefix = prefix.split('=')[-1].strip('"\'')
	prefix = os.path.expanduser(prefix)
	prefix = os.path.join(os.getcwd(), prefix)
	print "\nyou may need to add the following to your .bashrc"
	print "export PYTHONPATH=$PYTHONPATH:%s"%get_python_lib(True,False, prefix)
	print "export PATH=$PATH:%s/bin"%prefix
	print "or"
	print 'alias pyzgoubi="PYTHONPATH=%s python %s/bin/pyzgoubi.py"'%(get_python_lib(True,False, prefix), prefix)


