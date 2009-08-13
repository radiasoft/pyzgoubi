#!/usr/bin/env python
from __future__ import division

import ConfigParser
import os

home = os.path.expanduser('~')
config_dir = os.path.join(home, ".pyzgoubi")

if not os.path.exists(config_dir):
	os.mkdir(config_dir)


config_path = os.path.join(config_dir, "settings.ini")

config = ConfigParser.RawConfigParser()
config.add_section('pyzgoubi')
#default values
config.set('pyzgoubi', 'tmp_dir', '/tmp/')
config.set('pyzgoubi', 'extra_defs_files', '')
#config.set('pyzgoubi', 'pyzgoubi_dir', "/home/sam/pyzgoubi/")
#config.set('pyzgoubi', 'defs_dir', "/home/sam/pyzgoubi/defs")
config.set('pyzgoubi', 'zgoubi_path', "zgoubi")


if os.path.exists(config_path):
	config.read(config_path)
else:
	config.write(open(config_path, 'w'))

# check zgoubi binary
import subprocess
try:
	a = subprocess.Popen([config.get('pyzgoubi', 'zgoubi_path')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except OSError:
	print "can't zgoubi binary", config.get('pyzgoubi', 'zgoubi_path')
	print "please modify the 'zgoubi_path' entry in", config_path, "to the full path of zgoubi"
	print "or add the direcory containing zgoubi to the $PATH"
	import sys
	sys.exit()


zgoubi_settings = {}
zgoubi_settings['tmp_dir'] = config.get('pyzgoubi', 'tmp_dir')
#zgoubi_settings['pyzgoubi_dir'] = config.get('pyzgoubi', 'pyzgoubi_dir')
#zgoubi_settings['defs_dir'] = config.get('pyzgoubi', 'defs_dir')

zgoubi_settings['zgoubi_path'] = config.get('pyzgoubi', 'zgoubi_path')

zgoubi_settings['extra_defs_files'] = config.get('pyzgoubi', 'extra_defs_files').split(',')
if zgoubi_settings['extra_defs_files'] == ['']:
	zgoubi_settings['extra_defs_files'] = []


# create example defs file
example_defs_path = os.path.join(config_dir, "user_elements.defs")

if not os.path.exists(example_defs_path):
	example_defs = open(example_defs_path, 'w')
	example_defs.write("# define any elements you need here, the will be compiled to ~/.pyzgoubi/defs_cache\n")

