#!/usr/bin/env python

#
# run './install.py' to start installer
#
import sys, os
try:
	import readline
except ImportError:
	# Python on Mac has no readline
	pass

from glob import glob

try:
	"abc".rpartition('b')
	from operator import itemgetter
	print "Python 2.5 or greater found"
except AttributeError:
	print "pyzgoubi needs Python version 2.5 or greater"
	sys.exit()


try:
	import numpy
	print "Numpy found"
except ImportError:
	print "pyzgoubi needs numpy"
	sys.exit()

try:
	install_mode = sys.argv[1]
except IndexError:
	install_mode = 'install'

things_to_copy = ['zgoubi.py',
				'zgoubi_makedefs.py',
				'zgoubi_utils.py',
				'zgoubi_constants.py',
				'pyzgoubi.py',
				'defs/simple_elements.defs',
				'defs/static_defs.py'
				]
things_to_copy += glob('examples/*')

dirs_to_make = ['defs', 'examples']

def install_files(prefix, zgoubi_path, dirs_to_make, things_to_copy):
	if not os.path.isdir(prefix):
		os.mkdir(prefix)

	for dirname in dirs_to_make:
		new_dir = os.path.join(prefix, dirname)
		if not os.path.isdir(new_dir):
			os.mkdir(new_dir)

	import shutil
	for thing in things_to_copy:
		print "Installing", thing
		dest = os.path.join(prefix, thing)
		shutil.copyfile(thing, dest)
	
	settings_dict = dict(pyzgoubi_path=prefix, zgoubi_path=zgoubi_path)
	settings_text = open('zgoubi_settings.py.template').read()
	settings_file = open(os.path.join(prefix, 'zgoubi_settings.py'), 'w')
	settings_file.write(settings_text%settings_dict)


	bash_env_text = open('pyzgoubi_env.bash').read()
	bash_env_file = open(os.path.join(prefix, 'pyzgoubi_env.bash'), 'w')
	bash_env_file.write(bash_env_text%settings_dict)

if install_mode == 'install':
	print "Where do you want to install? (eg /opt/pyzgoubi, ~/pyzgoubi)"
	prefix = raw_input().strip()
	if not prefix.endswith(os.sep):
		prefix += os.sep

	prefix = os.path.expanduser(prefix)

	print "Install to:", prefix, "(y/n)"
	if raw_input().strip() != 'y':
		print "aborting install"
		sys.exit()

	
	zgoubi_path= raw_input('please enter the path to the zgoubi binary:\n')
	if not os.path.exists(zgoubi_path):
		print "no such file:", zgoubi_path
		print "exiting"
		sys.exit()
	
	install_files(prefix, zgoubi_path, dirs_to_make, things_to_copy)
	
	env_command = "source %spyzgoubi_env.bash"%prefix

	print "To set up your environment to run pyzgoubi you need to run:"
	print env_command
	print "If you add that line to your .bashrc file this will bedone automaticall whenever you start a terminal session"

	print "Do you want this done automatically (y/n)"

	if (raw_input().strip()=='y'):
		bashrc = open(os.path.expanduser('~/.bashrc'), 'a')
		bashrc.write('\n')
		bashrc.write('#set up pyzgoubi\n')
		bashrc.write(env_command + '\n')
		bashrc.close()
		print ".bashrc modified"
		
elif install_mode == 'update':

	#find old path
	import popen2
#	out, err = popen2.popen2('bash -lc "alias pyzgoubi"')
#	try:
#		pyzgoubi_alias =  out.read().strip().split("'")[1].rpartition('=')[2].strip()
#		old_path = pyzgoubi_alias.split(' ')[0]
#	except IndexError:
#		print "pyzgoubi appears not to be installed"
#		print "exiting"
		
	out, err = popen2.popen2('grep pyzgoubi %s'%os.path.join(os.path.expanduser('~'), '.bashrc'))
	old_path =  out.readlines()[-1].split(' ')[1].strip()
	if (old_path.endswith('/pyzgoubi_env.bash')):
		old_path = old_path.rpartition('/')[0]
	else:
		print "pyzgoubi appears not to be installed"
		print "exiting"
		sys.exit() 
	print "installing to:", old_path
		
	
	sys.path.insert(0,old_path)
	import zgoubi_settings
	try:
		prefix = zgoubi_settings.pyzgoubi_dir
	except AttributeError:
		print "pyzgoubi appears not to be installed"
		print "exiting"
		sys.exit()
	zgoubi_path = zgoubi_settings.zgoubi_path

	install_files(prefix, zgoubi_path, dirs_to_make, things_to_copy)

else:
	print "unknown option:", install_mode
	print "try:"
	print "./install"
	print "or"
	print "./install update"




print "done"
print "you may need to start a new terminal, or run 'source ~/.bashrc' before using pyzgoubi"

