#!/usr/bin/env python
import os
import shutil
import subprocess
from zgoubi import common

class ZgoubiBuildError(Exception):
	"Problem building Zgoubi"
	pass


zgoubi_build_dir = os.path.expanduser("~/.pyzgoubi/build")
zgoubi_build_dir2 = os.path.join(zgoubi_build_dir, "zgoubi-trunk")
zgoubi_install_dir = os.path.expanduser("~/.pyzgoubi/bin")
zgoubi_svn_address = "https://zgoubi.svn.sourceforge.net/svnroot/zgoubi/trunk"

def get_zgoubi_svn():
	"Download zgoubi from SVN"
	devnull = open("/dev/null", "w")
	try:
		ret = subprocess.call(['svn', '--version'], stdout=devnull)
	except OSError:
		raise ZgoubiBuildError("svn not found: install subversion")
	if os.path.isdir(zgoubi_build_dir):
		print "Zgoubi build folder already exists:", zgoubi_build_dir
		ret = subprocess.call(['svn', 'info'], cwd=zgoubi_build_dir2)
		if ret == 0:
			print "Zgoubi SVN already present at:", zgoubi_build_dir
			print "Reverting local changes"
			ret2 = subprocess.call(['svn', 'revert', '-R', '.'], cwd=zgoubi_build_dir2)
			return
		
		if ret != 0 or ret2 != 0:
			raise ZgoubiBuildError("%s already exists, but does not contain a working checkout")

	print "Downloading Zgoubi SVN:", zgoubi_svn_address
	common.mkdir_p(zgoubi_build_dir)
	ret = subprocess.call(['svn', 'co', zgoubi_svn_address, "zgoubi-trunk"], cwd=zgoubi_build_dir)
	if ret != 0:
		raise ZgoubiBuildError("SVN download failed")


def set_zgoubi_version(version=None):
	"Set downloaded SVN to a given version. or, if no version given, to latest version"
	ret = subprocess.call(['svn', 'revert','-R',  '.'], cwd=zgoubi_build_dir2)
	if version == None:
		ret = subprocess.call(['svn', 'update'], cwd=zgoubi_build_dir2)
	else:
		ret = subprocess.call(['svn', 'update', '-r', '%s'%version], cwd=zgoubi_build_dir2)


def apply_zgoubi_patches(patches):
	"Download and apply a set of patches"
	for patch in patches:
		print "applying", patch
		patchname = patch.rpartition("/")[2]
		ret = subprocess.call(['wget', patch, "-O", patchname], cwd=zgoubi_build_dir2)
		if ret != 0:
			raise ZgoubiBuildError("Patch download failed: %s" % patch)
		ret = subprocess.call('patch -p0 < %s' % patchname, cwd=zgoubi_build_dir2, shell=True)
		if ret != 0:
			raise ZgoubiBuildError("Patch application failed: %s" % patch)


def make_zgoubi(threads=2):
	"Build zgoubi source code"
	print "building zgoubi"
	ret = subprocess.call(['make', 'clean' ], cwd=zgoubi_build_dir2)
	if ret != 0:
		raise ZgoubiBuildError("Make clean failed")
	ret = subprocess.call(['make', '-j%d'%threads ], cwd=zgoubi_build_dir2)
	if ret != 0:
		raise ZgoubiBuildError("Building zgoubi failed")

def install_zgoubi(postfix=""):
	"Install zgoubi into ~/.pyzgoubi folder"
	common.mkdir_p(zgoubi_install_dir)

	shutil.copy(os.path.join(zgoubi_build_dir2, "zgoubi", "zgoubi"),
	                os.path.join(zgoubi_install_dir, "zgoubi%s"%postfix )  )
	shutil.copy(os.path.join(zgoubi_build_dir2, "zpop", "zpop"),
	                os.path.join(zgoubi_install_dir, "zpop%s"%postfix )  )


def install_zgoubi_all():
	"This currently install a version of zgoubi known to work with pyzgoubi"
	get_zgoubi_svn()
	set_zgoubi_version(261)
	patches=[
	"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/mcmodel-fix.diff",
	"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/zgoubi_parallel_build.diff",
	"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/objet3.diff",
	"http://www.hep.man.ac.uk/u/sam/pyzgoubi/zgoubipatches/kobj301.diff",
	]

	apply_zgoubi_patches(patches)
	make_zgoubi()

	install_zgoubi("_261+patches")
	print "\nInstalled zgoubi into ", zgoubi_install_dir
	print "Add the following line to ~/.pyzgoubi/settings.ini"
	print "zgoubi_path = %s/zgoubi%s" % (zgoubi_install_dir, "_261+patches")


if __name__ == "__main__":
	install_zgoubi_all()
