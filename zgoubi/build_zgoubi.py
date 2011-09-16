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
	common.mkdir_p(zgoubi_build_dir)
	try:
		ret = subprocess.call(['svn', '--version'])
	except OSError:
		raise ZgoubiBuildError("svn not found: install subversion")
	ret = subprocess.call(['svn', 'co', zgoubi_svn_address, "zgoubi-trunk"], cwd=zgoubi_build_dir)
	if ret != 0:
		raise ZgoubiBuildError("svn download failed")


def set_zgoubi_version(version=None):
	ret = subprocess.call(['svn', 'revert','-R',  '.'], cwd=zgoubi_build_dir2)
	if version == None:
		ret = subprocess.call(['svn', 'update'], cwd=zgoubi_build_dir2)
	else:
		ret = subprocess.call(['svn', 'update', '-r', '%s'%version], cwd=zgoubi_build_dir2)


def apply_zgoubi_patches(patches):
	for patch in patches:
		print "applying", patch
		patchname = patch.rpartition("/")[2]
		ret = subprocess.call(['wget', patch, "-O", patchname], cwd=zgoubi_build_dir2)
		ret = subprocess.call('patch -p0 < %s' % patchname, cwd=zgoubi_build_dir2, shell=True)


def make_zgoubi(threads=2):
	print "building zgoubi"
	ret = subprocess.call(['make', 'clean' ], cwd=zgoubi_build_dir2)
	ret = subprocess.call(['make', '-j%d'%threads ], cwd=zgoubi_build_dir2)

def install_zgoubi(postfix=""):
	common.mkdir_p(zgoubi_install_dir)

	shutil.copy(os.path.join(zgoubi_build_dir2, "zgoubi", "zgoubi"),
	                os.path.join(zgoubi_install_dir, "zgoubi_%s"%postfix )  )
	shutil.copy(os.path.join(zgoubi_build_dir2, "zpop", "zpop"),
	                os.path.join(zgoubi_install_dir, "zpop_%s"%postfix )  )



if __name__ == "__main__":
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

	install_zgoubi("261+patches")


