#!/usr/bin/env python
"""Common functions that can be used anywhere in pyzgoubi

"""

import errno
import os

def mkdir_p(path):
	"""Make a directory. Make parents if needed. Don't raise error if it already exists
	From: http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python code by tzotzioy http://stackoverflow.com/users/6899/
	"""
	try:
		os.makedirs(path)
	except OSError as exc:
		if exc.errno == errno.EEXIST:
			pass
		else:
			raise

def open_file_or_name(forn, mode="r", mkdir=False):
	"""Pass either a filename or file handle like object. Returns a file like object
	If mode is 'w' or 'a' then can set mkdir=True to make sure the directory is created if needed.
	"""
	if hasattr(forn, 'readline'):
		return forn
	else:
		if mkdir and mode in ["w", "a"]:
			mkdir_p(os.path.dirname(forn))
		return open(forn, mode)


