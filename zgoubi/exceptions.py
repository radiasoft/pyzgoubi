#!/usr/bin/env python

"Exceptions raised by PyZgoubi"

class NoTrackError(Exception):
	"No particle track"
	pass

class BadLineError(Exception):
	"Line not sufficient for this function"
	pass

class OldFormatError(Exception):
	"Old file format passed to function that expects new format"
	pass

class BadFormatError(Exception):
	"Can't parse file format"
	pass

class EmptyFileError(Exception):
	"File contains no records"
	pass

class ZgoubiRunError(Exception):
	"Zgoubi hit a fatal error"
	pass
