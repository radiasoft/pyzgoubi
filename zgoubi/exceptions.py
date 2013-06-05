#!/usr/bin/env python

"Exceptions raised by pyzgoubi"

class NoTrackError(Exception):
	"No particle track"
	pass

class BadLineError(Exception):
	"Line on sufficient for this function"
	pass

class OldFormatError(Exception):
	"Old file format passed to function that expects new format"
	pass

class BadFormatError(Exception):
	"Can't parse file format"
	pass

class EmptyFileError(Exception):
	"File contians not records"
	pass
