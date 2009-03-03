#!/usr/bin/env python
import sys
import struct

def la_lc(la):
	la.sort()
	lc = []
	next_byte=0
	for l in la:
		gap = l[0] - next_byte
		if gap < 0:
			raise ValueError, "Overlap at:%s"%l
		elif gap > 0:
			for g in xrange(gap):
				lc.append('x')
				next_byte +=1
			
		
		lc.append(l[1])
		#print "#", l, next_byte
		#print ''.join(lc)
		next_byte+=struct.calcsize(l[1])
	return lc
			
		
def lc_la(lc):
	la = []
	next_byte=0
	for l in lc:
		if l == 'x':
			next_byte += 1
			continue
		else:
			la.append([next_byte, l])
			next_byte += struct.calcsize(l)
	return la


def file_la(path):
	fh = open(path)
	la = []
	for l in fh:
		if l.strip() == '' or l.startswith('#'):
			continue
		bits = l.split()
		la.append([int(bits[0]), bits[1], bits[2]])
	return la



if __name__ == '__main__':
	fmt_la1 = [[0,'c','name'],[1,'d'],[10,'i'], [15,'3c']]

	lc1 =  la_lc(fmt_la1)
	print lc1
	la1 = lc_la(lc1)
	print la1
	
