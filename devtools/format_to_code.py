#!/usr/bin/env python
import sys
import struct

import bin_fmt



def write_dec(fmt_la):
	fmt_lc = bin_fmt.la_lc(fmt_la)
	fmt = '='+''.join(fmt_lc)
	code = ''
	code += "fmt = '%s'\n"%fmt
	code += "chunk_part = chunk[0: struct.calcsize(fmt)]\n"
	code += "bits = struct.unpack(fmt, chunk_part)\n"
	for x in xrange(len(fmt_la)):
		if fmt_la[x][2] != '?' and not fmt_la[x][2].startswith('rec_len'):
			code += "particle['%s']=bits[%i]\n"%(fmt_la[x][2], x)
	return code

if __name__ == '__main__':
	fmt_path = sys.argv[1]
#	bf_path = sys.argv[2]
#	head_len = int(sys.argv[3])
#	chunk_len = int(sys.argv[4])
#	chunk_num = int(sys.argv[5])

	la = bin_fmt.file_la(fmt_path)

#	bf = open(bf_path)
#	bf.seek(head_len + chunk_len*chunk_num)
#	chunk = bf.read(chunk_len)
	print write_dec(la)
	

