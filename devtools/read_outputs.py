#!/usr/bib/env python
#import io2 as io

from zgoubi import io

io.store_def_all() ; exit()

af =  io.read_file("ascii.fai")
#print af
print af.shape
print af[0]

bf =  io.read_file("binary.fai")
#print bf
print bf.shape
print bf[0]

#for n, p in zip( af[0], bf[0]): print p



print "plt"

af =  io.read_file("ascii.plt")
#print af
print af.shape
print af[0]

bf =  io.read_file("binary.plt")
#print bf
print bf.shape
print bf[0]

#for p in zip( af[0], bf[0]): print p





import timeit

#print timeit.timeit('io.define_file("ascii.fai")', setup="import io2 as io", number=10)
#print timeit.timeit('io.define_file("ascii.plt")', setup="import io2 as io", number=10)
#print timeit.timeit('io.define_file("ascii.fai", allow_lookup=False)', setup="import io2 as io", number=100)

#print timeit.timeit('io.read_file("ascii.fai")', setup="import io2 as io", number=100)
#print timeit.timeit('io.read_file("binary.fai")', setup="import io2 as io", number=100)

