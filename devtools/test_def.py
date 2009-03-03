#!/usr/bin/env python

import zgoubi_makedefs

defs = """
FFAG
FFAG
IL : I
N : I
!N*{
A : I
!}

"""


defs = [x for x in defs.split('\n') if x.strip() != ""]
code = zgoubi_makedefs.make_element_class(defs)

print '###########'
print code

