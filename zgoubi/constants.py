#!/usr/bin/env python
from __future__ import division
from math import log

ELECTRON_MASS = 0.51099892e6 #eV
ELECTRON_CHARGE = -1.60217653e-19 #C
ELECTRON_ANOM_MAG_MOM = (2.0023193043622 -2)/2 # from http://en.wikipedia.org/wiki/G-factor
ELECTRON_MEAN_LIFE = 0 #does not decay

PROTON_MASS = 938.27203e6 #eV
PROTON_CHARGE = 1.602176487e-19
PROTON_ANOM_MAG_MOM = (5.585694701 -2)/2 # from http://en.wikipedia.org/wiki/G-factor
PROTON_MEAN_LIFE = 0

MUON_MASS = 105.65837e6 #eV
MUON_CHARGE = -1.60217653e-19 #C
MUON_ANOM_MAG_MOM = (2.0023318396 -2)/2 # from http://en.wikipedia.org/wiki/G-factor
MUON_MEAN_LIFE = 2.197029e-6  #

SPEED_OF_LIGHT = 299792458


