#!/usr/bin/env python
from __future__ import division
import math

"""All values in 
energy eV,
momentum eV/c,
mass ev/c**2,
speed m/s,
charge elementry charge
"""

SPEED_OF_LIGHT = 299792458 # m/s

def ke_to_gamma(mass, ke):
	return (mass+ke) / mass

def gamma_to_ke(mass, gamma):
	return (gamma-1)*mass

def ke_to_te(mass, ke):
	return mass + ke

def te_to_ke(mass, te):
	return te-mass

def beta_to_gamma(beta):
	return 1 / math.sqrt(1-beta**2)

def gamma_to_beta(gamma):
	return math.sqrt(1 - (1/gamma**2))

def ke_to_beta(mass, ke):
	return gamma_to_beta(ke_to_gamma(mass, ke))

def beta_to_ke(mass, beta):
	return gamma_to_ke(mass, beta_to_gamma(beta))

def te_to_mom(mass, te):
	return math.sqrt(  te**2 - mass**2)

def mom_to_te(mass, mom):
	return math.sqrt(mom**2 + mass**2)

def ke_to_mom(mass, ke):
	return te_to_mom(mass, ke_to_te(mass, ke))

def mom_to_ke(mass, mom):
	return te_to_ke(mass, mom_to_te(mass, mom))

def mom_to_rigidity(mom, charge):
	return mom/ SPEED_OF_LIGHT / charge

def rigidity_to_mom(rigidity, charge):
	return rigidity * SPEED_OF_LIGHT * charge

def ke_to_rigidity(mass, ke, charge):
	return mom_to_rigidity(ke_to_mom(mass, ke),charge)

def rigidity_to_ke(mass, rigidity, charge):
	return mom_to_ke(mass, rigidity_to_mom(rigidity, charge))



