#!/usr/bin/env python

"""All values in 
energy eV,
momentum eV/c,
mass ev/c**2,
speed c,
charge elementry charge
"""
from __future__ import division
import math

SPEED_OF_LIGHT = 299792458 # m/s

def ke_to_gamma(mass, ke):
	"Convert kinetic energy to lorentz factor"
	return (mass+ke) / mass

def gamma_to_ke(mass, gamma):
	"Convert lorentz factor to kinetic energy"
	return (gamma-1)*mass

def ke_to_te(mass, ke):
	"Convert kinetic energy to total energy"
	return mass + ke

def te_to_ke(mass, te):
	"Convert total energy to kinetic energy"
	return te-mass

def beta_to_gamma(beta):
	"Convert speed to lorentz factor"
	return 1 / math.sqrt(1-beta**2)

def gamma_to_beta(gamma):
	"Convert lorentz factor to speed"
	return math.sqrt(1 - (1/gamma**2))

def ke_to_beta(mass, ke):
	"Convert kinetic energy to speed"
	return gamma_to_beta(ke_to_gamma(mass, ke))

def beta_to_ke(mass, beta):
	"Convert speed to kinetic energy"
	return gamma_to_ke(mass, beta_to_gamma(beta))

def te_to_mom(mass, te):
	"Convert total energy to momentum"
	return math.sqrt(  te**2 - mass**2)

def mom_to_te(mass, mom):
	"Convert momentum to total energy"
	return math.sqrt(mom**2 + mass**2)

def ke_to_mom(mass, ke):
	"Convert kinetic energy to momentum"
	return te_to_mom(mass, ke_to_te(mass, ke))

def mom_to_ke(mass, mom):
	"Convert momentum to kinetic energy"
	return te_to_ke(mass, mom_to_te(mass, mom))

def mom_to_rigidity(mom, charge):
	"Convert momentum to magnetic rigitity"
	return mom/ SPEED_OF_LIGHT / charge

def rigidity_to_mom(rigidity, charge):
	"Convert magnetic rigitity to momentum"
	return rigidity * SPEED_OF_LIGHT * charge

def ke_to_rigidity(mass, ke, charge):
	"Convert kinetic energy to magnetic rigitity"
	return mom_to_rigidity(ke_to_mom(mass, ke), charge)

def rigidity_to_ke(mass, rigidity, charge):
	"Convert magnetic rigitity to kinetic energy"
	return mom_to_ke(mass, rigidity_to_mom(rigidity, charge))



