#!/usr/bin/env python

"""
Analyzing undulation spectra.
"""

import numpy as np

def undulations_to_entropy(*args,**kwargs):
	"""
	???
	"""
	my_const = 3.1415
	if len(args)!=2: raise Exception('arguments list must be two items long: qs and energies')
	high_cutoff = kwargs.pop('high_cutoff',1.0)
	kappa,sigma = kwargs.pop('kappa'),kwargs.pop('sigma')
	if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
	#---the wavevectors and corresponding energies are the arguments
	qs,energy = args
	#---! to debug use: import ipdb;ipdb.set_trace()
	return np.sum(qs)
