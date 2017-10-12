#!/usr/bin/env python

"""
Analyzing undulation spectra.
"""

import numpy as np

def undulations_to_entropy(*args,**kwargs):
	"""
	???
	"""
	kB = 1.38E-23
	temp = 300.00
	mass = 5.5E-22
	Plank = 6.3E-34
	Area_mem = 0.25E-12
	my_const = 3.1415

	if len(args)!=2: raise Exception('arguments list must be two items long: qs and energies')
	high_cutoff = kwargs.pop('high_cutoff',1.0)
	kappa,sigma = kwargs.pop('kappa'),kwargs.pop('sigma')
	if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
	#---the wavevectors and corresponding energies are the arguments
	qs,energy = args

	#---! ryan tested some code here
	#---get indices where qs are larger than zero
	qs_valid = np.where(qs>0)[0]
	#---a list comprehension is the same as a for loop
	freq = np.array([np.sqrt(Area_mem*(kappa*q**4+sigma*q**2))/np.pi for q in qs[qs_valid]])
	#---alternately you can do a vectorized form without any explicit loop
	freq_alt = np.sqrt(Area_mem*(kappa*qs[qs_valid]**4+sigma*qs[qs_valid]**2))/np.pi
	oper = (freq*Plank/(kB*temp)).astype(np.float64)
	#---! it looks like np.exp(oper) is just 1.0s here
	Entropy_term = np.sum((oper/(np.exp(oper)-1.0))-np.log(1.0-np.exp(-oper)))
	Entropy = Entropy_term*kB*6.022140857E23 
	#---! the above method also gets info because oper is too small
	#---! also beware shitespace issue below. due to tabs/spaces you have a double (nested) loop below, not 2 for loops

	freq, oper = np.zeros(qs.shape,'d'), np.zeros(qs.shape,'d')
	nqs = len(qs)
	for j in range(0,nqs):
		if qs[j] > 0.0:
                   freq[j] = np.sqrt(Area_mem*(kappa*qs[j]**4+sigma*qs[j]**2))/np.pi
                   oper[j] = freq[j]*Plank/(kB*temp)

        Entropy_term = 0.0
        for j in range (0,nqs):
            if qs[j] > 0.0:
	       Entropy_term += (oper[j]/(np.exp(oper[j])-1.0))-np.log(1.0-np.exp(-oper[j]))

        Entropy = Entropy_term*kB*6.022140857E23    # un j/K/mol

	#---! to debug use: import ipdb;ipdb.set_trace()
	return Entropy
