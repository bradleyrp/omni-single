#!/usr/bin/env python

import time
import numpy as np
import scipy
import scipy.spatial
from base.tools import status
import numpy as np

"""
Plotting the 3D lipid areas was taking a while because it was loading the mesh objects.
This function computes the areas and stores them for faster plotting.
"""

def lipid_areas3d(**kwargs):
	"""
	Compute bilayer midplane structures for studying undulations.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	dat = kwargs['upstream']['lipid_mesh']
	i2s = lambda mn,fr,key : '%d.%d.%s'%(mn,fr,key)
	nmols = [int(dat[i2s(mn,0,'nmol')]) for mn in range(2)]
	nframes = int(dat['nframes'])
	resnames = dat['resnames']
	results = {}

	areas = [np.zeros((nframes,nmols[mn])) for mn in range(2)]
	for mn in range(2):
		for fr in range(nframes):
			simps = dat[i2s(mn,fr,'simplices')]
			vertex_areas = np.array([sum(dat[i2s(mn,fr,'areas')]
				[np.where(np.any(simps==i,axis=1))[0]])/3. for i in range(nmols[mn])])
			areas[mn][fr] = vertex_areas

	#---sort lipids from each monolayer by resname
	imono = dat['monolayer_indices']
	reslist = np.unique(resnames)
	monolayer_residues = [[np.concatenate([np.where(np.where(imono==mn)[0]==i)[0] 
		for i in np.where(resnames==r)[0]]) for r in reslist] for mn in range(2)]
	#---we must reshape below if lipid flip flop
	mu = [np.mean(np.concatenate([np.reshape(areas[mn][:,monolayer_residues[mn][rr]],-1) 
		for mn in range(2) if len(monolayer_residues[mn][rr])!=0])) for rr,r in enumerate(reslist)]
	sigma = [np.std(np.concatenate([np.reshape(areas[mn][:,monolayer_residues[mn][rr]],-1) 
		for mn in range(2) if len(monolayer_residues[mn][rr])!=0])) for rr,r in enumerate(reslist)]

	#---pack
	attrs,result = {},{}
	result['areas0'] = np.array(areas[0])
	result['areas1'] = np.array(areas[1])
	result['nframes'] = np.array(nframes)
	result['vecs'] = dat['vecs']
	result['monolayer_indices'] = dat['monolayer_indices']
	result['resnames'] = dat['resnames']
	result['mu'],result['sigma'] = mu,sigma
	return result,attrs	
