#!/usr/bin/python

import numpy as np
import MDAnalysis
from codes.diffusion import *

def diffusion_lipids(**kwargs):

	"""
	Compute bilayer midplane structures for studying undulations.
	"""

	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	dat = kwargs['upstream']['lipid_abstractor']
	resnames = dat['resnames']
	monolayer_indices = dat['monolayer_indices']
	nframes = dat['nframes']	
	diffusion_skip = calc['specs']['diffusion_skip']
	trajectory = dat['points']
	vecs = dat['vecs']

	#---get times 
	#! unfortunately we don't have the timeseries in this version of the calculation
	#! generalize this method so times are easier to get
	uni = MDAnalysis.Universe(kwargs['structure'],kwargs['trajectory'])
	times = []
	for fr in range(len(uni.trajectory)):
		uni.trajectory[fr]
		times.append(uni.trajectory.time)
	times = np.array(times)
	#! times = work.slice(sn)[calc['slice_name']][
	#!   'all' if not kwargs['group'] else kwargs['group']]['timeseries']
	result,attrs = diffusion(trajectory,times,vecs,skip=diffusion_skip)
	attrs['diffusion_skip'] = diffusion_skip
	result['resnames'] = resnames
	return result,attrs

