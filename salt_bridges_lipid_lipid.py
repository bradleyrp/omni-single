#!/usr/bin/env python

def salt_bridges_lipid_lipid(**kwargs):

	"""
	Identify salt bridges. Mimics the beginning of the hydrogen bond
	"""

	#---unpack
	grofile = kwargs['grofile']
	trajfile = kwargs['trajfile']
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)

	#---settings. distance cutoff is larger for salt bridges than hydrogen bonds
	distance_cutoff = calc['specs']['distance_cutoff']

	#---prepare universe	
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)

	import ipdb;ipdb.set_trace()
