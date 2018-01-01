#!/usr/bin/env python

import time
import numpy as np
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory

from codes.mesh import makemesh_regular
from base.timer import checktime
from base.tools import status,framelooper
from base.store import datmerge

def undulations(**kwargs):

	"""
	Compute bilayer midplane structures for studying undulations.
	"""

	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	upname = 'lipid_abstractor'
	grid_spacing = calc['specs']['grid_spacing']
	vecs = datmerge(kwargs,upname,'vecs')
	nframes = int(np.sum(datmerge(kwargs,upname,'nframes')))
	trajectory = datmerge(kwargs,upname,'points')
	attrs,result = {},{}
	#---! hacking through error with monolayer separation
	try: monolayer_indices = kwargs['upstream'][upname+'0']['monolayer_indices']
	except: monolayer_indices = kwargs['upstream'][upname]['monolayer_indices']
	#---choose grid dimensions
	grid = np.array([round(i) for i in np.mean(vecs,axis=0)/grid_spacing])[:2]
	#---! removed timeseries from result for new version of omnicalc
	#---parallel
	mesh = [[],[]]
	for mn in range(2):
		start = time.time()
		mesh[mn] = Parallel(n_jobs=work.nprocs,verbose=0)(
			delayed(makemesh_regular,has_shareable_memory)(
				trajectory[fr][np.where(monolayer_indices==mn)],vecs[fr],grid)
			for fr in framelooper(nframes,start=start,text='monolayer %d, frame'%mn))
	checktime()

	#---pack
	result['mesh'] = np.array(mesh)
	result['grid'] = np.array(grid)
	result['nframes'] = np.array(nframes)
	result['vecs'] = vecs
	attrs['grid_spacing'] = grid_spacing
	return result,attrs	
