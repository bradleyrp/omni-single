#!/usr/bin/env python

"""
Current optimized version of the curvature-undulation coupling method.
"""

import os
import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from codes.readymade_meso_v1 import import_membrane_mesh

import time
import numpy as np
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory

from codes.mesh import makemesh_regular
from base.timer import checktime
from base.tools import status,framelooper
from base.store import datmerge

from codes.readymade_meso_v1 import import_membrane_mesh

def import_readymade_meso_v1_membrane(**kwargs):
	"""
	Compute bilayer midplane structures for studying undulations.
	Adapted from `undulations.py`.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---import mesh points
	points = import_membrane_mesh(sn=sn,calc=calc,work=work)
	#---ensure there are the same number of points
	points_shapes = list(set([p.shape for p in points]))
	if len(points_shapes)!=1: 
		raise Exception('some frames have a different number of points: %s'%points_shapes)
	else: npoints,ncols = points_shapes[0]
	if ncols!=4: raise Exception('expecting 4-column input on incoming data')
	#---with a consistent number of points everything is an array
	points = np.array(points)[:,:,:3]
	#---previously checked that the minimum points were identically zero but this was not always true
	#---box vectors are just the maximum points
	#---! check that this assumption makes sense
	vecs = points.max(axis=1)[:,:3]
	#---debug the shapes in 3D
	if False:
		from codes import review3d
		fr = 0
		review3d.pbcbox(vecs[fr])
		review3d.review3d(points=[points[fr][:,:3]],radius=10)
	grid_spacing = calc['specs']['grid_spacing']
	nframes = len(points)
	#---choose grid dimensions
	grid = np.array([round(i) for i in np.mean(vecs,axis=0)/grid_spacing])[:2]
	#---compute in parallel
	start = time.time()
	mesh = [[]]
	mesh[0] = Parallel(n_jobs=work.nprocs,verbose=0)(
		delayed(makemesh_regular,has_shareable_memory)(points[fr],vecs[fr],grid)
		for fr in framelooper(nframes,start=start,text='frame'))
	checktime()
	#---pack
	attrs,result = {},{}
	result['mesh'] = np.array(mesh)
	result['grid'] = np.array(grid)
	result['nframes'] = np.array(nframes)
	result['vecs'] = vecs
	attrs['grid_spacing'] = grid_spacing
	#---introduce a dated validator string to ensure that any changes to the pipeline do not overwrite
	#---...other data and are fully propagated downstream
	attrs['validator'] = '2017.08.16.1930'
	return result,attrs	
