#!/usr/bin/env python

"""
Current optimized version of the curvature-undulation coupling method.
"""

import os
import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from codes.readymade_meso_v1 import import_curvature_inducer_points,import_membrane_mesh

def curvature_undulation_coupling_dex_XXXXXXXXXXXXXXXXXXXX(**kwargs):
	"""
	Supervise the curvature-undulation coupling calculation.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---the original version of this function was designed to import mesoscale data and beging a curvature
	#---...coupling calculation. however the mesoscale data were not on a regular grid so we moved 
	#---...this one step earlier in the pipeline so we can interpolate the positions here
	#---in this section we retrieve the membrane positions
	#---! COMMENT NEEDS FIXED BELOW!
	#---note that in contrast to the standard undulations.py we retrieve the 
	#---...upstream data in a separate function via readymade_meso_v1. in the standard method, both the
	#---...protein_abstractor and undulations data come into this function via e.g.
	#---...kwargs['upstream']['protein_abstractor'] which gets sent to InvestigateCurvature
	#---note that users who wish to manipulate the incoming data can do so in a custom copy of 
	#---...curvature_undulation_coupling.py (i.e. this script) or in the loader functions
	membrane_mesh = import_membrane_mesh(calc=calc,work=work,sn=sn)
	curvature_inducer_points = import_curvature_inducer_points(calc=calc,work=work,sn=sn)
	#---instantiate the calculation	
	ic = InvestigateCurvature(sn=sn,work=kwargs['workspace'],
		design=kwargs['calc']['specs'].get('design',{}),
		protein_abstractor=curvature_inducer_points,
		undulations=membrane_mesh)
	#---repackage the data
	attrs,result = ic.finding['attrs'],ic.finding['result']
	return result,attrs

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
	#---infer the box vectors from the data
	#if not np.all(points.min(axis=0).min(axis=0)[:2]==np.array([0,0])):
	#	print('minimum points must be the origin')
	#	import ipdb;ipdb.set_trace()
	#---box vectors are just the maximum points
	#---! check that this assumption makes sense
	vecs = points.max(axis=0)[:,:2]
	#"""
	#points[0][np.where(np.abs(points[0][:,2]-151.575007)<10**-3)]
	os.environ['USE_SYSTEM_VTK'] = "OFF"
	from codes import review3d
	review3d.review3d(points=[points[0][:,:3]],radius=10)
	#return
	#looks great
	#"""
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
	return result,attrs	
