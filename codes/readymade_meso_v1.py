#!/usr/bin/env python

"""
IMPORT SAMANEH'S DATA!
"""

import os,sys,glob
from base.tools import status

def frame_to_mesh():
	"""
	Convert a list of XYZ points into the mesh format with several checks for consistency.
	"""
	xs,ys = [np.unique(i) for i in dat[:,:2].T]


def import_membrane_mesh(**kwargs):
	"""
	"""
	sn = kwargs.pop('sn',None)
	calc = kwargs.pop('calc',None)
	work = kwargs.pop('work',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---location data can be found in the slices dictionary
	#---! note that the slice name is hard-coded here: "current"
	location = work.slices[sn]['readymade_meso_v1']['current']
	#---get all of the filenames
	fns = sorted(glob.glob(os.path.join(location['path'],location['directory'],
		location['membrane_xyz_glob'])))
	#---! save timestamps here if desired. ensure synchronicity with the nanogel inputs
	#---read each file
	nframes = len(fns)
	for fnum,fn in enumerate(fns):
		status('reading %s'%fn,i=fnum,looplen=nframes,tag='load')
		with open(fn) as fp: text = fp.read()
		lines = text.splitlines()
		#---first line is metadata
		#---! should we drop the first line?
		topline = lines[0]
		#---data is xyz in columns plus fourth column
		dat = np.array([[float(j) for j in line.split()] for line in lines[1:]])
		#---! UNDER CONSTRUCTION. need to load these into a mesh object
		import ipdb;ipdb.set_trace()
	reform = {}
	#---reformulate the data in the manner InvestigateCurvature expects i.e. upstream omnicalc format
	reform['nframes'] = len(mesh)
	#---! UNDER CONSTRUCTION
	reform['grid_spacing'] = np.array(-1.0)
	reform['vecs'] = np.array([-1.0,-1.0,-1.0])
	reform['grid'] = np.array([nx,ny])

def import_curvature_inducer_points(**kwargs):
	"""
	"""
	calc = kwargs.pop('calc',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	import ipdb;ipdb.set_trace()
