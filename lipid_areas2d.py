#!/usr/bin/env python

import time
import numpy as np
import scipy
import scipy.spatial
from base.tools import status

def area_voronoi2d(pts,nmol):
	"""
	Compute the areas of Voronoi cells in 2D from a 3D mesh.
	"""
	vmap = scipy.spatial.Voronoi(pts[:,0:2])
	vertices = [[vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]] for p in range(nmol)]
	#---we zip the vertices to make pairs and then compute their areas, for each cell
	areas = [abs(sum(x1*y2-y1*x2 for (x1, y1), (x2, y2) in 
		zip(vert, vert[1:]+vert[0:1])
		)/2) for vert in vertices]
	return areas
	
def lipid_areas2d(**kwargs):
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
	#---! could not run in parallel?
	start = time.time()
	areas = [[],[]]
	for mn in range(2):
		for fr in range(nframes):
			status('voronoi areas monolayer %s'%mn,i=fr,looplen=nframes,start=start,tag='compute')
			areas[mn].append(area_voronoi2d(dat[i2s(mn,fr,'points')],nmols[mn]))
	#---pack
	attrs,result = {},{}
	result['areas0'] = np.array(areas[0])
	result['areas1'] = np.array(areas[1])
	result['nframes'] = np.array(nframes)
	result['vecs'] = dat['vecs']
	result['monolayer_indices'] = dat['monolayer_indices']
	result['resnames'] = dat['resnames']
	return result,attrs	
