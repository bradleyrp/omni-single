#!/usr/bin/python

import time
from numpy import *
import MDAnalysis
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
from base.timer import checktime
from codes.mesh import *

def protein_abstractor(grofile,trajfile,**kwargs):

	"""
	PROTEIN ABSTRACTOR
	Compute the centroids of proteins in a simulation.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	parallel = kwargs.get('parallel',False)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.
	
	#---get protein coms here
	uni = MDAnalysis.Universe(grofile,trajfile)
	sel = uni.select_atoms(work.vars['selectors']['protein_selection'])
	nprots = work.meta[sn]['nprots']
	beads_per_protein = len(sel.resids)/nprots
	nframes = len(uni.trajectory)
	inds = [arange(i*beads_per_protein,(i+1)*beads_per_protein) for i in range(nprots)]
	trajectory,trajectory_all,vecs = [],[],[]
	start = time.time()
	for fr in range(nframes):
		status('collecting protein centroids',i=fr,looplen=nframes,start=start)
		uni.trajectory[fr]
		#---center of geometry not centroid because masses are all 72 in martini
		pts = sel.coordinates()[array(inds).astype(int)]/lenscale
		pts_mean = pts.mean(axis=0)
		trajectory.append(pts_mean)
		trajectory_all.append(pts)
		vecs.append(sel.dimensions[:3])

	#---pack
	attrs,result = {},{}
	result['resnames'] = array(sel.residues.resnames)
	result['names'] = array(sel.atoms.names)
	result['vecs'] = array(vecs)/lenscale
	result['nframes'] = array(nframes)
	result['points'] = array(trajectory)
	result['points_all'] = array(trajectory_all)
	return result,attrs	

