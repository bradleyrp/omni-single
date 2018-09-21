#!/usr/bin/python

import time
from numpy import *
import numpy as np
import MDAnalysis
from joblib import Parallel,delayed
# from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
import re

from codes.binding import *

def ion_binding(grofile,trajfile,**kwargs):

	"""
	Analyze bound ion distances to the nearest lipids.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	nrank = kwargs['calc']['specs']['nrank']
	#---! note that the parallel code is severely broken and caused wierd spikes in the distances!!!
	#! on 2018.07.13 trying to revive parallel
	compute_parallel = True
	
	#---prepare universe	
	#! is this deprecated? grofile,trajfile = [work.slice(sn)['current']['all'][i] for i in ['gro','xtc']]
	#! uni = MDAnalysis.Universe(work.postdir+grofile,work.postdir+trajfile)
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.

	#---compute masses by atoms within the selection
	sel_lipids = uni.select_atoms(' or '.join('resname %s'%r 
		for r in work.vars['selectors']['resnames_lipid_chol']))
	sel_ions = uni.select_atoms(work.vars['selectors']['cations'])

	#---load lipid points into memory
	trajectory_ions = zeros((nframes,len(sel_ions),3))
	trajectory_lipids = zeros((nframes,len(sel_lipids),3))
	vecs = zeros((nframes,3))
	for fr in range(nframes):
		status('loading frame',tag='load',i=fr,looplen=nframes)
		uni.trajectory[fr]
		#! mdanalysis removed coordinates(), using positions with cast to be sure
		trajectory_lipids[fr] = np.array(sel_lipids.positions)/lenscale
		trajectory_ions[fr] = np.array(sel_ions.positions)/lenscale
		vecs[fr] = sel_lipids.dimensions[:3]/lenscale

	monolayer_indices = kwargs['upstream']['lipid_abstractor']['monolayer_indices']
	resids = kwargs['upstream']['lipid_abstractor']['resids']
	monolayer_residues = [resids[where(monolayer_indices==mn)[0]] for mn in range(2)]
	group_lipid = uni.select_atoms(' or '.join(['resid '+str(i) for mononum in range(2) 
		for i in monolayer_residues[mononum]]))
	lipid_resids = array([i.resid for i in group_lipid])
	if work.meta[sn]['composition_name'] != 'asymmetric': lipid_resid_subselect = slice(None,None)
	#---hack to account for asymmetric bilayer by analyzing only the first (top) monolayer 
	else: lipid_resid_subselect = where([i.resid in monolayer_residues[0] for i in group_lipid])[0]

	#---parallel partnerfinder
	start = time.time()
	lipid_resids = array([i.resid for i in group_lipid])
	if not compute_parallel:
		incoming = []
		start = time.time()
		for fr in range(nframes):
			status('frame',i=fr,looplen=nframes,start=start,tag='compute')
			ans = partnerfinder(trajectory_lipids[fr],trajectory_ions[fr],vecs[fr],
				lipid_resids,nrank,includes=lipid_resid_subselect)
			incoming.append(ans)
	else:
		incoming = Parallel(n_jobs=4,verbose=0,require='sharedmem')(
			delayed(partnerfinder)
				(trajectory_lipids[fr],trajectory_ions[fr],vecs[fr],
					lipid_resids,nrank,includes=lipid_resid_subselect)
			for fr in framelooper(nframes,start=start))
	n_ion_atoms,n_lipid_atoms = len(sel_ions),len(sel_lipids)
	lipid_distances = zeros((nframes,n_ion_atoms,nrank))
	partners_atoms = zeros((nframes,n_ion_atoms,nrank))

	#---unpack
	start = time.time()
	for fr in range(nframes):
		status('[UNPACK] frame',i=fr,looplen=nframes,start=start)
		lipid_distances[fr] = incoming[fr][0]		
		partners_atoms[fr] = incoming[fr][1]
		
	result,attrs = {},{}
	attrs['nrank'] = nrank
	result['lipid_distances'] = lipid_distances
	result['partners_atoms'] = partners_atoms.astype(int)
	result['names'] = array([i.name for i in group_lipid])
	result['resnames'] = array([i.resname for i in group_lipid])
	result['resids'] = array([i.resid for i in group_lipid])
	result['nframes'] = array(nframes)
	return result,attrs
