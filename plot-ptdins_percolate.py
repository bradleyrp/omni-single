#!/usr/bin/env python

import sys
import time
import numpy as np
from joblib import Parallel,delayed
#! had to touch omni/__init__.py for the following !?
from omni.base.compute_loop import framelooper
from codes.binding import partnerfinder

@autoload(plotrun)
def load():
	"""Load everything for the plot only once."""
	data = plotload(plotname)
	sn = 'membrane-v566'
	dat = data.this[sn]
	"""
	if you want the binding combinator use plot:
	  percolate: 
	    script: plot-ptdins_percolate.py
	    collections: long
	    plotload_version: 2
	    autoplot: True
	    calculation: 
	      ion_binding: {}
	      ion_binding_combinator:
	        # cannot multiplex the bonds here
	        zonecut: 
	          loop:
	            [3.0,3.05,2.3,2.6]
	problem is that the ion binding combinator does not give you the actual ion indices
	#! data.set(name='ion_binding_combinator',select={('zonecut',):2.6})
	"""

if __name__=='__main__':

	if 'uni' not in globals():
		import MDAnalysis
		structure,trajectory = [work.postdir+data.extras[sn]['slice_path']+'.'+i for i in ['gro','xtc']]
		uni = MDAnalysis.Universe(structure,trajectory)
		lenscale = 10.0 # standard MDAnalysis Angstroms in desired units: nm
		nframes = uni.trajectory.n_frames

		nrank = 3
		compute_parallel = True

		# compute masses by atoms within the selection
		sel_lipids = uni.select_atoms(' or '.join('resname %s'%r 
			for r in work.vars['selectors']['resnames_lipid_chol']))
		sel_ions = uni.select_atoms(work.vars['selectors']['cations'])

		# load lipid points into memory
		trajectory_ions = np.zeros((nframes,len(sel_ions),3))
		trajectory_lipids = np.zeros((nframes,len(sel_lipids),3))
		vecs = np.zeros((nframes,3))
		for fr in range(nframes):
			status('loading frame',tag='load',i=fr,looplen=nframes)
			uni.trajectory[fr]
			#! mdanalysis removed coordinates(), using positions with cast to be sure
			trajectory_lipids[fr] = np.array(sel_lipids.positions)/lenscale
			trajectory_ions[fr] = np.array(sel_ions.positions)/lenscale
			vecs[fr] = sel_lipids.dimensions[:3]/lenscale

		#! monolayer_indices = kwargs['upstream']['lipid_abstractor']['monolayer_indices']
		monolayer_indices = dat['monolayer_indices']
		#! resids = kwargs['upstream']['lipid_abstractor']['resids']
		resids = dat['monolayer_indices']
		monolayer_residues = [resids[np.where(monolayer_indices==mn)[0]] for mn in range(2)]
		group_lipid = uni.select_atoms(' or '.join(['resid '+str(i) for mononum in range(2) 
			for i in monolayer_residues[mononum]]))
		lipid_resids = np.array([i.resid for i in group_lipid])

	if 1:
		if work.meta[sn]['composition_name'] != 'asymmetric': lipid_resid_subselect = slice(None,None)
		# hack to account for asymmetric bilayer by analyzing only the first (top) monolayer 
		else: lipid_resid_subselect = np.where([i.resid in monolayer_residues[0] for i in group_lipid])[0]

		#---parallel partnerfinder
		start = time.time()
		lipid_resids = np.array([i.resid for i in group_lipid])
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

		# unpack
		start = time.time()
		for fr in range(nframes):
			status('[UNPACK] frame',i=fr,looplen=nframes,start=start)
			lipid_distances[fr] = incoming[fr][0]		
			partners_atoms[fr] = incoming[fr][1]

