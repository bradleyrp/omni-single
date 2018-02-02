#!/usr/bin/env python

import time
import numpy as np
import MDAnalysis
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
from base.timer import checktime
import codes.mesh

def lipid_abstractor(grofile,trajfile,**kwargs):

	"""
	LIPID ABSTRACTOR
	Reduce a bilayer simulation to a set of points.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	parallel = kwargs.get('parallel',False)
	#---prepare universe	
	#---note that the universe throws a UserWarning on coarse-grained systems
	#---...which is annoying to elevate to error stage and handled below without problems
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.
	#---select residues of interest
	selector = kwargs['calc']['specs']['selector']
	nojumps = kwargs['calc']['specs'].get('nojumps','')

	#---center of mass over residues
	if 'type' in selector and selector['type'] == 'com' and 'resnames' in selector:
		resnames = selector['resnames']
		selstring = '('+' or '.join(['resname %s'%i for i in resnames])+')'
	elif 'type' in selector and selector['type'] == 'select' and 'selection' in selector:
		if 'resnames' not in selector: raise Exception('add resnames to the selector')
		selstring = selector['selection']
	else: raise Exception('\n[ERROR] unclear selection %s'%str(selector))

	#---compute masses by atoms within the selection
	sel = uni.select_atoms(selstring)
	if len(sel)==0: raise Exception('empty selection')
	mass_table = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'P':30.974,'S':32.065}
	missing_atoms_aamd = list(set([i[0] for i in sel.atoms.names if i[0] not in mass_table]))
	if any(missing_atoms_aamd): 
		print('[WARNING] missing mass for atoms %s so we assume this is coarse-grained'%missing_atoms_aamd)
		#---MARTINI masses
		mass_table = {'C':72,'N':72,'P':72,'S':45,'G':72,'D':72,'R':72}
		missing_atoms_cgmd = list(set([i[0] for i in sel.atoms.names if i[0] not in mass_table]))
		if any(missing_atoms_cgmd):
			raise Exception('we are trying to assign masses. if this simulation is atomistic then we are '+
				'missing atoms "%s". if it is MARTINI then we are missing atoms "%s"'%(
				missing_atoms_aamd,missing_atoms_cgmd))
		else: masses = np.array([mass_table[i[0]] for i in sel.atoms.names])
	else: masses = np.array([mass_table[i[0]] for i in sel.atoms.names])
	
	resids = sel.resids
	#---create lookup table of residue indices
	if len(sel.resids)==len(np.unique(sel.resids)):
		divider = [np.where(resids==r) for r in np.unique(resids)]
	#---note that redundant residue numbering requires special treatment
	else:
		if ('type' in selector) and (selector['type'] in ['com','select']) and ('resnames' in selector):
			#---note that MDAnalysis sel.residues *cannot* handle redundant numbering
			#---note also that some test cases have redundant residues *and* adjacent residues with
			#---...the same numbering. previously we tried a method that used the following sequence:
			#---......divider = [np.where(np.in1d(np.where(np.in1d(
			#---..........uni.select_atoms('all').resnames,resnames))[0],d))[0] for d in divider_abs]
			#---...however this method is flawed because it uses MDAnalysis sel.residues and in fact
			#---...since it recently worked, RPB suspects that a recent patch to MDAnalysis has broken it
			#---note that rpb started a method to correct this and found v inconsistent MDAnalysis behavior
			#---the final fix is heavy-handed: leaving nothing to MDAnalysis subselections
			allsel = uni.select_atoms('all')
			lipids = np.where(np.in1d(allsel.resnames,np.array(selector['resnames'])))[0]
			resid_changes = np.concatenate(([-1],np.where(
				allsel[lipids].resids[1:]!=allsel[lipids].resids[:-1])[0]))
			residue_atomcounts = resid_changes[1:]-resid_changes[:-1]
			#---get the residue names for each lipid in our selection by the first atom in that lipid
			#---the resid_changes is prepended with -1 in the unlikely (but it happened) event that 
			#---...a unique lipid leads this list (note that a blase comment dismissed this possibility at 
			#---...first!) and here we correct the resnames list to reflect this. resnames samples the last 
			#---...atom in each residue from allsel 
			resnames = np.concatenate((allsel[lipids].resnames[resid_changes[1:]],
				[allsel[lipids].resnames[-1]]))
			guess_atoms_per_residue = np.array(zip(resnames,residue_atomcounts))
			#---get consensus counts for each lipid name
			atoms_per_residue = {}
			for name in np.unique(resnames):
				#---get the most common count
				counts,obs_counts = np.unique(guess_atoms_per_residue[:,1][
					np.where(guess_atoms_per_residue[:,0]==name)[0]].astype(int),return_counts=True)
				atoms_per_residue[name] = counts[obs_counts.argmax()]
			#---faster method
			resid_to_start = np.transpose(np.unique(allsel.resids,return_index=True))
			resid_to_start = np.concatenate((resid_to_start,[[resid_to_start[-1][0]+1,len(lipids)]]))
			divider = np.array([np.arange(i,j) 
				for i,j in np.transpose((resid_to_start[:,1][:-1],resid_to_start[:,1][1:]))])
			#---make sure no molecules have the wrong number of atoms
			if not set(np.unique([len(i) for i in divider]))==set(atoms_per_residue.values()):
				status('checking lipid residue indices the careful way',tag='warning')
				#---the following method is slow on large systems. we use it when the fast method above fails
				#---iterate over the list of lipid atoms and get the indices for each N-atoms for each lipid
				counter,divider = 0,[]
				while counter<len(lipids):
					status('indexing lipids',
						i=counter,looplen=len(lipids),tag='compute')
					#---until the end, get the next lipid resname
					this_resname = allsel.resnames[lipids][counter]
					if selector['type']=='select':		
						#---the only way to subselect here is to select on each divided lipid (since 
						#---...the procedure above has correctly divided the lipids). we perform the 
						#---...subselection by pivoting over indices
						#---! this method needs checked
						this_inds = np.arange(counter,counter+atoms_per_residue[this_resname])
						this_lipid = allsel[lipids][this_inds]
						this_subsel = np.where(np.in1d(this_lipid.indices,this_lipid.select_atoms(
							selector['selection']).indices))[0]
						divider.append(this_inds[this_subsel])
					else: divider.append(np.arange(counter,counter+atoms_per_residue[this_resname]))
					counter += atoms_per_residue[this_resname]
				#---in the careful method the sel from above is broken but allsel[lipids] is correct
				sel = allsel[lipids]
				masses = np.array([mass_table[i[0]] for i in sel.atoms.names])
		else: raise Exception('residues have redundant resids and selection is not the easy one')

	#---load trajectory into memory	
	trajectory,vecs = [],[]
	for fr in range(nframes):
		status('loading frame',tag='load',i=fr,looplen=nframes)
		uni.trajectory[fr]
		trajectory.append(sel.positions/lenscale)
		vecs.append(sel.dimensions[:3])
	vecs = np.array(vecs)/lenscale

	checktime()
	#---parallel
	start = time.time()
	if parallel:
		coms = Parallel(n_jobs=work.nprocs,verbose=0)(
			delayed(codes.mesh.centroid)(trajectory[fr],masses,divider)
			for fr in framelooper(nframes,start=start))
	else:
		coms = []
		for fr in range(nframes):
			status('computing centroid',tag='compute',i=fr,looplen=nframes,start=start)
			coms.append(codes.mesh.centroid(trajectory[fr],masses,divider))

	#---identify leaflets
	status('identify leaflets',tag='compute')
	separator = kwargs['calc']['specs'].get('separator',{})
	leaflet_finder_trials = separator.get('trials',3)
	#---preselect a few frames, always including the zeroth
	selected_frames = [0]+list(np.random.choice(np.arange(1,nframes),leaflet_finder_trials,replace=False))
    #---alternate lipid representation is useful for separating monolayers
	if 'lipid_tip' in separator:
		tip_select = separator['lipid_tip']
		sel = uni.select_atoms(tip_select)
		atoms_separator = []
		for fr in selected_frames:
			uni.trajectory[fr]
			atoms_separator.append(sel.positions/lenscale)
	#---default is to use the centers of mass to distinguish leaflets
	else: atoms_separator = [coms[fr] for fr in selected_frames]
	#---pass frames to the leaflet finder, which has legacy and cluster modes
	leaflet_finder = codes.mesh.LeafletFinder(
		atoms_separator=atoms_separator,
		#---pass along the corresponding vectors for topologize
		vecs=[vecs[i] for i in selected_frames],
		cluster=separator.get('cluster',False),
		cluster_neighbors=separator.get('cluster_neighbors',None),
		topologize_tolerance=separator.get('topologize_tolerance',None))
	#---get the indices from the leaflet finder
	monolayer_indices = leaflet_finder.monolayer_indices

	checktime()
	coms_out = np.array(coms)
	#---remove jumping in some directions if requested
	if nojumps:
		nojump_dims = ['xyz'.index(j) for j in nojumps]
		nobjs = coms_out.shape[1]
		displacements = np.array([(coms_out[1:]-coms_out[:-1])[...,i] for i in range(3)])
		for d in nojump_dims:
			shift_binary = (np.abs(displacements)*(1.-2*(displacements<0))/
				(np.transpose(np.tile(vecs[:-1],(nobjs,1,1)))/2.))[d].astype(int)
			shift = (np.cumsum(-1*shift_binary,axis=0)*np.transpose(np.tile(vecs[:-1,d],(nobjs,1))))
			coms_out[1:,:,d] += shift

	#---pack
	attrs,result = {},{}
	attrs['selector'] = selector
	attrs['nojumps'] = nojumps
	result['resnames'] = np.array(sel.residues.resnames)
	result['monolayer_indices'] = np.array(monolayer_indices)
	result['vecs'] = vecs
	result['nframes'] = np.array(nframes)
	result['points'] = coms_out
	result['resids'] = np.array(np.unique(resids))
	result['resids_exact'] = resids
	attrs['separator'] = kwargs['calc']['specs']['separator']
	return result,attrs	
