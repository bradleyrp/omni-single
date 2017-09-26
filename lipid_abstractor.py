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
		selstring = selector['selection']
	else: raise Exception('\n[ERROR] unclear selection %s'%str(selector))

	#---compute masses by atoms within the selection
	sel = uni.select_atoms(selstring)
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
		if 'type' in selector and selector['type'] == 'com' and 'resnames' in selector:
			#---note that MDAnalysis sel.residues *cannot* handle redundant numbering
			#---note also that the "ocean" test case has redundant residues *and* adjacent residues with
			#---...the same numbering. previously we tried a method that used the following sequence:
			#---......divider = [np.where(np.in1d(np.where(np.in1d(
			#---..........uni.select_atoms('all').resnames,resnames))[0],d))[0] for d in divider_abs]
			#---...however this method is flawed because it uses MDAnalysis sel.residues and in fact
			#---...since it recently worked, RPB suspects that a recent patch to MDAnalysis has broken it
			#---note that rpb started a method to correct this and found v inconsistent MDAnalysis behavior
			#---the final fix is heavy-handed: leaving nothing to MDAnalysis subselections
			allsel = uni.select_atoms('all')
			lipids = np.where(np.in1d(allsel.resnames,np.array(selector['resnames'])))[0]
			resid_changes = np.where(allsel[lipids].resids[1:]!=allsel[lipids].resids[:-1])[0]
			residue_atomcounts = resid_changes[1:]-resid_changes[:-1]
			#---! this method fails
			if False:
				#---use the number of atoms between resid changes to represent the number of atoms in that residue
				guess_atoms_per_residue = np.array(list(set(zip(allsel.resnames[lipids][resid_changes],
					residue_atomcounts))))
				#---figure out the number of atoms in each lipid type by consensus
				atoms_per_residue = {}
				lipid_resnames_obs = np.unique(allsel.resnames[lipids])
				for name in lipid_resnames_obs:
					subs = np.where(guess_atoms_per_residue[:,0]==name)[0]
					consensus_count = guess_atoms_per_residue[np.argsort(
						guess_atoms_per_residue[subs][:,1].astype(int))[0]][1].astype(int)
					atoms_per_residue[name] = consensus_count
			#---get the residue names for each lipid in our selection by the first atom in that lipid
			#---! this probably skip the last residue but if everything hinges on that well oh well
			resnames = allsel[lipids].resnames[resid_changes]
			guess_atoms_per_residue = np.array(zip(resnames,residue_atomcounts))
			#---get consensus counts for each lipid name
			atoms_per_residue = {}
			for name in np.unique(resnames):
				#---get the most common count
				counts,obs_counts = np.unique(guess_atoms_per_residue[:,1][
					np.where(guess_atoms_per_residue[:,0]==name)[0]].astype(int),return_counts=True)
				atoms_per_residue[name] = counts[obs_counts.argmax()]
			#---iterate over the list of lipid atoms and get the indices for each N-atoms for each lipid type
			counter,divider = 0,[]
			while counter<len(lipids):
				#---until the end, get the next lipid resname
				this_resname = allsel.resnames[lipids][counter]
				divider.append(np.arange(counter,counter+atoms_per_residue[this_resname]))
				counter += atoms_per_residue[this_resname]
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

    #---alternate lipid representation is useful for separating monolayers
	monolayer_cutoff = kwargs['calc']['specs']['separator']['monolayer_cutoff']
	monolayer_cutoff_retry = kwargs['calc']['specs']['separator'].get('monolayer_cutoff',True)
	if 'lipid_tip' in kwargs['calc']['specs']['selector']:
		tip_select = kwargs['calc']['specs']['selector']['lipid_tip']
		sel = uni.select_atoms(tip_select)
		atoms_separator = []
		for fr in range(nframes):
			status('loading lipid tips',tag='load',i=fr,looplen=nframes)
			uni.trajectory[fr]
			atoms_separator.append(sel.positions/lenscale)
	else: atoms_separator = coms
	#---identify monolayers
	status('identify leaflets',tag='compute')
	#---randomly select frames for testing monolayers
	random_tries = 3
	for fr in [0]+[np.random.randint(nframes) for i in range(random_tries)]:
		finder_args = dict(monolayer_cutoff=monolayer_cutoff,monolayer_cutoff_retry=monolayer_cutoff_retry)
		top_tol = kwargs['calc']['specs']['separator'].get('topologize_tolerance',None)
		if top_tol: finder_args.update(topologize_tolerance=top_tol)
		monolayer_indices = codes.mesh.identify_lipid_leaflets(atoms_separator[fr],vecs[fr],**finder_args)
		if type(monolayer_indices)!=bool: break

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
