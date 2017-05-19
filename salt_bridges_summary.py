#!/usr/bin/env python

import os,glob,re,time
import numpy as np

def uniquify(array):
    """Get unique rows in an array. Duplicated from art_ptdins.py"""
    #---contiguous array trick
    alt = np.ascontiguousarray(array).view(
        np.dtype((np.void,array.dtype.itemsize*array.shape[1])))
    unique,idx,counts = np.unique(alt,return_index=True,return_counts=True)
    #---sort by count, descending
    idx_sorted = np.argsort(counts)[::-1]
    return idx[idx_sorted],counts[idx_sorted]

def count_hydrogen_bonds_redux(bonds,obs,resnames,nmols,resnames_PIP2):
	"""
	Innermost loop of count_hydrogen_bonds function in plot-hydrogen_bonding.py for use on both 
	hydrogen bonding data and the salt bridge data. 
	"""
	#---filter out intralipid hydrogen bonds
	resids_d = bonds[:,1].astype(int)
	resids_a = bonds[:,4].astype(int)
	inter_inds = np.where(resids_d!=resids_a)[0]
	inter = bonds[inter_inds]
	#---catalog bonding pairs by residue
	resnames_combos = np.transpose((inter[:,0],inter[:,3]))
	inds,counts = uniquify(resnames_combos)
	pairs = resnames_combos[inds]
	#---discard sol
	lipid_only = np.where(np.all(pairs!='SOL',axis=1))[0]
	nmols = dict(zip(resnames,nmols))
	#---cut cholesterol in half because it is in both leaflets and POPC does not make hbonds
	nmols_leaflet = 400 
	if 'CHL1' in nmols: nmols['CHL1'] /= 2.0
	#---get the proportions relative to combos
	#---! CHECK THE MEANING OF THIS NUMBER, COMBINATORICALLY-WISE
	#---subsample the obs matrix (frames by combos) to pick out each unique resname combo
	#---summing over the first axis adds up all unique instances of a particular combo
	counts_by_pair = [obs[:,np.where(resnames_combos==p)[0]].sum(axis=1) for p in pairs[lipid_only]]
	#---counts are normalized by the number of each species (note that 400 choose 2 is 79800)
	norm_combos = lambda i,j: (nmols[i]*nmols[j])
	props_mean = np.array([float(counts_by_pair[k].mean())/norm_combos(i,j)
		for k,(i,j) in enumerate(pairs[lipid_only])])
	props_std = np.array([float(counts_by_pair[k].std())/norm_combos(i,j)
		for k,(i,j) in enumerate(pairs[lipid_only])])
	#---convert pairs to PtdIns alias
	pip2_alias = lambda x: 'PtdIns' if x in resnames_PIP2 else x
	aliased_names = np.array([[pip2_alias(i) for i in j] for j in pairs[lipid_only]])
	result = dict(zip([tuple(i) for i in aliased_names],props_mean))
	result_err = dict(zip([tuple(i) for i in aliased_names],props_std))
	return result,result_err

def salt_bridges_summary(**kwargs):

	"""
	Generic hydrogen bonding code.
	Revamped on 2017.4.28 to generate a more uniform data structure.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	resnames_PIP2 = work.vars['selectors']['resnames_PIP2']
	calc = kwargs['calc']

	#---replicate the count_hydrogen_bonds function in plot-hydrogen_bonding.py
	dat = kwargs['upstream']['salt_bridges']
	bonds,obs = dat['bonds'],dat['observations']
	mean_hbonds,std_hbonds = count_hydrogen_bonds_redux(bonds,obs,dat['resnames'],dat['nmols'],resnames_PIP2)

	#---extension for salt bridges
	#---! sorry for the inconsistent naming scheme
	bonds,obs = dat['bonds_salt'],dat['counts_per_frame_salt']
	mean_salt,std_salt = count_hydrogen_bonds_redux(bonds,obs,dat['resnames'],dat['nmols'],resnames_PIP2)

	import ipdb;ipdb.set_trace()