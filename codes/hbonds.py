#!/usr/bin/env python

import numpy as np
import sklearn
import sklearn.neighbors
import scipy
import scipy.spatial

vecnorm = lambda vec: vec/np.linalg.norm(vec)
vecangle = lambda v1,v2 : np.arccos(np.dot(vecnorm(v1),vecnorm(v2)))*(180./np.pi)

def hbonder_framewise(fr,**kwargs):

	"""
	Compute hydrogen bonds inside a joblib parallel call, over periodic boundary conditions.
	See hydrogen_bonding.py.
	Uses the old, brute force method for inferring hyrogen bonds.
	Deprecated by the non-brute, explicit hydrogen method from the ITP parser code.
	"""

	distance_cutoff = kwargs['distance_cutoff']
	angle_cutoff = kwargs['angle_cutoff']
	lenscale = kwargs.get('lenscale',10.0)
	vec = vecs[fr]

	#---ensure that the points are inside the box
	boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
	#---convert back to advanced indexing
	aind = lambda x : tuple(x.T)

	pts_back_unstuffed = all_donor_coords[fr]
	pts_fore_unstuffed = all_acceptor_coords[fr]
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)

	#---! why does vec need to be twice as long? (tested that the limits work though)
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return {'donors':np.array([]),'acceptors':np.array([])}
	close,nns = tree.query(pts_fore,k=10,distance_upper_bound=distance_cutoff/lenscale)

	#---index pairs within the cutoff distance
	close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
	#---list of close donors
	close_donors = nns[aind(close_pairs)]
	close_acceptors = close_pairs[:,0]

	#---now that we have the close matches we can get the subsample of the close hydrogens
	pts_h_unstuffed = all_h_coords[fr][close_donors]
	pts_h = boxstuff(pts_h_unstuffed,vec)
	v1long,v2long = pts_fore[close_acceptors]-pts_h,pts_back[close_donors]-pts_h
	#---periodic correction to the hydrogen distances
	v1,v2 = (v1long-(v1long>vec/2.0)*vec),(v2long-(v2long>vec/2.0)*vec)
	#---combinations of coords_donors indices and assoc_h columns serve as the master indexer
	v1n = v1/np.tile(np.linalg.norm(v1,axis=1),(3,1)).T
	v2n = v2/np.tile(np.linalg.norm(v2,axis=1),(3,1)).T
	#---to avoid large objects try inner1d although I think this is already fast
	#---round the dotted values to avoid arccos problems
	dotted = (v1n.T*v2n.T).sum(axis=0).round(6)
	angles = np.arccos(dotted)*180./np.pi
	valid = np.where(angles>angle_cutoff)[0]
	valid_angles = angles[valid]
	valid_dists = close[aind(close_pairs)][valid]
	valid_donors = close_donors[valid]
	valid_acceptors = close_acceptors[valid]
	#---discard the precise angles and distances since we are only trying to count bonds
	return {'donors':valid_donors,'acceptors':valid_acceptors}

def salt_bridges_framewise(fr,**kwargs):

	"""
	Compute hydrogen bonds inside a joblib parallel call, over periodic boundary conditions.
	See hydrogen_bonding.py.
	Uses the old, brute force method for inferring hyrogen bonds.
	Deprecated by the non-brute, explicit hydrogen method from the ITP parser code.
	"""

	distance_cutoff = kwargs['distance_cutoff']
	lenscale = kwargs.get('lenscale',10.0)
	vec = vecs[fr]

	#---ensure that the points are inside the box
	boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
	#---convert back to advanced indexing
	aind = lambda x : tuple(x.T)

	pts_back_unstuffed = all_donor_coords[fr]
	pts_fore_unstuffed = all_acceptor_coords[fr]
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)
	pts_cations = boxstuff(all_cation_coords[fr],vec)

	try: tree = scipy.spatial.ckdtree.cKDTree(pts_cations,boxsize=np.concatenate((vec,vec)))
	#---sometimes the tree fails?
	except: return []

	#---query both the accceptors and donors against the cations tree (fore/back is just historical name)
	close_d,nns_d = tree.query(pts_back,k=10,distance_upper_bound=distance_cutoff/lenscale)
	close_a,nns_a = tree.query(pts_fore,k=10,distance_upper_bound=distance_cutoff/lenscale)
	#---get the valid close pairs of acceptor-cation and donor-cation where the first column is the index
	#---...in the input to the query while the second is the rank in the nearest neighbors
	close_pairs_d = np.transpose(np.where(np.all((close_d<distance_cutoff/lenscale,close_d>0),axis=0)))
	close_pairs_a = np.transpose(np.where(np.all((close_a<distance_cutoff/lenscale,close_a>0),axis=0)))
	#---get the ions with close donors or acceptors
	close_ions_to_donors = nns_d[aind(close_pairs_d)]
	close_ions_to_acceptors = nns_a[aind(close_pairs_a)]
	#---use some memory to organize these for an easy comparison by row (ions)
	cations_to_acceptors = np.zeros((len(pts_cations),len(pts_fore)))
	cations_to_donors = np.zeros((len(pts_cations),len(pts_back)))
	cations_to_acceptors[tuple((close_ions_to_acceptors,close_pairs_a[:,0]))] += 1
	cations_to_donors[tuple((close_ions_to_donors,close_pairs_d[:,0]))] += 1

	#---previously filtered SOL out here but this had an indexing inconsistency 
	#---...so we moved that filter to salt_bridges.py
	valid_bridges = np.where(np.all((cations_to_acceptors.sum(axis=1),
		cations_to_donors.sum(axis=1)),axis=0))[0]

	#---take all combinations of donor and acceptor for each ion, particularly since they might be on the 
	#---...same molecule and later we want to get the intermolecular ones
	master_bridge_listing = []
	for vb in valid_bridges:
		combos = np.array(
			np.meshgrid(np.where(cations_to_acceptors[vb])[0],[vb],np.where(cations_to_donors[vb])[0])
			).T.reshape(-1,3)
		master_bridge_listing.extend(combos)
	master_bridge_listing = np.array(master_bridge_listing)

	#---! NAMES ARE REAL BAD HERE. this was a mistake
	if False:
		#---! in the following "close_donors" is e.g. the close ions not the heavy atoms
		close_d,nns_d = tree.query(pts_back,k=10,distance_upper_bound=distance_cutoff/lenscale)
		#---index pairs within the cutoff distance
		close_pairs_d = np.transpose(np.where(np.all((close_d<distance_cutoff/lenscale,close_d>0),axis=0)))
		#---list of close donors
		close_donors = nns_d[aind(close_pairs_d)]
		close_ions_d = close_pairs_d[:,0]
		#---now get the acceptors
		close_a,nns_a = tree.query(pts_fore,k=10,distance_upper_bound=distance_cutoff/lenscale)
		#---index pairs within the cutoff distance
		close_pairs_a = np.transpose(np.where(np.all((close_a<distance_cutoff/lenscale,close_a>0),axis=0)))
		#---list of close donors
		close_acceptors = nns_a[aind(close_pairs_a)]
		close_ions_a = close_pairs_a[:,0]

		#---waste some memory to organize these
		cations_to_acceptors = np.zeros((len(pts_cations),len(pts_fore)))
		cations_to_donors = np.zeros((len(pts_cations),len(pts_back)))
		cations_to_acceptors[tuple((close_acceptors,close_ions_a))] += 1
		cations_to_donors[tuple((close_donors,close_ions_d))] += 1

		#---previously filtered SOL out here but this had an indexing inconsistency 
		#---...so we moved that filter to salt_bridges.py
		valid_bridges = np.where(np.all((cations_to_acceptors.sum(axis=1),
			cations_to_donors.sum(axis=1)),axis=0))[0]

		#---take all combinations of donor and acceptor for each ion, particularly since they might be on the 
		#---...same molecule and later we want to get the intermolecular ones
		master_bridge_listing = []
		for vb in valid_bridges:
			combos = np.array(
				np.meshgrid(np.where(close_acceptors==vb)[0],[vb],np.where(close_donors==vb)[0])
				).T.reshape(-1,3)
			master_bridge_listing.extend(combos)
		master_bridge_listing = np.array(master_bridge_listing)

	#---return a list of acceptor, ion, donor triplets
	return master_bridge_listing
