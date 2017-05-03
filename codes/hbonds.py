#!/usr/bin/env python

import numpy as np
import sklearn
import sklearn.neighbors
import scipy
import scipy.spatial

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

	#---how that we have the close matches we can get the subsample of the close hydrogens
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
