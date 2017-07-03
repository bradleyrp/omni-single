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

	"""
	picking out the right position:
		all_donor_coords[fr][np.where(np.all((donors_side.resids==797,donors_side.names=='O6'),axis=0))[0][0]]
		array([ 11.59600067,   4.88000011,   4.16400003], dtype=float32)
		vmd says:
			Info) molecule id: 0
			Info) trajectory frame: 686
			Info) name: O6
			Info) type: O6
			Info) index: 94085
			Info) residue: 1018
			Info) resname: PI2P
			Info) resid: 797
			Info) chain: X
			Info) segname: 
			Info) x: 26.800001
			Info) y: 135.780014
			Info) z: 98.960007
		looks like we need to do distance
	trying distance
		np.linalg.norm(all_donor_coords[fr][np.where(np.all((donors_side.resids==797,donors_side.names=='O6'),axis=0))[0][0]]-all_donor_coords[fr][np.where(np.all((donors_side.resids==800,donors_side.names=='O5'),axis=0))[0][0]])
		giving large distances like 11A
	back to the main code to confirm distances with VMD
		donors_side.positions[np.where(np.all((donors_side.resids==797,donors_side.names=='O6'),axis=0))[0][0]]
			array([  26.80000114,  135.78001404,   98.96000671], dtype=float32)
			and the match is now exact
		ipdb> np.linalg.norm(donors_side.positions[np.where(np.all((donors_side.resids==797,donors_side.names=='O6'),axis=0))[0][0]]-acceptors_side.positions[np.where(np.all((acceptors_side.resids==800,acceptors_side.names=='O5'),axis=0))[0][0]])
			4.1640506
			and the bond length is exact
	tracking these positions through the code
		realized that the all_donor_coords is indexed by donors_inds
	went back and reindexed inside the hbonds function
		ipdb> all_donor_coords[fr][np.where(np.all((donors_side[donors_inds[0]].resids==797,donors_side[donors_inds[0]].names=='O6'),axis=0))[0][0]]
		array([  2.68000007,  13.57800102,   9.89600086], dtype=float32)
		ipdb> np.linalg.norm(all_acceptor_coords[fr][np.where(np.all((acceptors_side.resids==800,acceptors_side.names=='O5'),axis=0))[0][0]]-all_donor_coords[fr][np.where(np.all((donors_side[donors_inds[0]].resids==797,donors_side[donors_inds[0]].names=='O6'),axis=0))[0][0]])
			0.41640487
	so at this point we have access to the correct points and positions that I see in VMD
		the next steps are the box-stuffer
		followed by the bulk of the routine
	now that we have found the pathological bond in VMD and the top of the core compute function, we can begin to search for it in the results from the top down
		first note that the special 4.16A bond we measured above cannot be found in the close array
	construct two where functions to search for the wrong bond in the close acceptors and donors
		np.where(np.all((acceptors_side[close_acceptors].resids==800,acceptors_side[close_acceptors].names=='O5'),axis=0))[0]
		np.where(np.all((donors_side[donors_inds[0]][close_donors].resids==797,donors_side[donors_inds[0]][close_donors].names=='O6'),axis=0))[0]
		note first that these are repeats since there are several places where each gets close enough to participate in a bond
		we can check this by putting the wheres back into the sub-indexed (twice in the case of donors) x_side arrays
		for now we should combine them into one where function to find the particular combination we are looking for, and then ask why it is getting recorded by checking its apparent angle and distance
	the joint search for the bond we are printing in VMD via vmdmake turns up nothing using the following command
		ipdb> np.where(np.all((donors_side[donors_inds[0]][close_donors].resids==797,donors_side[donors_inds[0]][close_donors].names=='O6',acceptors_side[close_acceptors].resids==800,acceptors_side[close_acceptors].names=='O5'),axis=0))
		(array([19914]),)
		in a sense this is good because it means that this bond is not being counted
	the problem now is working backwards from the end of the hbonder_framewise function to see why this was counted, since the problem must be in the interpretation of the results of that function
		but first, while we are here, let's drop the name requirement and find all bonds between these two resdiues (797 and 800 on frame 686)
		there is one at 19914
		donors_side[donors_inds[0]][close_donors][19914],acceptors_side[close_acceptors][19914] says the bond is between OP52 and OP54
		and this is the one we are plotting already! which is good.
		this suggests that there might be a problem way downstream when we are getting multiple bonds
			since we have direct confirmation that the measurement is only identifying one, then we should go to the vmdmake interface and see why it's plotting two
	continued in plot-ptdins_snapshots after moving the data back into position in post
		starting at the definition of map_red_to_obs
		we find that the top hit is the resname matchup for our 800-797 pairs (as expected)
		interestingly, none of the pairs are between O5 and O6 which means that these two atoms which we are picking up *never* make a close bond
	after plotting the rank 1 (second-ranked) item, I realize we are not at the top-ranked item in the plot I am debugging alongside this in VMD. also, interestingly:
		>>> bonds[map_red_to_obs[0]]
		array([['PI2P', '800', 'O2', 'PI2P', '797', 'O4', 'HO2'],
			['PI2P', '800', 'O2', 'PI2P', '797', 'O5', 'HO2'],
			['PI2P', '800', 'O2', 'PI2P', '797', 'OP43', 'HO2'],
			['PI2P', '800', 'O3', 'PI2P', '797', 'O3', 'HO3'],
			['PI2P', '800', 'O3', 'PI2P', '797', 'O4', 'HO3'],
			['PI2P', '800', 'O3', 'PI2P', '797', 'OP43', 'HO3'],
			['PI2P', '800', 'O6', 'PI2P', '797', 'O13', 'HO6'],
			['PI2P', '800', 'OP52', 'PI2P', '797', 'O13', 'H52'],
			['PI2P', '800', 'OP52', 'PI2P', '797', 'O5', 'H52'],
			['PI2P', '800', 'OP52', 'PI2P', '797', 'O6', 'H52']], 
			dtype='|S4')
		>>> bonds[map_red_to_obs[ranking[rank_num]]]
		array([['PI2P', '797', 'O3', 'PI2P', '800', 'O3', 'HO3'],
			['PI2P', '797', 'O6', 'PI2P', '800', 'O2', 'HO6'],
			['PI2P', '797', 'O6', 'PI2P', '800', 'O5', 'HO6'],
			['PI2P', '797', 'O6', 'PI2P', '800', 'OP52', 'HO6'],
			['PI2P', '797', 'O6', 'PI2P', '800', 'OP53', 'HO6'],
			['PI2P', '797', 'O6', 'PI2P', '800', 'OP54', 'HO6'],
			['PI2P', '797', 'OP52', 'PI2P', '800', 'O5', 'H52'],
			['PI2P', '797', 'OP52', 'PI2P', '800', 'OP52', 'H52'],
			['PI2P', '797', 'OP52', 'PI2P', '800', 'OP53', 'H52'],
			['PI2P', '797', 'OP52', 'PI2P', '800', 'OP54', 'H52']], 
			dtype='|S4')
		>>> bonds_red[bonds_inds][subsel]
		array([['PI2P', '800', 'PI2P', '797'],
			['PI2P', '797', 'PI2P', '800'], ...
		and this means that we can check for the bond again
			and there it is, when we check:
				bonds[map_red_to_obs[ranking[rank_num]]]
			so the above confusion about the bond never being there is wrong and now we are in the right spot
	the kernel of our selection procedure, after whittling things down a bunch, is as follows
		map_red_to_obs[ranking[rank_num]]
		bonds[map_red_to_obs[ranking[rank_num]]]
		the bonds version tells us there is an O6-O5 bond between 797-800 in position 2 the third
		let's see why it's being added to the bond_spec for rendering
	having selected two residues, then gotten the bonds with these residues using map_red_to_obs, we now select a frame
		np.argmax(obs.T[map_red_to_obs[ranking[rank_num]]].sum(axis=0))
		this tells us that the most number of bonds for the 797-800 pairing is 2
		however most are 0,1,2 so there might be redundancies
		the argmax is at frame 686 which we have been intensely debugging
		obs[fr][map_red_to_obs[ranking[rank_num]]]
		np.where(obs[fr][map_red_to_obs[ranking[rank_num]]])
		this tells us that at frame 686 we are observing bonds 2 and 9
		these are the two bonds we are plotting
		>>> bonds[map_red_to_obs[ranking[rank_num]]][np.where(obs[fr][map_red_to_obs[ranking[rank_num]]])[0]]
		array([['PI2P', '797', 'O6', 'PI2P', '800', 'O5', 'HO6'],
			['PI2P', '797', 'OP52', 'PI2P', '800', 'OP54', 'H52']],dtype='|S4')
	at this point it appears that the calculation is correctly recording a single OP52-OP54 bond from the core code at frame 686, *and* that the snapshotter is picking up that bond and another one
		this leaves few places for the error to occur, so hopefully I am getting close
		the problem might be in obs falsely recording that bond at this frame
		which means that the error is in the tabulator
		which will be harder to debug because it requires a full run and not a single frame
	ran the tabulator and pulled up the snapshotter alongside it by moving the post data around so now we can compare
		first we get the row number for our target bonds from the snapshotter
			>>> obs[fr][map_red_to_obs[ranking[rank_num]]]
			array([ 0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  1.])
			>>> map_red_to_obs[ranking[rank_num]]
			array([18593, 18602, 18603, 18604, 18605, 18606, 18616, 18617, 18618, 18619])
		then we check they are the right bonds in the calculation
			ipdb> bonds[18603]
			array(['PI2P', 797, 'O6', 'PI2P', 800, 'O5', 'HO6'], dtype=object)
			ipdb> bonds[18619]
			array(['PI2P', 797, 'OP52', 'PI2P', 800, 'OP54', 'H52'], dtype=object)
		as expected, both of the bonds are there
		now let's check the stuff that went into the tabulator via "incoming[686]"
		perform the same joint search as above on the incoming variable
		but first note that the donors in incoming is the valid donors
		and valid_donors = close_donors[valid]
		the original joint search is
			np.where(np.all((donors_side[donors_inds[0]][close_donors].resids==797,donors_side[donors_inds[0]][close_donors].names=='O6',acceptors_side[close_acceptors].resids==800,acceptors_side[close_acceptors].names=='O5'),axis=0))
		but we are using close donors not valid donors
		let's reconstruct two versions of it
		first we put the donors from incoming in as close_donors
			np.where(np.all((donors_side[donors_inds[0]][incoming[686]['donors']].resids==797,donors_side[donors_inds[0]][incoming[686]['donors']].names=='O6',acceptors_side[incoming[686]['acceptors']].resids==800,acceptors_side[incoming[686]['acceptors']].names=='O5'),axis=0))
			and we get nothing
		now we check the bond that we know is there is in the list
			np.where(np.all((donors_side[donors_inds[0]][incoming[686]['donors']].resids==797,donors_side[donors_inds[0]][incoming[686]['donors']].names=='OP52',acceptors_side[incoming[686]['acceptors']].resids==800,acceptors_side[incoming[686]['acceptors']].names=='OP54'),axis=0))
			also nothing
		more conservatively we try
			np.where(np.all((donors_side[donors_inds[0]][incoming[686]['donors']].resids==797,acceptors_side[incoming[686]['acceptors']].resids==800),axis=0))
			also nothing
		looks like we have the wrong indexing
			(donors_side[incoming[686]['donors']].resids==797).sum() returns 5
		dropping the donors_inds[0] indexing and trying the joint resid-resid check again gives nothing
			np.where(np.all((donors_side[incoming[686]['donors']].resids==797,acceptors_side[incoming[686]['acceptors']].resids==800),axis=0))
		checking them individually
			ipdb> acceptors_side[acceptors_side[incoming[686]['acceptors']].resids==800].resnames
				array(['CHL1', 'CHL1', 'CHL1', 'CHL1', 'CHL1', 'CHL1', 'CHL1', 'CHL1',
				'CHL1', 'CHL1', 'CHL1'], dtype=object)			
		ipdb> donors_side[donors_side[incoming[686]['donors']].resids==797].resnames
		/run/media/rpb/store-omicron/factory/env/envs/py2/lib/python2.7/site-packages/MDAnalysis/core/groups.py:447: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 182387 but corresponding boolean dimension is 42482
			return self._derived_class(self.ix[item], self.universe)
			array(['PI2P', 'SOL', 'SOL', 'SOL', 'SOL'], dtype=object)
		so this is definitely a problem
		trying many different ways of indexing donors_side with donors_inds and incoming[686]['donors']
			but none of them are making sense e.g.
				ipdb> donors_side[donors_side[incoming[686]['donors']].resids==797].resnames
				/run/media/rpb/store-omicron/factory/env/envs/py2/lib/python2.7/site-packages/MDAnalysis/core/groups.py:447: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 182387 but corresponding boolean dimension is 42482
				return self._derived_class(self.ix[item], self.universe)
				array(['PI2P', 'SOL', 'SOL', 'SOL', 'SOL'], dtype=object)
				ipdb> donors_side[donors_inds[0]][donors_side[donors_inds[0]][incoming[686]['donors']].resids==797].resnames
				/run/media/rpb/store-omicron/factory/env/envs/py2/lib/python2.7/site-packages/MDAnalysis/core/groups.py:447: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 122138 but corresponding boolean dimension is 42482
				return self._derived_class(self.ix[item], self.universe)
				array(['SOL', 'SOL'], dtype=object)
				ipdb> donors_side[donors_side[incoming[686]['donors'][donors_inds[0]]].resids==797].resnames
				*** IndexError: index 42482 is out of bounds for axis 1 with size 42482
				ipdb> donors_side[donors_side[donors_inds[0][incoming[686]['donors']]].resids==797].resnames
				array(['SOL', 'SOL'], dtype=object)
				ipdb> donors_side[donors_inds[0]][donors_side[donors_inds[0][incoming[686]['donors']]].resids==797].resnames
				array(['SOL', 'SOL'], dtype=object)
		first we cannot index donors side directly because it is showing hydrogens as heavy donors
			ipdb> donors_side[incoming[686]['donors']].names
			array(['HW2', 'OW', 'OW', ..., 'HW2', 'HW2', 'HW1'], dtype=object)
		it makes more sense to index donors_side via donors_inds before indexing over incoming
			ipdb> donors_side[donors_inds[0]][incoming[686]['donors']].names
			array(['OW', 'OW', 'OW', ..., 'OW', 'OW', 'OW'], dtype=object)
			ipdb> donors_side[donors_inds[1]][incoming[686]['donors']].names
			array(['HW2', 'HW1', 'HW2', ..., 'HW1', 'HW1', 'HW1'], dtype=object)
		this is already the way that donor_cat is constructed, so it seems correct
	snip out the frame's bond list from tabulation just to check it
		bonds_this = tabulation[frame_lims[686]:frame_lims[686+1]]
		bonds_this[np.where(bonds_this[:,1]==797)]
			shows our nefarious evil bond
	note that the above bonds_this call is like the obs that gets added to the total number of observations during the counting step, so the following shows us how the wrong bond gets from tabulation to the final observation list
			ipdb> bonds_this[np.where(bonds_this[:,1]==797)]
			array([['PI2P', 797, 'O3', 'PI2P', 797, 'OP42', 'HO3'],
				['PI2P', 797, 'O3', 'PI2P', 797, 'OP42', 'HO3'],
				['PI2P', 797, 'O6', 'PI2P', 800, 'O5', 'HO6'],
				['PI2P', 797, 'O6', 'PI2P', 800, 'O5', 'HO6'],
				['PI2P', 797, 'OP52', 'PI2P', 800, 'OP54', 'H52']], dtype=object)
			ipdb> np.array([bonds_to_idx[tuple(o)] for o in bonds_this[np.where(bonds_this[:,1]==797)]])
			array([18591, 18591, 18603, 18603, 18619])
		so the remaining questions are basically how does the bond get into the tabulation list for this frame, and more importantly, can we prove it's not in the incoming list?
		also note that some of the resids for SOL are strings but for lipids they are integers?
			this is a trick we play with water, changing everything to residue "1"
	possible cause: the valid_frames has caused us to be off by a certain number of frames, even though I checked +/- 1 frame, it might just be the case that it's off by more
	switched bonds_this to valid_frames[686] which is really frame 689
		ipdb> bonds_this = tabulation[frame_lims[689]:frame_lims[689+1]]
		ipdb> np.array([bonds_to_idx[tuple(o)] for o in bonds_this[np.where(bonds_this[:,1]==797)]])
		array([18591, 18591, 18588, 18588, 18619])
		ipdb> bonds[np.array([bonds_to_idx[tuple(o)] for o in bonds_this[np.where(bonds_this[:,1]==797)]])]
		array([['PI2P', 797, 'O3', 'PI2P', 797, 'OP42', 'HO3'],
		       ['PI2P', 797, 'O3', 'PI2P', 797, 'OP42', 'HO3'],
		       ['PI2P', 797, 'O2', 'PI2P', 797, 'OP53', 'HO2'],
		       ['PI2P', 797, 'O2', 'PI2P', 797, 'OP53', 'HO2'],
		       ['PI2P', 797, 'OP52', 'PI2P', 800, 'OP54', 'H52']], dtype=object)
		now we only see the one correct bond so maybe 686 is really 689
		but I am concerned this might be exactly backwards
		went back to VMD (thankfully I am at home) and jumped to frame 689 and BOTH OF THE BONDS ARE THERE AND THEY LOOK GREAT HUZZAH THIS IS THE PROBLEM. debugged this from 0755-0945
	end result is that all of the hydrogen bonds are wrong so they need to be recomputed with a fix
	what should the fix be?
	"""
	#import ipdb;ipdb.set_trace()
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
