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

vecnorm = lambda vec: vec/np.linalg.norm(vec)
vecangle = lambda v1,v2 : np.arccos(np.dot(vecnorm(v1),vecnorm(v2)))*(180./np.pi)

def hbonder_and_salt_bridges_framewise_DEPRECATED(fr,**kwargs):

	#---DEPRECATED BECAUSE too slow. see one-tree method below

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
	#---! note that ten is excessive! also note that nns entries that are len(pts_fore) when distance is inf
	#---index pairs within the cutoff distance
	close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
	#---list of close donors
	close_donors = nns[aind(close_pairs)]
	close_acceptors = close_pairs[:,0]
	#---procedure is the same as hydrogen_bonding to this point
	#---search for nearby hydrogens and cations

	# ...

	"""
	pseudocode
		we have a giant list of possible hydrogen bonding pairs
		and now we need to scan over hydrogens and cations
		note that the close_donors has some redundancy in it 415000 vs 86000 unique
	"""

	#---unique donor mapping
	close_donor_u,close_donor_u_inds = np.unique(close_donors,return_index=True)
	#---get unique close donors and assemble them as the background points for a new tree
	tree_h = scipy.spatial.ckdtree.cKDTree(pts_back[close_donor_u],boxsize=np.concatenate((vec,vec)))
	#---now find hydrogens that are close to these donors (which are close to other acceptors)
	pts_fore_h = boxstuff(all_h_coords[fr],vec)
	close_h,nns_h = tree_h.query(pts_fore_h,k=10,distance_upper_bound=distance_cutoff/lenscale)
	#---now we have hydrogens close to unique donors, each of which might have more than one possible 
	#---...acceptor, giving is two levels of possible redundancy. each donor has multiple close hydrogens
	#---...and multiple close acceptors. we need to assemble these into a list.
	#---first get the close hydrogens from the query
	close_pairs_h = np.transpose(np.where(np.all((close_h<distance_cutoff/lenscale,close_h>0),axis=0)))
	close_donors_h = nns_h[aind(close_pairs_h)]
	close_hydrogens = close_pairs_h[:,0]
	#---next we distribute these pairs among all donor-acceptor pairs

	"""
	pseudocode
		at this point we have done two lookups:
			get donors and acceptors
			get unique donors from the list above and get their hydrogens
		now for each donor-acceptor pair we must make a list of close hydrogens
		then check the angles
		basics:
			get a donor-acceptor pair
			get the donor number
			get donor-hydrogen pairs from the second tree
			...
		example
			an acceptor-donor pair is: close_donors[0],close_acceptors[0]
			a donor-hydrogen pair is: close_donors_h[0],close_hydrogens[0]
			we want the triplet where the donor is the same
		make the triplet the naive way
			#---allocate N_donors by ten array of hydrogens
			donors_to_hydrogens = [[] for i in pts_back]
			donors_to_acceptors = [[] for i in pts_back]
			for dd,donor in enumerate(close_donors): donors_to_acceptors[donor].append(close_acceptors[dd])
			for dd,donor in enumerate(close_donors_h): donors_to_hydrogens[donor].append(close_hydrogens[dd])
		problem is that this discards the index with the right distances, so we cannot recover the positions
			to check the angles
		try a more fine-grained approach with loops (it's OK if this is slow at first; optimize later)
			take the first donor-acceptor pair
	"""

	#---! needs vectorized
	if False:
		import time
		start_time = time.time()
		master_list = []
		#---loop over donors in the donor-acceptor close list
		for cdind in range(len(close_donors)):
			print('%d/%d'%(cdind,len(close_donors)))
			#---take first observed donor-acceptor pair
			cdind = 0
			#---get the donor and acceptor index for this pair
			dind,acind = close_donors[cdind],close_acceptors[cdind]
			#---find donors on the donor-to-hydrogen search
			# dhinds = np.where(close_donors_h==dind)[0]
			dhinds = np.where(close_donors_h==np.where(close_donor_u==dind)[0][0])[0]
			#---triplets in the box-unwrapped space
			#---! is this the right space do be in. does it handle PBCs?
			for dhind in dhinds:
				if False:
					#---assemble the triplet
					#---note that" all_donor_coords[fr][close_donors_h[dhind]] == pts_back[dind]
					#---!!! current problem is that the hydrogen point doesn't make any sense ...
					#---!!! root cause is that we have pts_back[close_donor_u]
					#---!!! so you can get the right donor-hydrogen in the triplet with: 
					#---!!! ...triplet = pts_back[close_donor_u][dind],pts_fore_h[close_hydrogens[dhind]],pts_fore[acind]
					#---!!! but now the acceptor is wrong
					"""
					construct the triplet first by choosing the donor index
					the donors are reindexed for hydrogens
					to get the right donor we need: np.where(close_donors_h==np.where(close_donor_u==dind)[0][0])[0]
					"""
				triplet = pts_back[dind],pts_fore_h[close_hydrogens[dhind]],pts_fore[acind]
				if any([np.linalg.norm(triplet[i]-triplet[j])>distance_cutoff for i,j in [(0,1),(0,2),(1,2)]]):
					print('distances are incorrect')
					import ipdb;ipdb.set_trace()
				else:
					if vecangle(triplet[0]-triplet[1],triplet[2]-triplet[1])>angle_cutoff:
						master_list.append([dind,acind,dhind])
		print('%.1f minutes'%((time.time()-start_time)/60.))
		import ipdb;ipdb.set_trace()

	"""
	development notes
		the above block obviously needs vectorized
			7.6 minutes to look over close_donors which is 415829 long
			this would be 61.4 days for all of the simulations
		vectorization scheme
			ultimately we want a list of triplets that correspond to the donor-hydrogen-acceptor indices
			everything pivots on the donors, since they find nearby hydrogens (this choice is arbitrary)
			we find it hard to believe that one donor makes more than a couple hydrogen bonds
	"""

	import ipdb;ipdb.set_trace()

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

def hbonder_and_salt_bridges_framewise_DEPRECATED2(fr,**kwargs):

	"""
	Compute hydrogen bonds inside a joblib parallel call, over periodic boundary conditions.
	See hydrogen_bonding.py.
	Uses the old, brute force method for inferring hyrogen bonds.
	Deprecated by the non-brute, explicit hydrogen method from the ITP parser code.
	"""

	#---! hacking this for hydrogen as the subject
	distance_cutoff = 1.7

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
	pts_h = boxstuff(all_h_coords[fr],vec)
	pts_combo = np.concatenate((pts_fore,pts_h))

	#---switching perspective to hydrogen
	if False:
		#---idea: try one-tree method to get close acceptors and hydrogens (also cations)
		#---! why does vec need to be twice as long? (tested that the limits work though)
		try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
		#---KDTree failures are blanked
		except: return {'donors':np.array([]),'acceptors':np.array([])}
		close,nns = tree.query(pts_combo,k=10,distance_upper_bound=distance_cutoff/lenscale)
		#---index pairs within the cutoff distance
		close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
		#---list of close donors and hydrogen/acceptors
		close_donors = nns[aind(close_pairs)]
		close_acceptors = close_pairs[:,0]

	"""
	at this point we have to get angles each donor and the close acceptors both above and below the length
		of the foreground points
	problem is that getting the close pairs in a big list means we need to match combiniations on this list
		according to the donor. and this would be slow
	other possible design changes would be:
		1. do everything from the perspective of the hydrogen/cation
		2. do a trick where you do regular hydrogen bonds with the old method, and separately do the salt 
			bridges. but this would be less elegant and accumulated clumsiness begets errors
		3. discard water
			best to wait until the algorithm works in the most elegant form, and is still too slow to do this
	switching perspective to hydrogen
		once we have the tree, we want to go through all pairs of hydrogen-others and match up the others
			this matching step "sounds" like it's going to be slow without a trick.
		half of the first dimension of nns is the donors, half is acceptors
		try a trick: take the closest donor and acceptor to each and see what happens
			this requires us to reverse the foreground background order
			made donors/acceptors into the back
	"""
	#---switching direction again. hydrogens are foreground
	if False:
		tree = scipy.spatial.ckdtree.cKDTree(pts_h,boxsize=np.concatenate((vec,vec)))
		pts_combo = np.concatenate((pts_fore,pts_back))
		close,nns = tree.query(pts_combo,k=10,distance_upper_bound=distance_cutoff/lenscale)
		#---index pairs within the cutoff distance
		close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
		#---list of close donors
		close_donors = nns[aind(close_pairs)]
		close_acceptors = close_pairs[:,0]

	pts_combo = np.concatenate((pts_fore,pts_back))
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_combo,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return {'donors':np.array([]),'acceptors':np.array([])}
	close,nns = tree.query(pts_h,k=10,distance_upper_bound=distance_cutoff/lenscale)
	#---index pairs within the cutoff distance
	close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
	#---list of close donors
	close_donors = nns[aind(close_pairs)]
	close_acceptors = close_pairs[:,0]
	#---get places where the two closest neighbors are close enough
	#valids = np.where(np.all(close[:,:2]<distance_cutoff/lenscale,axis=1))[0]
	#---get places where there is a donor and an acceptor in either prder
	#np.where(np.any([np.all((nns[:,i]<len(pts_fore),nns[:,j]>=len(pts_fore)),axis=0) for i,j in [(0,1),(1,0)]],axis=0))[0]
	#---get all rows where first two columns of nns have both an acceptor and a donor separately, and they
	#---...are close enough to the hydrogen (the subject)
	valids = np.where(np.all((np.any([np.all((nns[:,i]<len(pts_fore),nns[:,j]>=len(pts_fore)),axis=0) 
		for i,j in [(0,1),(1,0)]],axis=0),np.all(close[:,:2]<distance_cutoff/lenscale,axis=1)),axis=0))[0]
	#---now get the angles for each 
	triplets = np.transpose((pts_combo[nns[valids,0]],pts_h[valids],pts_combo[nns[valids,1]]))
	#---! NOTE THAT THIS IS ONLY THE CLOSEST ... !!! ??? !!!
	#---now that we have the triplets we resume the original method
	v1long,v2long = (triplets[0]-triplets[1],triplets[2]-triplets[1])
	#---periodic correction to the hydrogen distances
	v1,v2 = (v1long-(v1long>vec/2.0)*vec),(v2long-(v2long>vec/2.0)*vec)
	#---combinations of coords_donors indices and assoc_h columns serve as the master indexer
	#---! invalid division warnings for the next two lines
	v1n = v1/np.tile(np.linalg.norm(v1,axis=1),(3,1)).T
	v2n = v2/np.tile(np.linalg.norm(v2,axis=1),(3,1)).T
	#---to avoid large objects try inner1d although I think this is already fast
	#---round the dotted values to avoid arccos problems
	dotted = (v1n.T*v2n.T).sum(axis=0).round(6)
	angles = np.arccos(dotted)*180./np.pi
	#---end original method
	#---figure out proper hydrogen bonds/salt bridges
	#---! invalid division warning
	propers = np.where(angles>angle_cutoff)[0]
	#---compile lists of the indices for each identified bond
	inds_h = valids[propers]
	#---nns has 2 columns where each one could be either the acceptor or the donor and we need to detangle
	#---! wrote this in a manic frenzy. probably should check it. argsort works on guaranteed T/F pairs
	#---acceptors come first. note when writing this I forgot to use valids[propers] as the filter
	inds_acceptors = nns[valids[propers]][((np.arange(len(valids[propers])),
		np.argsort(nns[valids[propers],:2]>=len(pts_fore)).T[0]))]
	inds_donors = nns[valids[propers]][((np.arange(len(valids[propers])),
		np.argsort(nns[valids[propers],:2]>=len(pts_fore)).T[1]))]-len(pts_fore)

	"""
	deprecated again ... see convo with joe+gabriela
	0.0032409074540871445
	ipdb> np.mean(hydrogens_side[incoming[1]['mediators']].resnames=='NA')
	0.0029923587980692163
	ipdb> np.mean(hydrogens_side[incoming[2]['mediators']].resnames=='NA')
	0.0028847004689850457
	"""

	return {'donors':inds_donors,'acceptors':inds_acceptors,'mediators':inds_h}

def hbonder_salt_bridges_framewise(fr,**kwargs):

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

	#---! in the following "close_donors" is e.g. the close ions not the heavy atoms
	close_d,nns_d = tree.query(pts_back,k=10,distance_upper_bound=distance_cutoff/lenscale)
	#---index pairs within the cutoff distance
	close_pairs_d = np.transpose(np.where(np.all((close_d<distance_cutoff/lenscale,close_d>0),axis=0)))
	#---list of close donors
	close_donors = nns_d[aind(close_pairs_d)]
	close_ions_d = close_pairs_d[:,0]

	close_a,nns_a = tree.query(pts_fore,k=10,distance_upper_bound=distance_cutoff/lenscale)
	#---index pairs within the cutoff distance
	close_pairs_a = np.transpose(np.where(np.all((close_a<distance_cutoff/lenscale,close_a>0),axis=0)))
	#---list of close donors
	close_acceptors = nns_a[aind(close_pairs_a)]
	close_ions_a = close_pairs_a[:,0]


	#---! this is too large see optimization below
	if True:

		#---waste some memory to organize these
		cations_to_acceptors = np.zeros((len(pts_cations),len(pts_fore)))
		cations_to_donors = np.zeros((len(pts_cations),len(pts_back)))
		cations_to_acceptors[tuple((close_acceptors,close_ions_a))] += 1
		cations_to_donors[tuple((close_donors,close_ions_d))] += 1

		if False:

			#---this block is an addition to the original method to exclude waters
			#---both tabulators above have a row for each ion and a column for each donor or acceptor
			#---here we filter out the solvent. note that we could do this earlier before running the tree
			#---! do this before the tree? or keep it here for indexing ease?
			cols_lipids_d = np.where(np.in1d(donors_resnames,
				np.array([i for i in np.unique(donors_resnames) if i!='SOL'])))[0]
			cols_lipids_a = np.where(np.in1d(acceptors_resnames,
				np.array([i for i in np.unique(acceptors_resnames) if i!='SOL'])))[0]

			#---total number of salt bridges is here 180/400 which is neat
			#---this line was modified below to exclude waters using the non-solvent columns from above
			valid_bridges = np.where(np.all((cations_to_acceptors[:,cols_lipids_a].sum(axis=1),cations_to_donors[:,cols_lipids_d].sum(axis=1)),axis=0))[0]

		#---! no longer removing water in the block above, but way upstream
		else:

			valid_bridges = np.where(np.all((cations_to_acceptors.sum(axis=1),cations_to_donors.sum(axis=1)),axis=0))[0]


		#### GO HERE

		#---take all combinations of donor and acceptor for each ion, particularly since they might be on the 
		#---...same molecule and later we want to get the intermolecular ones
		master_bridge_listing = []
		for vb in valid_bridges:
			combos = np.array(
				np.meshgrid(np.where(close_acceptors==vb)[0],[vb],np.where(close_donors==vb)[0])
				).T.reshape(-1,3)
			master_bridge_listing.extend(combos)
		master_bridge_listing = np.array(master_bridge_listing)

		if False:

			"""
			pseudocode
				at this point we have two lists of close donors and acceptors to the ions
				some ions have multiple acceptors and donors close to them
					hence the second column of e.g. close_pairs_a may be 1,2, etc not just 0
						because it indexes the position in the nns array that has a valid close ion
				what should we do with multiple bonds?
					in the deprecated code above we take all combiniations and add them to a list
						but rpb is rewriting because that code looks wrong and we need it to be smaller anyway
					see what they are
						ship out the resnames and resids for development
						found that (duh!) there are a lot of bonds to water
							ipdb> acceptors_resnames[np.where(cations_to_acceptors[0])[0]]
							array(['PI2P', 'PI2P', 'SOL', 'SOL', 'SOL'], dtype=object)
							ipdb> donors_resnames[np.where(cations_to_donors[0])[0]]
							array(['SOL', 'SOL', 'SOL', 'SOL', 'SOL', 'SOL'], dtype=object)
					realize now that the newest definition of a salt bridge is less strict
						and along the way we forgot to throw away the solvent
					filtered the columns on the big memory tabulator 
						that matches (joins?) the cations to either the acceptors 
				a note on indexing
					either we exclude water upstream or we exclude it above
						as long as it's not too slow, let's stick with what we have
					excluding it before the tree will cause this code to diverge 
						from the hydrogen bonding code
						which will make more cognitive work
			"""

	#---return a list of acceptor, ion, donor triplets
	return master_bridge_listing
