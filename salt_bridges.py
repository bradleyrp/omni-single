#!/usr/bin/env python

import os,glob,re,time
import numpy as np
import sklearn
import sklearn.neighbors

vecnorm = lambda vec: vec/np.linalg.norm(vec)
vecangle = lambda v1,v2 : np.arccos(np.dot(vecnorm(v1),vecnorm(v2)))*(180./np.pi)

def tabulator(tabulation,valid_frames,obs_by_frames):

	"""
	Given a giant list of atomistic details, count up all of the unique entries for each frame.
	"""

	start_time = time.time()
	status('unique-ifying the tabulated bonds (estimated %ds)'%(len(tabulation)*1.3*10**-6),tag='compute')
	status('note: with 32GB memory, 33M observations works fine, but 46M hits the swap',tag='warning')
	#---note that unique is getting "axis" in np 1.13 but at some point on or before 1.12 they added some 
	#---...kind of a safety check on the following trick for unique rows, which check returns an error
	#---...message: "TypeError: Cannot change data-type for object array." which is solved by forcing 
	#---...the object to a string type. note that this method requires void and not a blank string, which
	#---...some examples will use. this changed must have happened in the <1 week since we wrote 
	#---...the hydrogen bonds code and tested it again on the factory
	#---uniquify the enormous list of all possible hydrogen bonds
	tabulation_reform = tabulation.astype(str)
	tabulation_unique = np.ascontiguousarray(tabulation_reform).view(
		np.dtype((np.void,tabulation_reform.dtype.itemsize*tabulation_reform.shape[1])))
	tabulation_view_unique,idx,counts = np.unique(tabulation_unique,return_index=True,return_counts=True)
	bonds = tabulation[idx]
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')

	start_time = time.time()
	#---preallocate bond counts per frame
	try: counts_per_frame = np.zeros((len(valid_frames),len(idx)))
	except:
		import ipdb;ipdb.set_trace()
	#---hash the binds over the indices
	bonds_to_idx = dict([(tuple(b),bb) for bb,b in enumerate(bonds)])
	frame_lims = np.concatenate(([0],np.cumsum(obs_by_frames)))
	for fr,i in enumerate(frame_lims[:-1]):
		status('counting observations per frame',i=fr,looplen=len(valid_frames),
			tag='compute',start=start_time)
		try: 
			obs = tabulation[frame_lims[fr]:frame_lims[fr+1]]
			counts_per_frame[fr][np.array([bonds_to_idx[tuple(o)] for o in obs])] += 1
		except:
			import ipdb;ipdb.set_trace()
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')
	status('done heavy lifting',tag='compute')
	return bonds,counts_per_frame

def salt_bridges(grofile,trajfile,**kwargs):

	"""
	Identify salt bridges. Mimics the beginning of the hydrogen bond
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)

	#---settings. distance cutoff is larger for salt bridges than hydrogen bonds
	distance_cutoff = calc['specs']['distance_cutoff']

	#---prepare universe	
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	lenscale = 10.
	start_job_time = time.time()

	#---save topology for later
	_,idx,counts = np.unique(uni.residues.resnames,return_index=True,return_counts=True)
	resnames = uni.residues.resnames[np.sort(idx)]
	rescounts = counts[np.argsort(idx)]

	import makeface
	#---get an automacs landscape
	#---! DEV. needs a clone and make to work
	try: mod = makeface.import_remote('amx/amx')
	except: raise Exception('please clone a copy of automacs next to omni in `amx`')
	mod['state'].force_field = 'charmm'
	Landscape = mod['Landscape']
	land = Landscape(cwd='amx/')
	#---use the landscape to get hydrogen bond donors and acceptors for lipids
	hydrogen_bond_ref = {}
	targets = land.objects_by_category('lipid')
	for resname in targets:
		mol = land.itps[land.objects[resname]['fn']][resname]
		#---collect all possible hydrogen bond acceptors
		#---! forbid water here for salt bridge calculation
		acceptor_names = [i['atom'] for i in mol['atoms'] if re.match('^(N|O|S)',i['atom']) and i!='OW']
		h_name = [i['atom'] for i in mol['atoms'] if re.match('^H',i['atom'])]	
		donor_candidates = [(j,k) for j,k in [(int(i['i']),int(i['j'])) 
			for i in mol['bonds']] 
			if any([re.match('^H',mol['atoms'][l-1]['atom']) for l in [j,k]]) 
			and any([mol['atoms'][m-1]['atom'] in acceptor_names for m in [j,k]])]
		donor_names = []
		for d in [(mol['atoms'][i-1]['atom'],mol['atoms'][j-1]['atom']) for i,j in donor_candidates]:
			if d not in donor_names: donor_names.append(d)
		hydrogen_bond_ref[resname] = {'acceptors':acceptor_names,'donors':donor_names}
	#---water-naming is hard-coded 
	hydrogen_bond_ref['water'] = {'donors':[('OW','HW1'),('OW','HW2')],'acceptors':['OW']}
	#---assemble the names
	donors_names = sorted(list(set([m for n in [zip(*i['donors'])[0] 
		for i in hydrogen_bond_ref.values() if i['donors']!=[]] for m in n])))
	hydrogens_names = sorted(list(set([m for n in [zip(*i['donors'])[1] 
		for i in hydrogen_bond_ref.values() if i['donors']!=[]] for m in n])))
	acceptors_names = sorted(list(set([m for n in [i['acceptors'] 
		for i in hydrogen_bond_ref.values() if i!=[]]for m in n])))
	#---generate atom groups
	donors = uni.select_atoms(' or '.join(['name %s'%i for i in donors_names]))
	acceptors = uni.select_atoms(' or '.join(['name %s'%i for i in acceptors_names]))
	hydrogens = uni.select_atoms(' or '.join(['name %s'%i for i in hydrogens_names]))

	#---! not necessary
	if False:

		#---we can either exclude water here or after the KD-Tree. tried the latter and there was index problem
		#---! note that we exclude water by ignoring OW and HW1 and HW2. they should not appear in other mol
		hydrogen_bond_ref = dict([(i,j) for i,j in hydrogen_bond_ref.items() if i!='water'])
		donors_h_pairs = [m for n in [i.get('donors',[]) for i in hydrogen_bond_ref.values()] for m in n]
		donors_h_pairs_flat = list(set([i for j in donors_h_pairs for i in j]))
		#---compared to the hydrogen bonding version, we only take the heavy atom the in donor pairs

		sel_d,sel_h = [uni.select_atoms(' or '.join(['name %s'%i 
			for i in list(set(zip(*donors_h_pairs)[j]))])) for j in range(2)]
		resids = np.unique(np.concatenate([sel_d.resids,sel_h.resids]))
		donors_side = uni.select_atoms(' or '.join(['name %s'%i for i in donors_h_pairs_flat]))
		donors_resids = np.unique(donors_side.resids)

	#---! use acceptor names for both donor and acceptor pairs for this calculation. that is, we do not 
	#---! ...require an intervening hydrogen

	donors_side = uni.select_atoms(' or '.join(['name %s'%i for i in acceptor_names]))
	donors_resids = np.unique(donors_side.resids)

	#---! not necessary
	if False:
		#---identifying residues with both a donor and a corresponding hydrogen
		both = np.zeros((len(resids),len(donors_h_pairs),2),dtype=bool)
		alles = np.array([[donors_side.names==i for i in zip(*donors_h_pairs)[j]] for j in range(2)])
		#---lookups for different atom names (fast because there are relatively few atom names)
		subsels = [[np.where(i)[0] for i in alles[j]] for j in range(2)]

		#---loop over heavy/light types
		for anum in range(2):
			#---loop over possible pairs
			for pnum in range(len(alles[anum])):
				#---crucial conversion back to zero-numbering from resids here
				both[donors_side.resids[subsels[anum][pnum]]-1,pnum,anum] = True

		#---use all to find out which residues have which opportunities for bonding
		bond_opps = np.transpose(np.where(np.all(both,axis=2)))

		#---some hydrogen bonds have the same donors for multiple hydrogens
		donors_inds = np.zeros((2,len(bond_opps))).astype(int)
		for anum in range(2):
			donors_names_u = np.unique(zip(*donors_h_pairs)[anum])
			#---for each bond opportunity, we list the heavy donor
			donors_side_names = np.array(donors_h_pairs).T[anum][bond_opps[:,1]]
			#---convert this into index (this is fast because it is over a short list of donor names)
			donors_side_inds = -1*np.ones(len(donors_side_names)).astype(int)
			for nn,n in enumerate(donors_names_u): donors_side_inds[np.where(donors_side_names==n)] = nn
			#---lookup from residue and unique heavy donor atom to absolute index in donors_side
			lookup = len(donors_names_u)*np.ones((len(donors_resids),len(donors_names_u)+1)).astype(int)
			#---convert this into index (this is fast because it is over a short list of donor names)		
			donors_side_names_inds = len(donors_names_u)*np.ones((len(donors_side.names))).astype(int)
			for nn,n in enumerate(donors_names_u): 
				donors_side_names_inds[np.where(donors_side.names==n)] = nn		
			lookup[tuple(np.transpose([donors_side.resids-1,donors_side_names_inds]).T)] = \
				np.arange(len(donors_side.resids))
			#---translate bond_opps from pair numbering to heavy donor numbering (which is unique)
			bond_opps_unique = np.array(bond_opps)
			bond_opps_unique[:,1] = donors_side_inds
			donors_inds[anum] = lookup[tuple(bond_opps_unique.T)]

	#---prepare the acceptors selections
	acceptors_names = np.unique([j for k in [i.get('acceptors',[]) 
		for i in hydrogen_bond_ref.values()] for j in k])
	#---! use acceptor_names not acceptors_names
	acceptors_side = uni.select_atoms(' or '.join(['name %s'%i for i in acceptor_names]))

	#---extend to include salt bridges
	#---some systems have two types of cations
	cation_names = work.meta[sn].get('cations',work.meta[sn].get('cation',None))
	if not cation_names: raise Exception('add "cations" to the meta dictionary for %s'%sn)
	if type(cation_names)!=list: cation_names = [cation_names]
	multiple_cations = len(cation_names)>1
	cations_side = uni.select_atoms(' or '.join(['name %s'%i for i in cation_names]))

	#---prepare coordinates for each frame
	st = time.time()
	vecs,all_donor_coords,all_acceptor_coords,all_h_coords,all_cation_coords = [],[],[],[],[]
	#---purposefully profligate with the memory so this goes quickly
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
		uni.trajectory[fr]
		vecs.append(uni.dimensions[:3]/lenscale)
		all_donor_coords.append(donors_side.positions/lenscale)
		all_h_coords.append(donors_side.positions/lenscale)
		all_acceptor_coords.append(acceptors_side.positions/lenscale)
		all_cation_coords.append(cations_side.positions/lenscale)
	status('completed caching in %.1f minutes'%((time.time()-st)/60.),tag='status')
	#---the preceding code is identical to the beginning of hydrogen_bonding

	#---export variables
	from codes import hbonds
	hbonds.hydrogen_bond_ref = hydrogen_bond_ref
	hbonds.all_donor_coords = all_donor_coords
	hbonds.all_acceptor_coords = all_acceptor_coords
	hbonds.all_h_coords = all_h_coords
	hbonds.all_cation_coords = all_cation_coords
	hbonds.vecs = vecs

	#---extra exports for development
	hbonds.acceptors_resids = acceptors_side.resids
	hbonds.acceptors_resnames = acceptors_side.resnames
	hbonds.donors_resids = donors_side.resids
	hbonds.donors_resnames = donors_side.resnames

	#---debug
	if debug:
		fr = 36
		incoming_salt = hbonds.salt_bridges_framewise(
			fr,distance_cutoff=distance_cutoff)
		import ipdb;ipdb.set_trace()
		sys.quit()

	start = time.time()
	out_args = {'distance_cutoff':distance_cutoff}
	if run_parallel:
		incoming_salt = Parallel(n_jobs=4,verbose=10 if debug else 0)(
			delayed(hbonds.salt_bridges_framewise,has_shareable_memory)(fr,**out_args) 
			for fr in framelooper(nframes,start=start))
	else: 
		incoming,incoming_salt = [],[]
		for fr in framelooper(nframes):
			incoming_salt.append(hbonds.salt_bridges_framewise(fr,**out_args))

	#---extension to salt bridges. tabulate each salt
	valid_frames_salt = np.array([ii for ii,i in enumerate(incoming_salt) if len(i)>0])
	obs_by_frames_salt = np.array([len(i) for ii,i in enumerate(incoming_salt) if len(i)>0]).astype(int)
	#---some simulations have no salt bridges
	if len(valid_frames_salt)==0: bonds_salt,counts_per_frame_salt = np.array([]),np.array([])
	else:
		salt_cat = np.concatenate([incoming_salt[i] for i in valid_frames_salt])
		status('tabulating all distinct salt bridges',tag='compute')
		status('sluggish sequence because there are {:,} bond observations'.format(len(salt_cat)),tag='warning')

		#---having excluded water we are small enough to add back the ion name
		tabulation_salt = np.transpose((
			acceptors_side.resnames[salt_cat[:,0]],
				acceptors_side.resids[salt_cat[:,0]],
				acceptors_side.names[salt_cat[:,0]],
			donors_side.resnames[salt_cat[:,2]],
				donors_side.resids[salt_cat[:,2]],
				donors_side.names[salt_cat[:,2]],
			#cations_side[salt_cat[:,1]].resids,
			#cations_side[salt_cat[:,1]].resnames,
			))
		#---send the hydrogen bonds to the tabulator
		tabulation_salt_out = tabulation_salt
		bonds_salt,counts_per_frame_salt = tabulator(tabulation_salt_out,
			valid_frames_salt,obs_by_frames_salt)

		#---!!! development note. data are too big so we discard cation data in the tabulation
		#---!!! ...which means that 

	#---package the dataset
	result,attrs = {},{}
	#---everything is indexed by idx
	result['bonds'] = bonds_salt
	result['observations'] = counts_per_frame_salt
	result['bonds_salt'] = bonds_salt
	result['valid_frames'] = valid_frames_salt
	result['nframes'] = np.array(nframes)
	result['resnames'] = resnames
	result['nmols'] = rescounts
	status('compute job lasted %.1fmin'%((time.time()-start_job_time)/60.),tag='time')
	return result,attrs

