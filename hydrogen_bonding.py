#!/usr/bin/env python

import os,glob,re,time
import numpy as np
import sklearn
import sklearn.neighbors

vecnorm = lambda vec: vec/np.linalg.norm(vec)
vecangle = lambda v1,v2 : np.arccos(np.dot(vecnorm(v1),vecnorm(v2)))*(180./np.pi)

def hydrogen_bonding(grofile,trajfile,**kwargs):

	"""
	Generic hydrogen bonding code.
	Revamped on 2017.4.28 to generate a more uniform data structure.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)

	#---settings
	distance_cutoff,angle_cutoff = [calc['specs'][i] for i in ['distance_cutoff','angle_cutoff']]
	#---cutoff for inferring hydrogens from a one-time distance search
	distance_h_cutoff = distance_cutoff

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
		acceptor_names = [i['atom'] for i in mol['atoms'] if re.match('^(N|O|S)',i['atom'])]
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
	hydrogens = uni.select_atoms(' or '.join(['name %s'%i for i in hydrogens_names]))#+' or name NA')

	donors_h_pairs = [m for n in [i.get('donors',[]) for i in hydrogen_bond_ref.values()] for m in n]
	donors_h_pairs_flat = list(set([i for j in donors_h_pairs for i in j]))
	sel_d,sel_h = [uni.select_atoms(' or '.join(['name %s'%i 
		for i in list(set(zip(*donors_h_pairs)[j]))])) for j in range(2)]
	resids = np.unique(np.concatenate([sel_d.resids,sel_h.resids]))
	donors_side = uni.select_atoms(' or '.join(['name %s'%i for i in donors_h_pairs_flat]))
	donors_resids = np.unique(donors_side.resids)

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
	acceptors_names = np.unique([j for k in 
		[i.get('acceptors',[]) for i in hydrogen_bond_ref.values()] for j in k])
	acceptors_side = uni.select_atoms(' or '.join(['name %s'%i for i in acceptors_names]))

	#---prepare coordinates for each frame
	st = time.time()
	vecs,all_donor_coords,all_acceptor_coords,all_h_coords = [],[],[],[]
	#---purposefully profligate with the memory so this goes quickly
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
		uni.trajectory[fr]
		vecs.append(uni.dimensions[:3]/lenscale)
		all_donor_coords.append(donors_side.positions[donors_inds[0]]/lenscale)
		all_h_coords.append(donors_side.positions[donors_inds[1]]/lenscale)
		all_acceptor_coords.append(acceptors_side.positions/lenscale)
	status('completed caching in %.1f minutes'%((time.time()-st)/60.),tag='status')

	#---export variables
	from codes import hbonds
	hbonds.hydrogen_bond_ref = hydrogen_bond_ref
	hbonds.all_donor_coords = all_donor_coords
	hbonds.all_acceptor_coords = all_acceptor_coords
	hbonds.all_h_coords = all_h_coords
	hbonds.vecs = vecs
		
	#---debug
	if debug:
		fr = 36
		incoming = hbonds.hbonder_framewise(fr,distance_cutoff=distance_cutoff,angle_cutoff=angle_cutoff)
		sys.quit()

	start = time.time()
	out_args = {'distance_cutoff':distance_cutoff,'angle_cutoff':angle_cutoff}
	if run_parallel:
		incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
			delayed(hbonds.hbonder_framewise,has_shareable_memory)(fr,**out_args) 
			for fr in framelooper(nframes,start=start))
	else: 
		incoming = []
		for fr in framelooper(nframes):
			incoming.append(hbonds.hbonder_framewise(fr,**out_args))

	#---get valid frames
	valid_frames = np.where([len(i['donors'])>0 for i in incoming])[0]
	#---concatenate the donor/acceptor indices across all frames
	donor_cat = np.concatenate([donors_inds[0][incoming[i]['donors']] for i in valid_frames]).astype(int)
	acceptor_cat = np.concatenate([incoming[i]['acceptors'] for i in valid_frames]).astype(int)
	obs_by_frames = np.array([len(incoming[i]['acceptors']) for i in valid_frames]).astype(int)

	#---compile the massive bond table
	if False:
		tabulation = []
		for ff,fr in enumerate(valid_frames):
			status('compiling large bond table',i=ff,looplen=len(valid_frames),tag='compute')
			iii = incoming[fr]
			donors_resnames = donors_side.resnames[donors_inds[0][iii['donors']]]
			donors_resids = donors_side.resids[donors_inds[0][iii['donors']]]
			donors_names1 = donors_side.names[donors_inds[0][iii['donors']]]
			donors_names2 = donors_side.names[donors_inds[1][iii['donors']]]
			acceptors_resnames = acceptors_side.resnames[iii['acceptors']]
			acceptors_resids = acceptors_side.resids[iii['acceptors']]
			acceptors_names = acceptors_side.names[iii['acceptors']]
			tabulation.append([
				donors_resnames,donors_resids,donors_names1,donors_names2,
				acceptors_resnames,acceptors_resids,acceptors_names])

	start_time = time.time()
	#---tabulate each bond observation
	status('sluggish sequence because there are {:,} bond observations'.format(len(donor_cat)),tag='warning')
	status('tabulating all distinct hydrogen bonds',tag='compute')
	tabulation = np.transpose((donors_side.resnames[donor_cat],donors_side.resids[donor_cat],
		donors_side.names[donor_cat],acceptors_side.resnames[acceptor_cat],
		acceptors_side.resids[acceptor_cat],acceptors_side.names[acceptor_cat],))
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')

	#---reduce tabulation by discarding all SOL-SOL bonds
	#---...note that this is necessary because we have 33M observations and almost all of them are "unique"
	#---...precisely because so many of them involve water
	#---actually, instead of discarding, let us change all waters to a single residue
	tabulation_explicit = tabulation
	tabulation = np.array(tabulation_explicit)
	for p in [0,3]:
		sols = np.where(tabulation[:,p]=='SOL')[0]
		tabulation[(sols,(np.ones((len(sols)))*(p+1)).astype(int))] = '1'

	start_time = time.time()
	status('unique-ifying the tabulated bonds (estimated %ds)'%(len(donor_cat)*1.3*10**-6),tag='compute')
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
	counts_per_frame = np.zeros((len(valid_frames),len(idx)))
	#---hash the binds over the indices
	bonds_to_idx = dict([(tuple(b),bb) for bb,b in enumerate(bonds)])
	frame_lims = np.concatenate(([0],np.cumsum(obs_by_frames)))
	for fr,i in enumerate(frame_lims[:-1]):
		status('counting observations per frame',i=fr,looplen=len(valid_frames),
			tag='compute',start=start_time)
		obs = tabulation[frame_lims[fr]:frame_lims[fr+1]]
		counts_per_frame[fr][np.array([bonds_to_idx[tuple(o)] for o in obs])] += 1
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')
	status('done heavy lifting',tag='compute')

	#---package the dataset
	result,attrs = {},{}
	#---everything is indexed by idx
	result['bonds'] = bonds
	result['observations'] = counts_per_frame
	result['valid_frames'] = valid_frames
	result['nframes'] = np.array(nframes)
	result['resnames'] = resnames
	result['nmols'] = rescounts
	status('compute job lasted %.1fmin'%((time.time()-start_job_time)/60.),tag='time')
	return result,attrs
