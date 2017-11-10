#!/usr/bin/env python

import os,sys,glob,re,time
import numpy as np
import sklearn
import sklearn.neighbors

str_types = [str,unicode] if sys.version_info<(3,0) else [str]

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

	#---prototyping an external module import to replace the sometimes tedious addition of 
	#---... many metadata to the meta dictionary in your YAML files.
	protein_itp_loader = work.vars.get('protein_itp_loader',None)
	#---default ITP paths are set for each simulation in the metadata
	if protein_itp_loader==None: 
		def protein_itp_loader(sn,**kwargs):
			itp = work.meta.get(sn,{}).get('protein_itp',None)
			if not itp: raise Exception('cannot find protein_itp in meta dictionary for %s '%sn+
				'note that you can also use the protein_itp_loader functionality to get the ITP file')
	#---custom ITP loader specified as an alternate_module
	else: protein_itp_loader = alternate_module(**protein_itp_loader)

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
	resnames_master = np.array(resnames)
	rescounts = counts[np.argsort(idx)]

	import makeface
	#---get an automacs landscape with a little help from the user
	try: mod = makeface.import_remote('amx/amx')
	except: raise Exception('please clone a copy of automacs next to omni in `amx`. '
		'you must also run `make setup all` from that directory to get force field files.')
	mod['state'].force_field = 'charmm'
	Landscape = mod['Landscape']
	land = Landscape(cwd='amx/')
	#---use the landscape to get hydrogen bond donors and acceptors for lipids
	hydrogen_bond_ref = {}
	#---this section relies on correct definitions from the Landscape
	targets = land.objects_by_category('lipid')+land.objects_by_category('sterol')
	#---!!!!!!!!!!!!!!!!!!!!!!!!!!! need to add all of the lipids here and in the charmm36.ff folder!
	#---METHODOLOGY NOTE: we catalog all hydrogen bonding opportunities ONLY BY NAME
	#---loop over lipid targets and scan them for hydrogen bond opportunities
	for resname in targets:
		#---each lipid ITP has an identical molecule with the same (residue) name
		itp = mod['GMXTopology'](land.objects[resname]['fn'])
		#---donor names come from a double-regex match over bonds
		donor_names = itp.get_bonds_by_regex(molname=resname,patterns=['^H','^(N|O|S)'])
		#---acceptor names have a single regex
		#---!!! check that this is the correct definition
		acceptor_names = [i['atom'] for i in 
			itp.molecules[resname]['atoms'] if re.match('^(N|O|S)',i['atom'])]
		hydrogen_bond_ref[resname] = {'acceptors':acceptor_names,'donors':donor_names}
	#---include any proteins as participants in the bonding
	if kwargs['calc']['specs'].get('protein',False):
		#---get the protein ITP from metadata
		itp_fn = work.meta[sn].get('protein_itp',protein_itp_loader(sn,work=work))
		if not itp_fn: raise Exception('add protein_itp to the meta for %s'%sn)
		#---get the sims spot path systematically
		if os.path.isfile(itp_fn): itp_fn_abs = itp_fn
		#---if path is relative then we consult the spots
		else:
			rootdir = work.raw.spots[(work.raw.spotname_lookup(sn),'structure')]['rootdir']
			sn_dir = os.path.join(rootdir,sn)
			#---user supplies step folder and path to the reference structure
			itp_fn_abs = os.path.join(sn_dir,itp_fn)
		protein_itp = mod['GMXTopology'](itp_fn_abs)
		for molname in protein_itp.molecules:
			#---mimic the procedure above for lipids
			#---donor names come from a double-regex match over bonds
			donor_resnames_names = protein_itp.get_bonds_by_regex(molname=molname,
				patterns=['^H','^(N|O|S)'],include_resname=True)
			#---organize hydrogen bonds by residue name
			resnames_all = list(set([i for j in zip(*donor_resnames_names)[0] for i in j]))
			for resname_focus in resnames_all:
				donor_list = []
				#---loop over resnames within the protein
				for resnames,names in donor_resnames_names:
					if resnames[0]!=resnames[1]: 
						raise Exception('invalid hydrogen bond spec %s,%s'%(resnames,names))
					elif resnames[0]==resname_focus: donor_list.append(names)
					else: continue
				#---acceptor names have a single regex
				#---!!! check that this is the correct definition
				acceptor_names = list(set([i['atom'] for i in 
					protein_itp.molecules[molname]['atoms'] if re.match('^(N|O|S)',i['atom'])
					and i['resname']==resname_focus]))
				hydrogen_bond_ref[(molname,resname_focus)] = {
					'acceptors':acceptor_names,'donors':donor_list}

	"""
	developing a new method for selecting the atoms correctly
	we need to get all possible donors into a big selection after which case the hbonds.hbonder_framewise
		does the rest and the plotting codes are decent at picking out inter-residue bonds and identities 
	consider the customer: hbonds_framewise 
		needs a list of donors, hydrogens, and acceptors
		since the donors and hydrogens must be related by bonds
			there is some redundancy in the donor list
			which redundancy hbonds_framewise accounts for
	construct the donors list
		start with the list of all hydrogens (somewhat counterintuitive)
		consult the hydrogen_bond_ref and loop over all residues and then grow a list of hydrogen indices
		for each hydrogen find the associated heavy atom and and add both indices to separate lists
		net result is two lists of indices over the hydrogens and donors which constitute a bond
	"""

	#---get the heavy atom side of the donors
	donors_heavy,donors_h = [uni.select_atoms(' or '.join(['(resname %s and (%s))'%(
		resname if type(resname) in str_types else resname[1],' or '.join(['name %s'%i[w] 
		for i in v['donors']])) for resname,v in hydrogen_bond_ref.items() 
		if v['donors']])) for w in range(2)]
	acceptors_heavy = uni.select_atoms(' or '.join(['(resname %s and (%s))'%(
		resname if type(resname) in str_types else resname[1],' or '.join(['name %s'%i 
		for i in v['acceptors']])) for resname,v in hydrogen_bond_ref.items() 
		if v['acceptors']]))

	#---check non-redundany residues
	if not len(donors_heavy.residues)==len(np.unique(donors_heavy.resids)):
		raise Exception('residue redundancy in the donor heavy list')
	if not len(donors_h.residues)==len(np.unique(donors_h.resids)):
		raise Exception('residue redundancy in the donor hydrogen list')
	#---constructing the donors side selection to preserve the bond relation
	donors_reindex = []
	for refkey,details in hydrogen_bond_ref.items():
		#---protein residues have the protein molecule name alongside
		resname = refkey if type(refkey) in str_types else refkey[1]
		for heavy,light in details['donors']:
			inds_heavy = np.where(np.all((
				donors_heavy.resnames==resname,donors_heavy.names==heavy),axis=0))[0]
			inds_light = np.where(np.all((
				donors_h.resnames==resname,donors_h.names==light),axis=0))[0]
			#---loop over resids and for each resid that has them, we add the indices to the list
			#---! these descending loops are clumsy but they should be fast and they make definitions precise
			for resid in np.unique(np.concatenate((donors_heavy[inds_heavy].resids,
				donors_h[inds_light].resids))):
				inds_heavy = np.where(np.all((
					donors_heavy.resnames==resname,donors_heavy.names==heavy,donors_heavy.resnums==resid
					),axis=0))[0]
				inds_light = np.where(np.all((
					donors_h.resnames==resname,donors_h.names==light,donors_h.resnums==resid),axis=0))[0]
				if len(inds_heavy)>1 or len(inds_light)>1: 
					raise Exception('serious error! one unique hydrogen bond in a single residue')
				if len(inds_heavy)==1 and len(inds_light)==1:
					donors_reindex.append((inds_heavy[0],inds_light[0]))
	#---the reindexed donors preserved the bond relation and covers all possible unique hydrogen bonds
	donors_reindex = np.array(donors_reindex)

	#---prepare coordinates for each frame
	st = time.time()
	vecs,all_donor_coords,all_acceptor_coords,all_h_coords = [],[],[],[]
	#---purposefully profligate with the memory so this goes quickly
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
		uni.trajectory[fr]
		vecs.append(uni.dimensions[:3]/lenscale)
		all_donor_coords.append(donors_heavy.positions[donors_reindex[:,0]]/lenscale)
		all_h_coords.append(donors_h.positions[donors_reindex[:,1]]/lenscale)
		all_acceptor_coords.append(acceptors_heavy.positions/lenscale)
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
		hbonds.donors_side = donors_side
		hbonds.donors_inds = donors_inds
		hbonds.donors_inds = donors_inds
		hbonds.acceptors_side = acceptors_side
		fr = 686 #---careful debugging at this frame
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
	donor_cat,donor_cat_h = [np.concatenate([donors_reindex.T[j][incoming[i]['donors']] 
		for i in valid_frames]).astype(int) for j in range(2)]
	acceptor_cat = np.concatenate([incoming[i]['acceptors'] for i in valid_frames]).astype(int)
	obs_by_frames = np.array([len(incoming[i]['acceptors']) for i in valid_frames]).astype(int)

	start_time = time.time()
	#---tabulate each bond observation
	status('sluggish sequence because there are {:,} bond observations'.format(len(donor_cat)),tag='warning')
	status('tabulating all distinct hydrogen bonds',tag='compute')
	tabulation = np.transpose((donors_heavy.resnames[donor_cat],donors_heavy.resids[donor_cat],
		donors_heavy.names[donor_cat],acceptors_heavy.resnames[acceptor_cat],
		acceptors_heavy.resids[acceptor_cat],acceptors_heavy.names[acceptor_cat],
		#---include the hydrogen identity here in the tabulation (note this might make things larger?)
		#---also note that the hydrogen atom name should be enough because we already have the donor resid
		donors_h.names[donor_cat_h],))
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
		status('counting observations',i=fr,looplen=len(valid_frames),
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
	result['resnames'] = resnames_master
	result['nmols'] = rescounts
	status('compute job lasted %.1fmin'%((time.time()-start_job_time)/60.),tag='time')
	return result,attrs
