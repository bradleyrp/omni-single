#!/usr/bin/env python

"""
Contact code derived from hydrogen_bonding and/or salt_bridges.
Note that the hydrogen bonding and salt bridges code explicitly codes for three-party interactions of 
a highly specific nature, whereas this code generalizes this to nearby contacts between two elements
with (1) the option to specify that the subject is e.g. the protein residues (instead of checking all 
contacts in the system, which would be costly) and (2) the same exact data structure as the two progenitor
codes so that later these contacts can be easily summed.
"""

import os,sys,glob,re,time
import numpy as np
import scipy
import scipy.spatial

str_types = [str,unicode] if sys.version_info<(3,0) else [str]

vecnorm = lambda vec: vec/np.linalg.norm(vec)
vecangle = lambda v1,v2 : np.arccos(np.dot(vecnorm(v1),vecnorm(v2)))*(180./np.pi)

def contacts_framewise(fr,**kwargs):

	"""
	"""

	global vecs,coords_targ,coords_subj
	distance_cutoff = kwargs['distance_cutoff']
	lenscale = kwargs.get('lenscale',10.0)
	vec = vecs[fr]

	#---ensure that the points are inside the box
	boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
	#---convert back to advanced indexing
	aind = lambda x : tuple(x.T)

	pts_back_unstuffed = coords_targ[fr]
	pts_fore_unstuffed = coords_subj[fr]
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)

	#---! why does vec need to be twice as long? (tested that the limits work though)
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return {'subjects':np.array([]),'targets':np.array([])}
	close,nns = tree.query(pts_fore,k=10,distance_upper_bound=distance_cutoff/lenscale)

	#---index pairs within the cutoff distance
	close_pairs = np.transpose(np.where(np.all((close<distance_cutoff/lenscale,close>0),axis=0)))
	#---list of close donors
	close_targets = nns[aind(close_pairs)]
	close_subjects = close_pairs[:,0]
	return {'subjects':close_subjects,'targets':close_targets}

def contacts(grofile,trajfile,**kwargs):

	"""
	GENERIC CONTACT CODE.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)

	#---settings
	lenscale = 10.0
	#---distance cutoff stays in angstroms until the compute function
	distance_cutoff = calc['specs']['cutoff']
	subject_selection = calc['specs'].get('subject','protein')
	object_flag = calc['specs'].get('object','lipid')

	#---prepare universe	
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
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

	#---get the subject of the calculation, the thing we wish to study the contacts of
	#---...typically the protein
	#---! need to add resid redundancy checks possibly
	subject = uni.select_atoms(subject_selection)
	#---get the objects
	if object_flag=='lipid':
		#---objects from the landscape returns resnames
		target_resnames = land.objects_by_category('lipid')
		#---explicitly ignore hydrogen contacts here
		targets = uni.select_atoms('(%s) and not name H*'%
			' or '.join(['resname %s'%i for i in target_resnames]))
	else: raise Exception('not set up for object %s'%object_flag)

	#---prepare coordinates for each frame
	st = time.time()
	global vecs,coords_subj,coords_targ
	vecs,coords_subj,coords_targ,times = [],[],[],[]
	#---purposefully profligate with the memory so this goes quickly
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
		uni.trajectory[fr]
		times.append(uni.trajectory.time)
		vecs.append(uni.dimensions[:3]/lenscale)
		coords_subj.append(subject.positions/lenscale)
		coords_targ.append(targets.positions/lenscale)
	status('completed caching in %.1f minutes'%((time.time()-st)/60.),tag='status')

	#---debug
	compute_function = contacts_framewise
	if debug:
		fr = 800
		incoming = compute_function(fr,distance_cutoff=distance_cutoff)
		sys.quit()

	#---compute loop
	start = time.time()
	out_args = {'distance_cutoff':distance_cutoff}
	if run_parallel:
		incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
			delayed(compute_function,has_shareable_memory)(fr,**out_args) 
			for fr in framelooper(nframes,start=start))
	else: 
		incoming = []
		for fr in framelooper(nframes):
			incoming.append(compute_function(fr,**out_args))

	#---chompdown
	#---get valid frames
	valid_frames = np.where([len(i['subjects'])>0 for i in incoming])[0]
	obs_by_frames = np.array([len(incoming[i]['subjects']) for i in valid_frames]).astype(int)
	#---concatenate the donor/acceptor indices across all frames
	subject_cat = np.concatenate([incoming[i]['subjects'] for i in valid_frames]).astype(int)
	target_cat = np.concatenate([incoming[i]['targets'] for i in valid_frames]).astype(int)

	start_time = time.time()
	#---tabulate each bond observation
	tabulation = np.transpose((subject.resnames[subject_cat],subject.resids[subject_cat],
		subject.names[subject_cat],targets.resnames[target_cat],targets.resids[target_cat],
		targets.names[target_cat],))
	status('stopwatch: %.1fs'%(time.time()-start_time),tag='compute')

	#---! move this somewhere more centrally located
	from art_ptdins import uniquify
	idx,counts = uniquify(tabulation.astype(str))
	bonds = tabulation[idx]

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
	result['resnames'] = resnames_master
	result['subject_residues_resnames'] = subject.residues.resnames
	result['subject_residues_resids'] = subject.residues.resids
	result['nmols'] = rescounts
	result['times'] = np.array(times)
	status('compute job lasted %.1fmin'%((time.time()-start_job_time)/60.),tag='time')
	return result,attrs
