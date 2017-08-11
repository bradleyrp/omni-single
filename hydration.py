#!/usr/bin/env python

import sys,time,re
import numpy as np
import scipy
import scipy.spatial
from codes.looptools import basic_compute_loop

#---ensure that the points are inside the box
boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
#---convert back to advanced indexing
aind = lambda x : tuple(x.T)
#---globals for shared memory
pts_ions,pts_water,pts_lipids,vecs,midplanes = None,None,None,None,None

def shell_counter(fr,cutoff):
	"""
	Compute objects near the subjects within the cutoff.
	"""
	global pts_ions,pts_water,vecs
	vec = vecs[fr]
	pts_back_unstuffed = pts_water[fr]
	pts_fore_unstuffed = pts_ions[fr]
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)
	#---! why does vec need to be twice as long? (tested that the limits work though)
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return np.array([])
	close,nns = tree.query(pts_fore,k=10,distance_upper_bound=cutoff)
	#---now it is easy to count the nearby waters
	return (close<=cutoff).sum(axis=1)

def minimum_distances(fr,**kwargs):
	"""
	Return the nearest object to each subject.
	"""
	global pts_ions,pts_lipids,vecs,midplanes
	distance_metric = kwargs.pop('distance_metric','r')
	if kwargs: raise Exception('unprocessed kwargs')
	vec = vecs[fr]
	pts_fore_unstuffed = pts_ions[fr]
	pts_fore = boxstuff(pts_fore_unstuffed,vec)
	if distance_metric=='r':
		pts_back_unstuffed = pts_lipids[fr]
		pts_back = boxstuff(pts_back_unstuffed,vec)
		#---! why does vec need to be twice as long? (tested that the limits work though)
		try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
		#---KDTree failures are blanked
		except: return np.array([])
		close,nns = tree.query(pts_fore,k=1)
	elif distance_metric=='z':
		#---no PBCs implemented here and we really just take the distance to the average-z
		close = np.abs(pts_fore[:,2]-midplanes[fr])
	else: raise Exception('unclear distance metric: %s'%distance_metric)
	return close

def hydration(grofile,trajfile,**kwargs):

	"""
	Hydration code revamped from simuluxe on 2017.6.21.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)
	start_job_time = time.time()

	#---prepare universe	
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	lenscale = 10.

	#---get selections
	if ',' in work.meta[sn]['cation']: cation = work.meta[sn]['cation_relevant']
	else: cation = work.meta[sn]['cation']
	sel_ions = uni.select_atoms('name %s'%cation)
	sel_lipids_str = ' or '.join(['resname %s'%i for i in work.vars['selectors']['resnames_lipid']])
	sel_lipids = uni.select_atoms(sel_lipids_str)
	#---atom subselection
	atom_filter = calc['specs'].get('atom_filter',None)
	if atom_filter:
		sel_lipids = uni.select_atoms('(%s) and (%s)'%(sel_lipids_str,' or '.join(['name %s'%i for i in 
			np.unique([n for n in sel_lipids.names if re.match(atom_filter,n)])])))
	#---handle the distance metric
	distance_metric = calc['specs'].get('distance_metric',None)
	#---pass the distance metric to the distance finder
	distance_args = {'distance_metric':distance_metric}
	#---we use the water oxygen only
	sel_water = uni.select_atoms('name OW')
	global pts_ions,pts_water,pts_lipids,vecs,midplanes
	#---the height distance metric needs the average z for each frame
	if distance_metric=='z': 
		midplanes = np.array([i.mean() for i in kwargs['upstream']['undulations']['mesh'].mean(axis=0)])
	#---cache the points
	pts_ions = np.zeros((nframes,len(sel_ions),3))
	pts_lipids = np.zeros((nframes,len(sel_lipids),3))
	pts_water = np.zeros((nframes,len(sel_water),3))
	vecs = np.zeros((nframes,3))
	start = time.time()
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=start)	
		uni.trajectory[fr]
		pts_ions[fr] = sel_ions.positions/lenscale
		pts_lipids[fr] = sel_lipids.positions/lenscale
		pts_water[fr] = sel_water.positions/lenscale
		vecs[fr] = uni.dimensions[:3]/lenscale

	#---prepare arguments for the compute functions
	hydration_cutoff = work.vars['hydration_cutoffs'][cation]/lenscale
	out_args = dict(cutoff=hydration_cutoff)
	if False:
		if debug:
			fr = 36
			shell_counts = shell_counter(fr,**out_args)
			near_lipids = minimum_distances(fr,**distance_args)
			import ipdb;ipdb.set_trace()
			sys.exit()

	#---compute the waters in the shell
	shell_counts = basic_compute_loop(
		compute_function=shell_counter,
		looper=[dict(fr=fr,**out_args) for fr in range(nframes)],
		run_parallel=run_parallel,debug=None)
	#---select valid frames
	valid_frames_shell_counts = np.array([i for i in range(len(shell_counts)) if len(shell_counts[i])>0])
	shell_counts = np.array(shell_counts)[valid_frames_shell_counts]
	#---for each ion get the minimum distance to any lipid
	near_lipids = basic_compute_loop(
		compute_function=minimum_distances,
		looper=[dict(fr=fr,**distance_args) for fr in range(nframes)],
		run_parallel=run_parallel,debug=None)
	#---select valid frames
	valid_frames_near_lipids = np.array([i for i in range(len(near_lipids)) if len(near_lipids[i])>0])
	near_lipids = np.array(near_lipids)[valid_frames_near_lipids]

	if False:
		#---for each ion count the waters within the shell
		start = time.time()
		if run_parallel:
			shell_counts = Parallel(n_jobs=8,verbose=10 if debug else 0)(
				delayed(shell_counter,has_shareable_memory)(fr,**out_args) 
				for fr in framelooper(nframes,start=start))
		else: 
			shell_counts = []
			for fr in framelooper(nframes):
				shell_counts.append(shell_counter(fr,**out_args))
		valid_frames_shell_counts = np.array([i for i in range(len(shell_counts)) if len(shell_counts[i])>0])
		if len(valid_frames_shell_counts)==0: 
			print('something is amiss you have no valid frames')
			import ipdb;ipdb.set_trace()
		shell_counts = np.array(shell_counts)[valid_frames_shell_counts]

	if False:
		#---! note some repetition in the debug/parallel/serial blocks in many functions
		#---for each ion get the minimum distance to any lipid
		start = time.time()
		if run_parallel:
			near_lipids = Parallel(n_jobs=8,verbose=10 if debug else 0)(
				delayed(minimum_distances,has_shareable_memory)(fr,**distance_args) 
				for fr in framelooper(nframes,start=start))
		else: 
			near_lipids = []
			for fr in framelooper(nframes):
				near_lipids.append(minimum_distances(fr,**distance_args))
		valid_frames_near_lipids = np.array([i for i in range(len(near_lipids)) if len(near_lipids[i])>0])
		near_lipids = np.array(near_lipids)[valid_frames_near_lipids]

	#---package the dataset
	result,attrs = {},{}
	#---everything is indexed by idx
	attrs['hydration_cutoff'] = hydration_cutoff
	result['nframes'] = np.array(nframes)
	result['shell_counts'] = shell_counts
	result['valid_frames_shell_counts'] = valid_frames_shell_counts
	result['near_lipids'] = near_lipids
	result['valid_frames_near_lipids'] = valid_frames_near_lipids
	status('compute job lasted %.1fmin'%((time.time()-start_job_time)/60.),tag='time')
	return result,attrs
