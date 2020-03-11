#!/usr/bin/env python

import sys,time,re
import numpy as np
import scipy
import scipy.spatial
import MDAnalysis
from codes.looptools import basic_compute_loop

def get_lipid_resnames():
	"""
	Brief utility for getting the lipid names from automacs.
	"""
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
	return targets

def torusnorm(pts1,pts2,vecs,paired=False):
	"""
	Compute distances between points on a torus. This is the paired version in contrast to the combination verison.
	"""
	#---adapted from torusnorm in codes/mesh.py and switched from 2 dimensions to 3. so it's a klein bottle (std PBCs)?
	if not paired: cd = np.array([scipy.spatial.distance.cdist(pts1[:,d:d+1],pts2[:,d:d+1]) for d in range(3)])
	else: cd = np.abs(pts1-pts2)
	for d in range(3): cd[d] -= (cd[d]>vecs[d]/2.)*vecs[d]
	cd2 = np.linalg.norm(cd,axis=1 if paired else 0)
	return cd2

def hydration_distribution_framewise(fr,knn):
	"""Compute nearest lipids and nearest waters."""
	global vecs,subject_coords,proxy_coords,water_coords
	vec = vecs[fr]
	pts_back_unstuffed = proxy_coords[fr]
	pts_fore_unstuffed = subject_coords[fr]
	#---standard distance method follows
	#---ensure that the points are inside the box
	boxstuff = lambda pts,vec : pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec
	pts_back = boxstuff(pts_back_unstuffed,vec)
	pts_fore = boxstuff(pts_fore_unstuffed,vec)
	#---! why does vec need to be twice as long? (tested that the limits work though)
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return [None,None,None]
	close,nns = tree.query(pts_fore,k=1)
	#---compute distances with PBC correction
	### torusnorm is broken!!! ld = torusnorm(proxy_coords[fr][nns],subject_coords[fr],vecs[fr],paired=True)
	ld_raw = np.abs(proxy_coords[fr][nns]-subject_coords[fr])
	ld = np.linalg.norm(ld_raw-(ld_raw>vec/2.0)*vec,axis=1)
	#---find the nearby waters using the distance method again with the same box stuff
	pts_back_unstuffed = water_coords[fr]
	pts_back = boxstuff(pts_back_unstuffed,vec)
	#---! why does vec need to be twice as long? (tested that the limits work though)
	try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
	#---KDTree failures are blanked
	except: return [None,None,None]
	close,nns = tree.query(pts_fore,k=knn)
	#---collect the water distances, in dimensions: ion number by proximity rank
	wd_raw = np.abs((np.tile(subject_coords[fr],(knn,1,1)).transpose((1,0,2))-water_coords[fr][nns]))
	#---numpy norm has an extremely useful axis keywrod now so the distance calculation is extremely elegant
	wd = np.linalg.norm(wd_raw-(wd_raw>vec/2.)*vec,axis=2)
	return wd,ld,fr

def hydration_distribution(grofile,trajfile,**kwargs):

	"""
	Compute the radial distribution function (RDF) a.k.a g(r) of water around ions but filter these distributions
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	debug = kwargs.get('debug',False)
	run_parallel = kwargs.get('run_parallel',True)
	start_job_time = time.time()
	#---nearest water distances to calculate
	knn = calc['specs'].get('k_nearest_waters',200)

	#---prepare universe	
	print((grofile,trajfile))
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	lenscale = 10.

	#---collect selection strings
	lipid_resnames = get_lipid_resnames()
	cation_names = work.meta[sn].get('cations',work.meta[sn].get('cation',None))
	if type(cation_names)!=list: cation_names = [cation_names]

	#---define selections
	sel_proxy = uni.select_atoms(' or '.join(['resname %s'%i for i in lipid_resnames]))
	sel_subject = uni.select_atoms(' or '.join(['name %s'%i for i in cation_names]))
	#---we use oxygen to denote the water
	sel_water = uni.select_atoms('resname SOL and name OW')

	#---prepare coordinates for each frame
	st = time.time()
	global vecs,subject_coords,proxy_coords,water_coords
	vecs,subject_coords,proxy_coords,water_coords = [],[],[],[]
	#---purposefully profligate with the memory so this goes quickly
	for fr in range(nframes):
		status('caching coordinates',tag='compute',i=fr,looplen=nframes,start=st)	
		uni.trajectory[fr]
		vecs.append(uni.dimensions[:3]/lenscale)
		subject_coords.append(sel_subject.positions/lenscale)
		proxy_coords.append(sel_proxy.positions/lenscale)
		water_coords.append(sel_water.positions/lenscale)
	status('completed caching in %.1f minutes'%((time.time()-st)/60.),tag='status')

	#---convert back to advanced indexing
	aind = lambda x : tuple(x.T)

	water_distances,lipid_distances,valid_frames = [],[],[]
	#---loop over frames
	st = time.time()
	looper = [dict(fr=fr,knn=knn) for fr in range(nframes)]
	incoming = basic_compute_loop(hydration_distribution_framewise,looper=looper)
	water_distances,lipid_distances,valid_frames = zip(*incoming)
	valid_frames = [fr for fr in valid_frames if fr!=None]
	water_distances = [water_distances[fr] for fr in valid_frames]
	lipid_distances = [lipid_distances[fr] for fr in valid_frames]

	#---package the dataset
	result,attrs = {},{}
	result['water_distances'] = np.array(water_distances)
	result['lipid_distances'] = np.array(lipid_distances)
	result['valid_frames'] = valid_frames
	result['nframes'] = np.array(nframes)
	result['cation_resids'] = sel_subject.resids
	return result,attrs

