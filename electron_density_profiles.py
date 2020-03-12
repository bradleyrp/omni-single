#!/usr/bin/env python

import time
import numpy as np
import MDAnalysis
from joblib import Parallel,delayed
# from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
from codes.looptools import basic_compute_loop
from codes.mesh import centroid
import re

def compute_edp_single(**kwargs):
	"""Compute electron density profiles for groups."""
	fr = kwargs['fr']
	global vecs,coords,nbins,groups,midpoint,charges
	vec = vecs[fr]
	zs = coords[fr,:,2]
	half_box_height = vec[2]/2.
	h_rel = (zs-midpoint[fr])
	h_rel_wrapped = (h_rel + ((-1*(h_rel>half_box_height)+(h_rel<-1*half_box_height)) *half_box_height*2))
	if not np.all(np.abs(h_rel_wrapped)<=half_box_height): 
		raise Exception('something leaked out of the box')
	dividers = np.linspace(-1*half_box_height,half_box_height,nbins+1)
	tabulated = np.zeros((len(groups),nbins))
	for gnum,group in enumerate(groups): 
		tabulated[gnum],_ = np.histogram(h_rel_wrapped[group],bins=dividers,weights=charges[group])
	return tabulated

def electron_density_profiles(**kwargs):
	"""
	Compute the electron density profiles.
	"""
	global vecs,coords,nbins,groups,midpoint,charges

	# hardcoded settings
	chargedict = {'^N(?!A$)':7,'^C[0-9]+':6,'^CL$':17,'^H':1,'^O':8,
		'^P':15,'^Cal':18,'^MG':10,'^NA':11,'^S':16,'K':18}
	# we consider residues then the following regular expressions
	group_regexes = kwargs['calc']['specs'].get('extra_regexes',['.+','^(OW)|(HW(1|2))$','^C[0-9]+'])
	# get the reference z for each frame
	bilayer_coms = kwargs['upstream']['lipid_abstractor']['points']
	imono = kwargs['upstream']['lipid_abstractor']['monolayer_indices']
	#! assume upstream lipid_abstractor is correct and the resulting points are not broken over PBCs
	midpoint = np.array([bilayer_coms[:,imono==mn][:,:,2].mean(axis=1) for mn in range(2)]).mean(axis=0)
	# get the trajectory
	grofile,trajfile = kwargs['structure'],kwargs['trajectory']
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	# MDAnalysis uses Angstroms not nm
	lenscale = 10.
	# choose a number of bins
	bin_size = kwargs['calc']['specs']['bin_size']
	vecs_upstream = kwargs['upstream']['lipid_abstractor']['vecs']
	# round the number of bins to ensure everything is flush
	nbins = np.round(vecs_upstream[:,2].mean()/bin_size).astype(int)
	# collect coordinates
	sel = uni.select_atoms('all')
	# assign charges
	namelist = uni.atoms.names
	resnamelist = list(set(uni.atoms.resnames))
	# charge dictionary for the atoms in this particular system
	chargedict_obs = dict([(name,[chargedict[key] for key in chargedict 
		if re.match(key,name)]) for name in np.unique(namelist)])
	unclear_charges = dict([(key,val) for key,val in chargedict_obs.items() if len(val)!=1])
	if any(unclear_charges): 
		raise Exception('charges for these atoms were not specified: %s'%unclear_charges)
	chargedict_obs = dict([(key,val[0]) for key,val in chargedict_obs.items()])
	charges = np.array([chargedict_obs[n] for n in namelist])
	# identify atoms for each residue type
	groups = [np.where(uni.atoms.resnames==r)[0] for r in resnamelist]
	groups += [np.array([i for i,j in enumerate(namelist) if re.match(reg,j)]) for reg in group_regexes]
	# cache the points
	coords = np.zeros((nframes,len(sel),3))
	vecs = np.zeros((nframes,3))
	for fr in range(nframes):
		status('loading frame',tag='load',i=fr,looplen=nframes)
		uni.trajectory[fr]
		vecs[fr] = uni.trajectory[fr].dimensions[:3]/lenscale
		coords[fr] = np.array(sel.positions)/lenscale
	# make sure vectors are the same
	if not np.all(vecs==vecs_upstream): 
		raise Exception('vectors do not match upstream lipid_abstractor')
	# compute
	looper = [dict(fr=fr) for fr in range(nframes)]
	incoming = basic_compute_loop(compute_edp_single,looper,run_parallel=True)
	# pack
	results,attrs = {},{}
	attrs['group_regexes'] = group_regexes
	for gnum,group in enumerate(groups):
		results['group_%d'%gnum] = group
	results['resnames'] = resnamelist
	results['tabulated'] = np.array(incoming)
	attrs['bin_size'] = bin_size
	results['midpoint'] = midpoint
	attrs['nbins'] = nbins
	results['vecs'] = vecs
	results['charges'] = charges
	return results,attrs

def electron_density_profiles_deprecated(grofile,trajfile,**kwargs):
	"""
	Compute the electron density profile
	"""
	#from calcs.codes.mesh import identify_lipid_leaflets
	bin_size = kwargs['calc']['specs']['bin_size']
	chargedict = {'^N(?!A$)':7,'^C[0-9]+':6,'^CL$':17,'^H':1,
		'^O':8,'^P':15,'^Cal':18,'^MG':10,
		'^NA':11,'^S':16,'K':18}
	group_regexes = ['.+','^(OW)|(HW(1|2))$','^C[0-9]+']

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	
	#---prepare universe	
	grofile,trajfile = kwargs['structure'],kwargs['trajectory']
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.
	
	#---select residues of interest
	selector = kwargs['calc']['specs']['selector']
	monolayer_cutoff = kwargs['calc']['specs']['selector']['monolayer_cutoff']

	#---center of mass over residues
	if 'type' in selector and selector['type'] == 'com' and 'resnames' in selector:
		resnames = selector['resnames']
		selstring = '('+' or '.join(['resname %s'%i for i in resnames])+')'
	else: raise Exception('\n[ERROR] unclear selection %s'%str(selector))
	
	#---compute masses by atoms within the selection
	sel = uni.select_atoms(selstring)
	mass_table = {'H':1.008,'C':12.011,'O':15.999,'N':14.007,'P':30.974}
	masses = np.array([mass_table[i[0]] for i in sel.atoms.names])
	resids = sel.resids
	#---create lookup table of residue indices
	divider = [np.where(resids==r) for r in np.unique(resids)]

	#---load trajectory into memory	
	trajectory,vecs = [],[]
	for fr in range(nframes):
		status('loading frame',tag='load',i=fr,looplen=nframes)
		uni.trajectory[fr]
		trajectory.append(sel.positions/lenscale)
		vecs.append(sel.dimensions[:3])

	#---parallel
	start = time.time()
	coms = Parallel(n_jobs=work.nprocs,verbose=0)(
		delayed(centroid)(trajectory[fr],masses,divider)
		for fr in framelooper(nframes,start=start))

	#---identify monolayers 
	#---! why though?
	#---note that this could just refer to the mesh object but this is very fast
	if False: monolayer_indices = identify_lipid_leaflets(coms[0],vecs[0],
		monolayer_cutoff=monolayer_cutoff)

	#---load trajectory into memory	
	allsel = uni.select_atoms('all')
	trajectory,vecs = [],[]
	for fr in range(nframes):
		status('loading frame',tag='load',i=fr,looplen=nframes)
		uni.trajectory[fr]
		trajectory.append(allsel.positions/lenscale)
		vecs.append(allsel.dimensions[:3])
	trajectory = np.array(trajectory)
	vecs = np.array(vecs)/lenscale
		
	#---center the mean of com positions at z=0
	midplane_heights = np.array([np.mean(coms[fr],axis=0)[2] for fr in range(nframes)])
	for fr in range(nframes): trajectory[fr,:,2] -= midplane_heights[fr]

	#---correct for periodic boundaries
	for fr in range(nframes):
		trajectory[fr,:,2] -= (trajectory[fr,:,2]>vecs[fr,2]/2.)*vecs[fr,2]
		trajectory[fr,:,2] += (trajectory[fr,:,2]<-1*vecs[fr,2]/2.)*vecs[fr,2]
	#---offset by one bin width here so the use of astype(int) is symmetric later on
	trajectory[...,2] += bin_size/2.

	#---assign charges
	namelist = uni.atoms.names
	resnamelist = list(set(uni.atoms.resnames))
	#---charge dictionary for the atoms in this particular system
	chargedict_obs = dict([(name,[chargedict[key] for key in chargedict if re.match(key,name)]) 
		for name in np.unique(namelist)])
	unclear_charges = dict([(key,val) for key,val in chargedict_obs.items() if len(val)!=1])
	if any(unclear_charges): 
		raise Exception('charges for these atoms were not specified: %s'%unclear_charges)
	chargedict_obs = dict([(key,val[0]) for key,val in chargedict_obs.items()])
	charges = [chargedict_obs[n] for n in namelist]

	#---identify atoms for each residue type
	groups = [[ii for ii,i in enumerate(uni.atoms) if i.resname==r] for r in resnamelist]
	groups += [np.array([i for i,j in enumerate(namelist) if re.match(reg,j)]) for reg in group_regexes]
	
	#---bin the heights according to bin_size
	xx = np.array(np.floor(trajectory[:,:,2]/bin_size)).astype(int)
	offset = xx.min()
	xx -= xx.min()
	bincounts_by_group = [np.zeros(xx.max()+1) for grp in groups]
	start = time.time()
	for fr in range(nframes):
		status('electron density',i=fr,looplen=nframes,tag='compute',start=start)
		for g,grp in enumerate(groups):
			for i,z in enumerate(xx[fr][grp]): 
				bincounts_by_group[g][z] += charges[grp[i]]
	xvals = bin_size*np.arange(0,len(bincounts_by_group[0])) - len(bincounts_by_group[0])*bin_size/2.	
	mvecs = np.mean(np.array([vecs[i] for i in range(nframes)]),axis=0)
	scaleconst = np.product(mvecs[:2])*(mvecs[2]/len(xvals))*nframes	
	results,attrs = {},{}
	results['midplane_heights'] = midplane_heights
	#---note that bincoungs_by_group is ordered first by resnamelist then group_regexes
	results['bincounts_by_group'] = np.array(bincounts_by_group)/scaleconst
	for index,group in enumerate(groups): results['groups%d'%index] = np.array(group)
	results['offset'] = np.array(offset)
	attrs['selector'] = selector
	attrs['resnamelist'] = resnamelist
	attrs['bin_size'] = bin_size
	attrs['group_regexes'] = group_regexes
	return results,attrs
