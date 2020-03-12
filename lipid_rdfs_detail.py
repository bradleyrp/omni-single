#!/usr/bin/env python

"""
Quick port from lipid_rdfs to get convergence.
"""

import numpy as np
import scipy
import scipy.spatial
from base.compute_loop import basic_compute_loop
import MDAnalysis
import itertools

def compute_rdf_imaged(coords,vec,groups,bins,dims=2):
	"""Use the "imaged" method to compute the RDF."""
	pd = np.array([scipy.spatial.distance.cdist(np.array([coords[:,d]]).T,
		np.array([coords[:,d]]).T) for d in range(dims)])
	if dims!=2: raise Exception('dev')
	# no conditionals here. instead we replicate in each direction. previously we 
	# ... corrected distances above a half-box but this makes tabulation difficult later
	shifts = np.array([(i*vec[0],j*vec[1]) for i in [-1,0,1] for j in [-1,0,1]]).T
	nmols = coords.shape[0]
	# tile the dimension-wise distances
	distances_breakout = np.tile(pd,((shifts.shape[1],1,1,1))).transpose((1,0,2,3))
	shifter = np.tile(shifts,(nmols,nmols,1,1)).transpose((2,3,0,1))
	distances = distances_breakout + shifter
	normed = np.linalg.norm(distances,axis=0)
	# for memory reasons we histogram each time the distances are calculated, hence we subsample by group
	reformed = dict([(key,normed[:,g1][...,g2].reshape(-1)) for key,(g1,g2) in groups.items()])
	# histogram each group
	grammed = dict([(key,np.histogram(val[val>0],bins=bins)[0]) for key,val in reformed.items()])
	# return a dictionary of histogrammed distances
	return grammed

def lipid_rdfs_detail(**kwargs):
	"""
	Compute 2D lipid radial distribution functions.
	"""
	dat = kwargs['upstream']['lipid_abstractor']
	imono = dat['monolayer_indices']
	points = dat['points']
	resnames = dat['resnames']
	vecs = dat['vecs']
	attrs,result = {},{}
	# global scanrange from specs to reduce the distances data
	cutoff = kwargs['calc']['specs']['cutoff']
	binsize = kwargs['calc']['specs']['binsize']
	scanrange = np.arange(0,cutoff,binsize)
	# loop over monolayers
	for mn in range(2):
		# prepare pairs
		resnames_u = np.unique(resnames[np.where(imono==mn)])
		pairs = ([([r],[r]) for r in resnames_u]+
			[([i],[j]) for i,j in itertools.combinations(resnames_u,2)]+
			[(resnames_u,resnames_u)])
		pairnames = [(i[0],j[0]) for i,j in pairs[:-1]]+[('all lipids','all lipids')]
		pairspec = dict(pairs=pairs,pairnames=pairnames,groups={})
		# get coordinates for this leaflet
		coords = points[:,np.where(imono==mn)[0],:2]
		# tabulate rows for each group
		for pair,pairname in zip(pairspec['pairs'],pairspec['pairnames']):
			group_1 = np.where(np.in1d(resnames[np.where(imono==mn)],pair[0]))[0]
			group_2 = np.where(np.in1d(resnames[np.where(imono==mn)],pair[1]))[0]
			pairspec['groups'][pairname] = (group_1,group_2)
		nframes = len(coords)
		looper = [dict(vec=vecs[fr],coords=coords[fr],bins=scanrange,groups=pairspec['groups']) 
			for fr in range(nframes)]
		incoming = np.array(basic_compute_loop(compute_rdf_imaged,looper))
		#! note that we no longer take the mean over axis=0 here, compared to lipid_rdfs
		#!   which allows us to later check the convergence of the RDFs
		obs = dict([(pair,np.array([i[pair] for i in incoming])) 
			for pair in pairspec['pairnames']])
		# package
		tag = '_mn%s'%mn
		attrs['cutoff'] = cutoff
		attrs['binsize'] = binsize
		attrs['pairs'+tag] = [[tuple(j) for j in i] for i in pairspec['pairs']]
		attrs['pairnames'+tag] = [tuple(i) for i in np.array(pairnames)]
		# save in sequence named by pairnames
		for pairname in attrs['pairnames'+tag]:
			result['counts'+tag+'_%s:%s'%tuple(pairname)] = obs[pairname]
	result['resnames'] = dat['resnames']
	result['monolayer_indices'] = dat['monolayer_indices']
	result['total_area'] = np.product(vecs.mean(axis=0)[:2])
	return result,attrs	
