#!/usr/bin/env python

import time,itertools
import numpy as np
import scipy
import scipy.stats
import numpy as np

from base.tools import status
from codes.looptools import basic_compute_loop

"""
Analyze the triples and pairs on the lipid mesh.
"""

i2s = lambda mn,fr,key : '%d.%d.%s'%(mn,fr,key)

def counter(mn,fr,nn,simplices=None,random=False):
	"""..."""
	global dat,results,rxmlook
	reslist = results['reslist']
	if type(simplices)==type(None): 
		simplices = dat[i2s(mn,fr,'simplices')]
		gids = dat[i2s(mn,fr,'ghost_ids')]
		simplices_gids = gids[simplices]
	# allow the calling function to supply the ghost-indexed simplices directly
	else: simplices_gids = simplices
	if nn==2:
		links = np.array(list(set([tuple(np.sort(j)) for j in 
			np.concatenate([np.transpose((simplices_gids[:,i],simplices_gids[:,(i+1)%3])) 
			for i in range(3)])])))
		#---uniquely identify each link type
		if not random:
			links_type = np.transpose([np.sum(rxmlook[mn][links]==i,axis=1) 
				for i in range(len(reslist))])
		else:
			#---sidestepping the probability questions and randomizing the identifier
			#---...which will effectively randomize the identities of the vertices
			#---...and then later we will take the average of these and then use it to see if observed
			#---...links are more or less likely to appear in our graph than in a random one
			#---modified for one-randomization-per-trial instead of per frame
			#rxmlook_rand = [np.random.permutation(r) for r in rxmlook]
			rxmlook_rand = results['rxmlook_rand']
			links_type = np.transpose([np.sum(rxmlook_rand[mn][links]==i,axis=1) 
				for i in range(len(reslist))])
		links_type_str = [''.join(['%s'%s for s in i]) for i in links_type]
		counts = dict(scipy.stats.itemfreq(links_type_str))
	elif nn==3:
		#---identify each triple type
		triples = simplices_gids
		triple_type = np.transpose([np.sum(rxmlook[mn][triples]==i,axis=1) for i in range(len(reslist))])
		if not random:
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type]
		else:
			#---randomize the triples
			#rxmlook_rand = [np.random.permutation(r) for r in rxmlook]
			rxmlook_rand = results['rxmlook_rand']
			triple_type_wrong = np.transpose([np.sum(rxmlook_rand[mn][triples]==i,axis=1) 
				for i in range(len(reslist))])
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type_wrong]
		counts = dict(scipy.stats.itemfreq(triple_type_str))
	else: raise Exception('counter only designed for nn in [2,3]')
	return np.array([int(counts[i]) if i in counts else 0 for i in results['combo_lookup_str_%d'%nn]])

def lipid_mesh_partners(**kwargs):
	"""
	Compute bilayer midplane structures for studying undulations.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---! deprecated random trials: n_trials = kwargs['calc']['specs']['n_trials']
	do_randomize = False
	#---globals for parallel
	global dat,results,rxmlook
	dat = kwargs['upstream']['lipid_mesh']
	nmols = [int(dat[i2s(mn,0,'nmol')]) for mn in range(2)]
	nframes = int(dat['nframes'])
	resnames = dat['resnames']
	attrs,results = {},{}

	#---code adapted from previous version at simuluxe and binding_combinator
	resnames = np.array(dat['resnames'])
	nframes = dat['nframes']
	lipids = np.array(list(resnames[np.sort(np.unique(resnames,return_index=True)[1])]))
	reslist = list(np.array(resnames)[np.sort(np.unique(resnames,return_index=True)[1])])
	results['reslist'] = reslist
	#---collect statistics for pairs and triples
	for nn in [2,3]:
		combos = np.array([''.join(j) for j in 
			itertools.product(''.join([str(i) for i in range(nn+1)]),repeat=len(lipids)) 
			if sum([int(k) for k in j])==nn])
		combonames = [tuple(v) for v in [np.concatenate([[lipids[ww]]*int(w) 
			for ww,w in enumerate(l)]) for l in combos]]
		results['combos_%d'%nn] = combos
		results['combonames_%d'%nn] = combonames
		combolookup = np.sum([np.array(combonames)==r for r in reslist],axis=2).T
		combolookup_str = [''.join(['%s'%s for s in i]) for i in combolookup]
		results['combo_lookup_%d'%nn] = combolookup
		results['combo_lookup_str_%d'%nn] = combolookup_str

	#---determine monolayer-specific residue indices
	imono = dat['monolayer_indices']
	nmols = [np.sum(dat['monolayer_indices']==i) for i in range(2)]
	resnames = np.array(dat['resnames'])
	rxm = [[
		np.array([np.where(np.where(imono==mn)[0]==i)[0][0] 
		for i in np.where(np.all((imono==mn,resnames==rn),axis=0))[0]])
		for rn in reslist] for mn in range(2)]
	rxmlook = [np.zeros(n) for n in nmols]
	for mn in range(2):
		for ri,r in enumerate(rxm[mn]):
			if r != []: rxmlook[mn][r] = ri

	#---count in parallel
	counts_trials = dict([(nn,[]) for nn in [2,3]])
	counts_observed = dict([(nn,None) for nn in [2,3]])
	for nn in [2,3]:
		status('observations for nn=%d'%(nn),tag='compute')
		looper = [dict(fr=fr,mn=mn,nn=nn) for mn in range(2) for fr in range(nframes)]
		incoming = basic_compute_loop(counter,looper)
		#---reindex data mn,fr,combo
		counts_observed[nn] = np.concatenate(incoming).reshape((2,nframes,len(results['combonames_%d'%nn])))
		if do_randomize:
			for trial in range(n_trials):
				results['rxmlook_rand'] = [np.random.permutation(r) for r in rxmlook]
				status('randomize trial for nn=%d trial=%d/%d'%(nn,trial+1,n_trials),tag='compute')
				looper = [dict(fr=fr,mn=mn,nn=nn,random=True) for mn in range(2) for fr in range(nframes)]
				incoming = basic_compute_loop(counter,looper)
				counts_trials[nn].append(
					np.concatenate(incoming).reshape((2,nframes,len(results['combonames_%d'%nn]))))
	if do_randomize: counts_random = dict([(nn,np.concatenate([counts_trials[nn]])) for nn in [2,3]])
	#---pack
	for nn in [2,3]:
		if do_randomize: results['counts_random_%d'%nn] = np.array(counts_random[nn])
		results['counts_observed_%d'%nn] = np.array(counts_observed[nn])
	results.pop('rxmlook_rand',None)
	#---save rxmlook for counting lipids
	for mn in range(2): results['monolayer_residues_%d'%mn]	= rxmlook[mn]
	return results,attrs	
