#!/usr/bin/python

import time
from numpy import *
import itertools
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.timer import checktime
from base.tools import status,framelooper

def ion_binding_norms(**kwargs):

	"""
	Compute normalization constants for weighting the observed lipid-ion-lipid bridges
	"""
	
	dat = kwargs['upstream']['lipid_mesh']
	dat_bridges = kwargs['upstream']['ion_binding_combinator']
	def i2s(*items,**kwargs): return kwargs.pop('delim','.').join([str(i) for i in items])
	nframes = dat['nframes']
	simps = [[dat[i2s(m,fr,'simplices')] for fr in range(nframes)] for m in range(2)]
	#---determine monolayer-specific residue indices
	imono = dat['monolayer_indices']
	nmols = [sum(imono==mn) for mn in range(2)]
	resnames = dat['resnames']
	reslist = list(array(resnames)[sort(unique(resnames,return_index=True)[1])])
	rxm = [[array([where(where(imono==mn)[0]==i)[0][0] 
		for i in where(all((imono==mn,resnames==rn),axis=0))[0]])
		for rn in reslist] for mn in range(2)]
	rxmlook = [zeros(n) for n in nmols]
	for mn in range(2):
		for ri,r in enumerate(rxm[mn]):
			if len(r)>0: rxmlook[mn][r] = ri

	#---count singletons
	#---note that we exclude cholesterol generally so we consult the mesh for singletons
	#---note that we have replaced an ill-advised use of unique (which is unordered)
	if 0: cns = [(v,) for v in unique(resnames)]
	if 0: tally1 = sum(array([[len(i) for ii,i in enumerate(rxm[mn])] for mn in range(2) 
		for fr in range(nframes)]),axis=0)
	cns = list(unique(reslist)[argsort(unique(reslist,return_index=True)[1])])
	tally1 = [0 for c in cns]
	for mn in range(2):
		for ri,rc in enumerate(rxm[mn]):
			if rc != []: tally1[cns.index(reslist[ri])] += float(len(rc)*nframes)
	combonames0 = array([cns]).T.astype(str)

	#---count possible pairs
	cns = dat_bridges[i2s(1,'combonames')]
	combolookup = sum([array(cns)==r for r in reslist],axis=2).T
	tally2 = zeros(len(combolookup))
	for mn in range(2):
		for fr in range(nframes):
			status('[COMPUTE] bridging norm nn=2',i=fr,looplen=nframes)
			
			#---removed neighborlist from mesh objects so we rely on simplices to generate the neighborlist
			if 0: nl = neighborunpack(dat[i2s(fr,mn,'nl')])
			simplices = dat[i2s(mn,fr,'simplices')]
			nl = [unique(simplices[where(any(simplices==v,axis=1))[0]]) 
				for v in range(len(dat[i2s(mn,fr,'points')]))]
			
			gi = dat[i2s(mn,fr,'ghost_ids')].astype(int)
			#---pair types is list of potential pairs from each vertex
			pair_types = array([transpose((tile(rxmlook[mn][ii],len(i)),rxmlook[mn][gi[i]])) 
				for ii,i in enumerate(nl[:nmols[mn]])])
			counts = concatenate([array([sum(t==r,axis=1) 
				for r in range(len(reslist))]).T for t in pair_types])
			ans = [where(all(c==combolookup,axis=1))[0][0] for c in counts]
			tally2[unique(ans)] += array([sum(ans==i)/2. for i in unique(ans)])
	combonames1 = array(cns)

	#---count the possible triple bridges
	cns = dat_bridges[i2s(2,'combonames')]
	combolookup = sum([array(cns)==r for r in reslist],axis=2).T
	tally3 = zeros(len(combolookup))
	for mn in range(2):
		for fr in range(nframes):
			status('[COMPUTE] bridging norm nn=3',i=fr,looplen=nframes)
			gi = dat[i2s(mn,fr,'ghost_ids')].astype(int)
			s = simps[mn][fr]
			counts = sum([rxmlook[mn][gi[s]]==r for r in range(len(reslist))],axis=2).T
			ans = [where(all(c==combolookup,axis=1))[0][0] for c in counts]
			tally3[unique(ans)] += array([sum(ans==i) for i in unique(ans)])
	combonames2 = array(cns)
		
	attrs,results = {},{'nn0':tally1,'nn1':tally2,'nn2':tally3,
		'cns0':combonames0,'cns1':combonames1,'cns2':combonames2}
	return results,attrs

