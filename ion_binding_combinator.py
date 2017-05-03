#!/usr/bin/python

import time
from numpy import *
import itertools
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.timer import checktime
from base.tools import status,framelooper

def ion_binding_combinator(**kwargs):

	"""
	Compute bridges.
	"""

	sn = kwargs['sn']
	dat = kwargs['upstream']['ion_binding']
	resnames = dat['resnames']
	pas = dat['partners_atoms']
	lipid_distances = dat['lipid_distances']
	nframes = dat['nframes']
	zonecut = kwargs['calc']['specs']['zonecut']
	results,attrs = {},{}
	attrs['zonecut'] = zonecut
	#---zonecut is angstroms while lipid_distances is nm
	zonecut = zonecut/10.
	i2s2 = lambda *items: '.'.join([str(i) for i in items])
	#---previous method may have created disorder downstream
	if 0: lipids = unique(resnames[unique([tuple(i) for fr in range(nframes) for i in pas[fr]])])
	lipids = array(list(resnames[sort(unique(resnames,return_index=True)[1])]))
	for nn in range(3):
		combos = array([''.join(j) for j in itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) if sum([int(k) for k in j])==nn+1])
		combonames = [tuple(v) for v in [concatenate([[lipids[ww]]*int(w) for ww,w in enumerate(l)]) for l in  combos]]
		#---! problematic method excised below
		#---! cind = lambda a : where(combos==''.join([str(sum(array(a)==i)) for i in lipids]))[0][0]

		wcs = zeros((nframes,len(combos)))
		st = time.time()
		status('[COMPUTE] combinator '+sn)
		#import pdb;pdb.set_trace()
		for fr in range(nframes):
			status('[COMPUTE] combinator nn='+str(nn+1),i=fr,looplen=nframes,start=st)
			parts = resnames[pas[fr,where(sum(lipid_distances[fr]<zonecut,axis=1)==nn+1)[0]]][:,:nn+1]
			#---! wcs[fr] = array([sum(array([cind(j) for j in parts])==i) for i in range(len(combos))])
			wcs[fr] = array([sum(array([where(combos==''.join([str(sum(array(j)==k)) for k in lipids]))[0][0] for j in parts])==i) for i in range(len(combos))])
		results[i2s2(nn,'wcs')] = wcs
		results[i2s2(nn,'combos')] = combos
		results[i2s2(nn,'combonames')] = array(combonames)
	return results,attrs

