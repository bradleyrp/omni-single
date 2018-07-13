#!/usr/bin/env python

import sys,time
from numpy import *
from joblib import Parallel,delayed
# from joblib.pool import has_shareable_memory

from codes.mesh import *
from base.timer import checktime
from base.tools import status,framelooper

def lipid_mesh(**kwargs):

	"""
	Compute monolayer mesh objects.
	"""

	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	dat = kwargs['upstream']['lipid_abstractor']
	resnames = dat['resnames']
	monolayer_indices = dat['monolayer_indices']
	nframes = dat['nframes']
	debug = kwargs.pop('debug',False)
	kwargs_out = dict(curvilinear=calc.get('specs',{}).get('curvilinear',False))

	#---parallel
	mesh = [[],[]]
	if debug: 
		mn,fr = 0,10
		makemesh(dat['points'][fr][where(monolayer_indices==mn)],dat['vecs'][fr],
			debug=True,**kwargs_out)
		sys.exit(1)
	for mn in range(2):
		start = time.time()
		mesh[mn] = Parallel(n_jobs=work.nprocs,verbose=0)(
			delayed(makemesh)(
				dat['points'][fr][where(monolayer_indices==mn)],dat['vecs'][fr],**kwargs_out)
			for fr in framelooper(nframes,start=start,text='monolayer %d, frame'%mn))
	checktime()

	#---pack
	attrs,result = {},{}
	result['nframes'] = array(nframes)
	result['vecs'] = dat['vecs']
	result['resnames'] = resnames
	result['monolayer_indices'] = monolayer_indices
		
	#---pack mesh objects
	#---keys include: vertnorms simplices nmol facenorms gauss points vec ghost_ids mean principals areas
	keylist = mesh[0][0].keys()
	for key in keylist:
		for mn in range(2):
			for fr in range(nframes): 
				result['%d.%d.%s'%(mn,fr,key)] = mesh[mn][fr][key]		
				
	return result,attrs	
