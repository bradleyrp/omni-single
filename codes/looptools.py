#!/usr/bin/env python

"""
Looping tools.
"""

import sys,time
import joblib
from base.tools import framelooper

def basic_compute_loop(compute_function,looper,run_parallel=True,debug=None):
	"""
	Canonical form of the basic compute loop.
	!!! remove this from contacts.py when it works
	"""
	#---send the frame as the debug argument
	if debug!=None and debug!=False:
		fr = debug
		incoming = compute_function(**looper[fr])
		import ipdb;ipdb.set_trace()
		sys.quit()
	start = time.time()
	if run_parallel:
		import joblib
		from joblib import Parallel,delayed
		if joblib.__version__<0.12:
			from joblib.pool import has_shareable_memory
			incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
				delayed(compute_function,has_shareable_memory)(**looper[ll]) 
				for ll in framelooper(len(looper),start=start))
		else:
			incoming = Parallel(n_jobs=8,verbose=10 if debug else 0,require='sharedmem')(
				delayed(compute_function)(**looper[ll]) 
				for ll in framelooper(len(looper),start=start))
		incoming = []
		for ll in framelooper(len(looper)):
			incoming.append(compute_function(**looper[ll]))
	return incoming
