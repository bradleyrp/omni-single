#!/usr/bin/env python

"""
CURVATURE COUPLING LOADER
Import GROMACS membrane data for the curvature coupling calculation.
This function is alone so it can be swapped out with other importers.
"""

import numpy as np
from base.tools import status
from curvature_coupling.tools import fft_field
from base.store import plotload

def curvature_coupling_loader_membrane(data): 
	"""
	Receive the undulation data and prepare the meshes for the curvature coupling calculation.
	"""
	#---point heights into "memory"
	status('populating memory',tag='load')
	memory = {}
	for sn in data['undulations'].keys():
		if (sn,'hqs') not in memory:
			dat = data['undulations'][sn]['data']
			vecs = dat['vecs']
			mesh = dat['mesh']
			midplane = mesh.mean(axis=0)
			zmeans = midplane.reshape((midplane.shape[0],-1)).mean(axis=1)
			midplane = np.array([i-zmeans[ii] for ii,i in enumerate(midplane)])
			hqs = fft_field(midplane)
			memory[(sn,'hqs')] = hqs
			memory[(sn,'vecs')] = vecs
	return memory

def curvature_coupling_loader_protein(data): 
	"""
	Receive the undulation data and prepare the meshes for the curvature coupling calculation.
	"""
	#---note that you can load other upstream data with the following command:
	#---... data_other,calc_other = work.plotload('undulations',status_override=True,
	#---... 	sns=data['protein_abstractor'].keys())
	#---note also that there are two ways to manipulate incoming data for InvestigateCurvature: (1) alter
	#---...the data here or (2) duplicate the curvature_undulation_coupling.py calculation with custom 
	#---...upstream data types and custom imports
	return data['protein_abstractor']
