#!/usr/bin/env python

"""
Current optimized version of the curvature-undulation coupling method.
"""

import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature

def curvature_undulation_coupling_controls(**kwargs):
	"""
	Supervise the curvature-undulation coupling calculation.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---special processing for matching the protein positions to the controls
	ortho_specs = work.calcs['curvature_undulation_coupling_controls']['ortho_specs']
	#---imaginary_fields should be a string or list of collections
	imaginaries = ortho_specs['imaginary_fields']
	if type(imaginaries)!=list: imaginaries = [imaginaries]
	protein_sources = work.get_simulations_in_collection(*imaginaries)
	attrs,result = {},{}
	#---save the comparison simulations in the order of the following loop
	attrs['compare_sns'] = protein_sources
	#---loop over all imaginary protein fields for the control
	for snum,protein_source in enumerate(protein_sources):
		status('running control %s with imaginary protein fields from %s'%(
			sn,protein_source),i=snum,looplen=len(protein_sources))
		data_prot,calc_prot = work.plotload('protein_abstractor',sns=[protein_source])
		#---repackage to simulate the case where the protein data came from this simulation
		protein_abstractor = {sn:dict(data=data_prot[protein_source]['data'])}
		#---reindex the protein points so there are enough of them
		keys_to_reindex = ['points_all','vecs']
		nframes_prot = len(protein_abstractor[sn]['data']['vecs'])
		nframes = kwargs['upstream']['undulations']['nframes']
		reindex = (np.arange(nframes)%nframes_prot)[:nframes]
		for key in keys_to_reindex:
			protein_abstractor[sn]['data'][key] = protein_abstractor[sn]['data'][key][reindex]
		#---instantiate the calculation	
		ic = InvestigateCurvature(sn=sn,work=kwargs['workspace'],
			design=kwargs['calc']['specs'].get('design',{}),
			protein_abstractor=data_prot[protein_source]['data'],
			undulations=kwargs['upstream']['undulations'])
		#---note the following incoming keys for attrs: spec, bundle
		#---note the following incoming keys for result: qs cf drop_gaussians_points cf_first ratios x jac
		#---repackage the data with indices
		for key,val in ic.finding['attrs'].items(): attrs['%d,%s'%(snum,key)] = val
		for key,val in ic.finding['result'].items(): result['%d,%s'%(snum,key)] = val
		del protein_abstractor,ic,data_prot,calc_prot
	return result,attrs
