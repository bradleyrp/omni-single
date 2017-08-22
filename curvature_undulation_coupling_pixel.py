#!/usr/bin/env python

"""
Current optimized version of the curvature-undulation coupling method.
"""

import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature

def curvature_undulation_coupling_pixel(**kwargs):
	"""
	Supervise the curvature-undulation coupling calculation.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---instantiate the calculation	
	ic = InvestigateCurvature(sn=sn,work=kwargs['workspace'],
		design=kwargs['calc']['specs'].get('design',{}),
		protein_abstractor='pixel_method_only',
		undulations=kwargs['upstream']['undulations'])
	#---repackage the data
	attrs,result = ic.finding['attrs'],ic.finding['result']
	return result,attrs
