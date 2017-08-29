#!/usr/bin/env python

"""
Calculation which interfaces to InvestigateCurvature which handles the curvature undulation coupling 
calculation. 
This version is modified for incoming mesoscale data.
"""

import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature

def curvature_undulation_coupling_dex(**kwargs):
	"""
	Supervise the curvature-undulation coupling calculation.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---retrieve the membrane and curvature-inducer locations
	#---note that in contrast to the standard curvature_undulation_coupling.py we retrieve the 
	#---...upstream data in a separate function via readymade_meso_v1. in the standard method, both the
	#---...protein_abstractor and undulations data come into this function via e.g.
	#---...kwargs['upstream']['protein_abstractor'] which gets sent to InvestigateCurvature
	#---note that users who wish to manipulate the incoming data can do so in a custom copy of 
	#---...curvature_undulation_coupling.py (i.e. this script) or in the loader functions
	#---note also that users could point to upstream data as we have done here (data which are different than
	#---...in the standard method this was derived from) or import alternate data directly here. we have 
	#---...chosen to implement an upstream calculation to package the data so it enters the pipeline as early
	#---...as possible however all of the methods described here are viable
	#---instantiate the calculation	
	ic = InvestigateCurvature(sn=sn,work=kwargs['workspace'],
		design=kwargs['calc']['specs'].get('design',{}),
		fitting=kwargs['calc']['specs'].get('fitting',{}),
		protein_abstractor=kwargs['upstream']['import_readymade_meso_v1_nanogel'],
		undulations=kwargs['upstream']['import_readymade_meso_v1_membrane'])
	#---repackage the data
	attrs,result = ic.finding['attrs'],ic.finding['result']
	return result,attrs
