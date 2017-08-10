#!/usr/bin/env python

"""
Current optimized version of the curvature-undulation coupling method.
"""

import numpy as np
import codes.curvature_coupling
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from codes.readymade_meso_v1 import import_curvature_inducer_points,import_membrane_mesh

def curvature_undulation_coupling_dex(**kwargs):
	"""
	Supervise the curvature-undulation coupling calculation.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	#---retrieve the membrane and curvature-inducer locations
	membrane_mesh = import_membrane_mesh(calc=calc,work=work,sn=sn)
	curvature_inducer_points = import_curvature_inducer_points(calc=calc,work=work,sn=sn)
	#---instantiate the calculation	
	ic = InvestigateCurvature(sn=sn,work=kwargs['workspace'],
		design=kwargs['calc']['specs'].get('design',{}),
		protein_abstractor=curvature_inducer_points,
		undulations=membrane_mesh)
	#---repackage the data
	attrs,result = ic.finding['attrs'],ic.finding['result']
	return result,attrs
