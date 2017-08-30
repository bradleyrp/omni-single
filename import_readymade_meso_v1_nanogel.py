#!/usr/bin/env python

import time
from numpy import *
import MDAnalysis
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
from base.timer import checktime
from codes.mesh import *
from codes.readymade_meso_v1 import import_nanogel_positions

def import_readymade_meso_v1_nanogel(**kwargs):
	"""
	PROTEIN ABSTRACTOR for readymade_meso_v1 data.
	Save the positions of the nanogel.
	Note that most of the work is done by an external function, which will be useful in the event that the 
	upstream data format changes slightly, in which case we can write a new function.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	mapping_style = kwargs['calc'].get('specs',{}).get('mapping_style',None)
	if not mapping_style: raise Exception('set the mapping style in the calculation specs')
	#---import and interpret the mesh points
	position_data = import_nanogel_positions(sn=sn,calc=calc,work=work)
	framenos,points = [position_data[j] for j in ['framenos','points']]
	nframes = len(framenos)
	#---pack
	attrs,result = {},{}
	#---we use a placeholder for this data type so downstream calculations are not confused
	result['resnames'] = 'readymade_meso_v1'
	result['names'] = 'readymade_meso_v1'
	result['vecs'] = 'readymade_meso_v1'
	result['nframes'] = array(nframes)
	result['frame_numbers'] = array(framenos)
	result['points'] = 'readymade_meso_v1'
	#---treat the nanogel as a single "protein" or curvature inducer. note that this method is compatible 
	#---...with two different downstream curvature coupling methods: the neighborhood method can apply to 
	#---...the entire nanogel or each bead can be labelled as a hotspot (or a separate step can identify the #---...most-contacting hotspots). any further specificity in the dropped gaussians would require 
	#---...additional explication
	if mapping_style=='single_inducer':
		#---reformulate the points so the dimensions are frames by 1 by beads by 3 (xyz)
		result['points_all'] = np.transpose([points],(1,0,2,3))
	else: raise Exception('invalid mapping_style %s'%mapping_style)
	return result,attrs	
