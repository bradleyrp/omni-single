#!/usr/bin/env python

"""
CURVATURE UNDULATION COUPLING DRILLDOWN

round 1: protein_dynamic_single vs protein_dynamic_single_uniform on 100 frames for two extents
round 2: four extents on protein_dynamic_single_uniform with 100 frames
round 3: comparison to all frames for the protein_dynamic_single_uniform
round 4: testing explicit binning
round 5: explicit binning, 100 frames, designed to compare to the legacy code

history:
	started with protein_dynamic_single_catchall and had to manually fix the coordinator
	added protein_dynamic_single_uniform
	added extents (v3,v4) and selected only protein_dynamic_single_uniform to save time
		the protein_dynamic_single runs for thousands of steps
	added v5,v6 to curvature.yaml and generated fields via `make compute` for all frames
	current batch of six tests per simulation is worth making a summary plot
		includes two extents for protein_dynamic_single and four for protein_dynamic_single_uniform
		otherwise lots of other combinations to try
		but currently working on the summary plots
		"""

#---see history notes above. the round name selects different optimmization tests
round_name = ['r1','r2','r3','r4','r2021v01'][-1]

#---all tags in the design loop in curvature.yaml
master_tags = {
	#---several rounds use blurry_explicit binner
	#---first round uses protein_dynamic_single
	'r1':['v1_fixed_extent_1','v2_fixed_extent_2'],
	#---second round uses protein_dynamic_single_uniform
	'r2':['v1_fixed_extent_1','v2_fixed_extent_2','v3_fixed_extent_0.5','v4_fixed_extent_4'],
	'r3':['v5_fixed_extent_all_frames_1','v6_fixed_extent_all_frames_2'],
	'r4':['v1_fixed_extent_1','v2_fixed_extent_2','v3_fixed_extent_0.5','v4_fixed_extent_4'],
	'r2021v01':['v9_fixed_extent_10'],}

from codes.hypothesizer import hypothesizer

def prepare_hypotheses():
	"""
	Generate a list of hypotheses for testing the optimizer.
	We generate the hypotheses in several "rounds" during development of the drilldown.
	"""
	#---start making hypotheses here with the blurry_explicit method
	if round_name in ['r1','r2','r3']:
	    #---default parameters
		opt_params = dict(design={
			'optimize_method':'Nelder-Mead',
			'store_instantaneous_fields':False,
			'store_instantaneous_fields_explicit':False},
			fitting={'initial_kappa':20.0})
		opt_params['design']['binner'] = 'explicit'
		opt_params['design']['weighting_scheme'] = 'blurry_explicit'
		if round_name=='r1':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single']},
				#---! include all elements of curvature_positions or it gets overwritten
				{'route':('design','curvature_positions','nframes'),'values':[100]},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r2':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('design','curvature_positions','nframes'),'values':[100]},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r3':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r4':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		else: raise Exception('fail')
	elif round_name in ['r4']:
		if round_name=='r4':
		    #---default parameters
			opt_params = dict(design={
				'optimize_method':'Nelder-Mead',
				'store_instantaneous_fields':False,
				'store_instantaneous_fields_explicit':False},
				fitting={'initial_kappa':20.0})
			opt_params['design']['binner'] = 'explicit'
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
	elif round_name in ['r2021v01']:
		if round_name=='r2021v01':
		    #---default parameters
			opt_params = dict(design={
				'optimize_method':'Nelder-Mead',
				'store_instantaneous_fields':False,
				'store_instantaneous_fields_explicit':False},
				fitting={'initial_kappa':20.0})
			opt_params['design']['binner'] = 'explicit'
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
	else: raise Exception('fail')
	return hypos
