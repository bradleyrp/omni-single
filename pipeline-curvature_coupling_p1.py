#!/usr/bin/env python

if False:
	execfile('./omni/base/header.py')
	import sys;sys.path.insert(0,'calcs')
	from base.store import plotload
	import codes.undulate as undulate
	execfile('./calcs/specs/figures.py')
	execfile('./calcs/specs/colors.py')
	import numpy as np
	import multiprocessing as mp
	import time,re
	import scipy
	import json
	from base.store import store
	from base.store import load
	import os
	if 'readline' not in globals():
		execfile(os.environ['PYTHONSTARTUP'])
		import os

#---SWEEPS
#-------------------------------------------------------------------------------------------------------------

#---track various versions of the calculation
tweak_menu = []

#---wish fulfillment? (note that 40 is the original 29)
#---turned out to be correct (thank the optimizer)
tweak_menu += [('v41',{
	'residual_form':'log',
	'sweep_code':'sweep-v1',
	'energy_form':energy_form,
	'curvature_field_spec':'cfv5',
	'inner_sign':-1.0,'fft_flag':'complex',
	'init':[20,0.0,2.0],'hicut':2.0,'opt':optimize_hypothesis,
	'cut_curve':(0.0,0.032),'error_ceiling':0.0012,
	'kappa_gamma_prefactor':{'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}})]

tweak_menu += [('v42',{
	'residual_form':'log',
	'sweep_code':'sweep-v2',
	'energy_form':energy_form,
	'curvature_field_spec':'cfv5',
	'inner_sign':-1.0,'fft_flag':'complex',
	'init':[20,0.0,2.0],'hicut':1.0,'opt':optimize_hypothesis,
	'cut_curve':(0.0,0.032),'error_ceiling':0.0014,
	'kappa_gamma_prefactor':{'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}})]

tweak_menu += [('v40',{
	'residual_form':'log',
	'sweep_code':'sweep-v1',
	'energy_form':energy_form,
	'curvature_field_spec':'cfv5',
	'inner_sign':-1.0,'fft_flag':'complex',
	'init':[20,0.0,2.0],'hicut':1.0,'opt':optimize_hypothesis,
	'cut_curve':(0.0,0.05),'error_ceiling':0.0014,
	'kappa_gamma_prefactor':{'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}})]

#---DATABASE
#-------------------------------------------------------------------------------------------------------------

database = {}
database['hypothesis'] = dict(
	table_name = 'hypothesis',
	data_map = (
		('sn',str),
		('timecode',str),		
		('curvature',float),
		('sigma_a',float),
		('sigma_b',float),
		('isotropy',float),
		('mapping',str),
		('motion',str),
		('theta',float),
		('vibe',float),
		('gamma',float),
		('kappa',float),
		('error',float),
		('fallback',str),
		),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein','multipoint'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],
		],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0 or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},
		],
	anticipate = ['kappa','gamma','vibe','error'],
	)

database['field'] = dict(
	table_name = 'field',
	data_map = (
		('sn',str),
		('timecode',str),		
		('curvature',float),
		('sigma_a',float),
		('sigma_b',float),
		('isotropy',float),
		('mapping',str),
		('motion',str),
		('theta',float),
		),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],
		],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0  or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},
		lambda x:{'curvature':1.0} if x['curvature']!=0.0 else {},
		lambda x:{'sn':x['fallback'],'fallback':'none'} if work.meta[x['sn']]['nprots']==0 
			and x['curvature'] != 0.0 and x['mapping']=='protein' else {},
		#---test specs throw errors when the fail but change specs simply "reduce" like hypotheses
		#---incorrect: lambda x:{'mapping':'single'} if work.meta[x['sn']]['nprots']==0 else {},
		],
	anticipate = ['kappa','gamma','vibe','error'],
	)

#---CLASS
#-------------------------------------------------------------------------------------------------------------

class CurvatureCoupler:

	#---note which "round" of the calculation we are on
	roundsig = 'hcv3'

	#---run "ln -s curvature_coupling-hcv3-cfv4 curvature_coupling-hcv3-cfv9" to use the current cfv4
	curvature_panels = {
		#---the cfv4 was first generated with hcv2 and linked in to the following codes because it is 32GB
		'cfv4':[{'motion':j,'mapping':i} for i in ['single','protein'] for j in ['static','dynamic']],
		#---ignore static for now because it is largely redundant
		'cfv5':[{'motion':j,'mapping':i} for i in ['single','protein'] for j in ['static','dynamic'][1:]],
		}

	def __init__(self,**kwargs):

		"""
		Perform a single curvature coupling calculation along with associated plots.
		"""

		self.hcv = kwargs['hcv']

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---parse arguments
re_hcv = '^hcv=(.+)$'
try: hcv, = [re.findall(re_hcv,i)[0] for i in sys.argv if re.match(re_hcv,i)]
except: raise Exception('you must supply a hard-coded version e.g. hcv="v42"')

#---instantiate the calculation
cc = CurvatureCoupler(hcv=hcv)
