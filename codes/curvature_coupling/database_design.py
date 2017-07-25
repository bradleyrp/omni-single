#!/usr/bin/env python

"""
"""

database_def = {}
database_def['hypothesis'] = dict(
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
		('fallback',str),),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein','multipoint'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0 or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},],
	anticipate = ['kappa','gamma','vibe','error'],)

database_def['field'] = dict(
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
		('theta',float),),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0  or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},
		lambda x:{'curvature':1.0} if x['curvature']!=0.0 else {},
		lambda x:{'sn':x['fallback'],'fallback':'none'} if work.meta[x['sn']].get('nprots',1)==0 
			and x['curvature'] != 0.0 and x['mapping']=='protein' else {},],
		#---test specs throw errors when the fail but change specs simply "reduce" like hypotheses
		#---incorrect: lambda x:{'mapping':'single'} if work.meta[x['sn']]['nprots']==0 else {},
	anticipate = ['kappa','gamma','vibe','error'],)
