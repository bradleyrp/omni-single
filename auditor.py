#!/usr/bin/env python

"""
AUDITOR
Catalogs and provides the status for your calculations.
"""

import os,collections,glob,re
from datapack import asciitree
from base.tools import status,str_types
import yaml

# trick for OrderedDict in  YAML
# via: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
#! solution above may not be durable in Python 3
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_representer(dumper, data):
	return dumper.represent_dict(data.iteritems())
def dict_constructor(loader, node):
	return collections.OrderedDict(loader.construct_pairs(node))
yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

class CalcsAuditor:
	"""
	Read an audit file and print a report.
	"""
	_ignore_files = ['.gitignore']
	def __init__(self,**kwargs):
		self.calcs_dn = kwargs.pop('calcs_dn','calcs')
		self.ledger_fn = kwargs.pop('ledger','audit.yaml')
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		status('welcome to the auditor')
		self.ledger = os.path.join(self.calcs_dn,self.ledger_fn)
		if not os.path.isfile(self.ledger): raise Exception('cannot find %s'%self.ledger)
		else: 
			with open(self.ledger) as fp: self.raw = yaml.load(fp.read())
		# print everything
		asciitree(self.raw)
		# only two default classes
		self.classes = dict(extrana=[],calculations=[])
		self.interpret()

	def interpret(self):
		# prepare glob patterns for organizing the files
		self.patterns = self.raw.get('meta',{}).get('patterns',{})
		# collect all files
		self.files = []
		# sources can be overridden in meta but they must be relative to self.calcs_dn
		for source in self.raw.get('meta',{}).get('sources',['.','specs']):
			for (path,dns,fns) in os.walk(os.path.join(self.calcs_dn,source)):
				self.files.extend([os.path.relpath(os.path.join(path,f),self.calcs_dn) for f in fns])
				break
		# we use regex instead of globs because the sources list controls input directories
		#! we could still use a glob to filter the sources? seems like too many options
		#self.patterns_matched = dict([(k,glob.glob(os.path.join(self.calcs_dn,v))) 
		#	for k,v in self.patterns.items()])
		for key in self.patterns: 
			if key not in self.classes: self.classes[key] = []
		self.strays = []
		for fn in self.files:
			matches = [key for key,val in self.patterns.items() if re.match(val,fn)]
			if len(matches)==0: self.strays.append(fn)
			elif len(matches)>1: raise Exception('file %s matches multiple globs: %s'%matches)
			else: self.classes[matches[0]].append(fn)
		if self.classes['extrana']: raise Exception('cannot populate extrana: %s'%self.classes['extrana'])
		self.classes['extrana'] = self.raw.get('extrana',[])
		self.strays = [fn for fn in self.strays if fn not in self.classes['extrana'] and fn not in self._ignore_files]
		if self.strays:
			for i in sorted(self.strays): print(i)
			raise Exception('stray files listed above')
		for kind,fns in self.classes.items():
			matched_not_named = [f for f in fns if f not in self.raw.get(kind,[])]
			if any(matched_not_named):
				for k in matched_not_named: print(k)
				raise Exception('the above files were not matched to %s'%kind)
		import ipdb;ipdb.set_trace()
