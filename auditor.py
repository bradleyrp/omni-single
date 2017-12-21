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
		self.root = 'calcs'
		self.ledger_fn = kwargs.pop('ledger','audit.yaml')
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		status('welcome to the auditor')
		self.ledger = os.path.join(self.root,self.ledger_fn)
		if not os.path.isfile(self.ledger): raise Exception('cannot find %s'%self.ledger)
		else: 
			with open(self.ledger) as fp: self.raw = yaml.load(fp.read())
		# print everything
		asciitree(self.raw)
		self.interpret()

	def interpret(self):
		"""
		Collect all files on disk and classify them.
		"""
		self.files = []
		for (path,dns,fns) in os.walk(self.root):
			self.files.extend([os.path.relpath(os.path.join(path,f),self.root) for f in fns])
		for kind,details in self.raw.get('meta',{}).get('kinds',{}).items():
			if set(details.keys())=={'pattern','do'}:
				if details['do']=='ignore':
					self.files = [i for i in self.files if not re.match(details['pattern'],i)]
				elif details['do']=='account':
					accounted = [i for i in self.files if re.match(details['pattern'],i)]
					unaccounted = [i for i in self.files 
						if i in accounted and i not in self.raw.get(kind,[])]
					if unaccounted: raise Exception('class %s is missing: %s'%(kind,unaccounted))
					self.files = [i for i in self.files if i not in accounted]
				else: raise Exception('do %s?'%details['do'])
			elif set(details.keys())=={'from','do'}:
				if details['do']=='ignore':
					self.files = [i for i in self.files if i not in self.raw.get(kind,[])]
				else: raise Exception('do %s?'%details['do'])
			else: raise Exception('cannot interpret %s,%s'%(kind,details))
		if self.files:
			for f in sorted(self.files): print(f)
			raise Exception('unaccounted files are listed above')

