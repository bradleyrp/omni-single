#!/usr/bin/env python

"""
AUDITOR
Catalogs and provides the status for your calculations.
"""

import os,collections,glob,re
from datapack import asciitree
from base.tools import status,str_types
import yaml

### YAML extensions
#! recommend putting these in a central place
# trick for no duplicate keys via: https://gist.github.com/pypt/94d747fe5180851196eb
# note that this must precede the ordered dictionary trick
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def no_dups(loader,node,deep=False):
	mapping = {}
	for key_node, value_node in node.value:
		key = loader.construct_object(key_node, deep=deep)
		value = loader.construct_object(value_node, deep=deep)
		if key in mapping:
			raise yaml.constructor.ConstructorError(
				'while mapping',node.start_mark,'found duplicate key %s'%key,key_node.start_mark)
		mapping[key] = value
	return loader.construct_mapping(node, deep)
yaml.add_constructor(_mapping_tag,no_dups)
# trick for OrderedDict in  YAML
# via: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
#! solution above may not be durable in Python 3
def dict_representer(dumper, data):
	return dumper.represent_dict(data.iteritems())
def dict_constructor(loader, node):
	return collections.OrderedDict(loader.construct_pairs(node))
yaml.add_representer(collections.OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag,dict_constructor)

class CalcsAuditor:
	"""
	Read an audit file and print a report.
	"""
	_ignore_files = ['.gitignore']
	def __init__(self,**kwargs):
		self.root = 'calcs'
		self.ledger_fn = kwargs.pop('ledger','audit.yaml')
		self.debug = kwargs.pop('debug',False)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		status('welcome to the auditor')
		self.ledger = os.path.join(self.root,self.ledger_fn)
		if not os.path.isfile(self.ledger): raise Exception('cannot find %s'%self.ledger)
		else: 
			with open(self.ledger) as fp: self.raw = yaml.load(fp.read())
		# print everything
		asciitree(self.raw)
		self.interpret()
		if self.debug:
			import ipdb
			ipdb.set_trace()

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
					file_list = self.raw.get(kind,[])
					if type(file_list) in str_types: 
						raise Exception('from list "%s" is a string: %s'%(kind,file_list))
					# first check if the list names non-existent files
					missing_fns = [i for i in file_list if i not in self.files]
					if missing_fns:
						raise Exception(
							'the %s list has files that do not exist or have been ignored: %s'%(
								kind,missing_fns))
					self.files = [i for i in self.files if i not in file_list]
				elif details['do']=='account':
					file_list = self.raw.get(kind,[])
					if type(file_list)==dict: file_list = list(file_list.keys())
					# first check if the list names non-existent files
					missing_fns = [i for i in file_list if i not in self.files]
					if missing_fns:
						raise Exception(
							'the %s list has files that do not exist or have been ignored: %s'%(
								kind,missing_fns))
					self.files = [i for i in self.files if i not in file_list]
				else: raise Exception('do %s?'%details['do'])
			else: raise Exception('cannot interpret %s,%s'%(kind,details))
		if self.files:
			for f in sorted(self.files): print('  - %s'%f)
			raise Exception('unaccounted files are listed above')
