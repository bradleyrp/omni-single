#!/usr/bin/env python

"""
AUDITOR
Catalogs and provides the status for your calculations.
"""

import os,collections
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
