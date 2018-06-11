#!/usr/bin/env python

"""
Scene file which returns video instructions from an interpreter function.
"""

import os

film_cuts = {
	'bilayer.side':{'does':['bonder'],'zoom':1.8,
		'resolution':(720,480),'resolution_snapshot':(720*2,480*2),'nframes_max':600,
		'kwargs':{'cgbondspath':os.path.join(os.getcwd(),'amx-martini/bin/cg_bonds.tcl',)},
		'selections':[dict(
			lipids=' or '.join(['resname %s'%i for i in work.vars['selectors']['resnames_lipid']]),
			style='Licorice 1.0 12.0 12.0',smooth=True,goodsell=True)]},}

def interpreter(): 
	global film_cuts
	return film_cuts
