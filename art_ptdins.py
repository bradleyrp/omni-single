#!/usr/bin/env python

"""
ART: PtdIns project
"""

#---forgo the use of the colors/labels since this project has complicated rendering
colors = None
labels = None

#---CUSTOM additions to globals
_art_words = ['colorize']

#---canonical colors for this project from brewer2mpl
import brewer2mpl

palette_colors = dict([(val,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[key]) for key,val in enumerate('red blue green purple orange yellow brown pink grey'.split())])
	
def colorize(metadat,comparison='',resname=''):
	"""
	Master listing of colors for different PtdIns comparisons.
	"""
	if resname != '':
		colordict = {
			'DOPC':'grey',
			'DOPS':'red',
			'PIP2':'purple',}
		return palette_colors[colordict[resname]]
	else:
		colordict = {
			'ions':{
				('NA',):'green',
				('MG',):'red',
				('Cal',):'blue',},
			'phosphate_position':{
				('NA','PI2P'):'green',
				('MG','PI2P'):'red',
				('Cal','PI2P'):'blue',
				('MG','P35P'):'purple',
				('Cal','P35P'):'orange',},
			'protonation':{
				('MG','PI2P'):'red',
				('Cal','PI2P'):'blue',
				('MG','P35P'):'purple',
				('Cal','P35P'):'orange',
				('NA','PI2P'):'green',
				('NA','PIPU'):'brown',
				('NA','PIPP'):'grey',},}
		compare_to_value = {
			'symmetric':['cation'],	
			'asymmetric':['cation'],
			'phosphate_position':['cation','ptdins_resname'],
			'protonation':['cation','ptdins_resname']}		
		#---aliases
		colordict['position'] = colordict['phosphate_position']
		compare_to_value['position'] = compare_to_value['phosphate_position']
		colordict['charge'] = colordict['protonation']
		compare_to_value['charge'] = compare_to_value['protonation']
		colordict['symmetric'] = dict(colordict['ions'])
		colordict['asymmetric'] = dict(colordict['ions'])
		descriptor = tuple([metadat[i] for i in compare_to_value[comparison]])
		return palette_colors[colordict[comparison][descriptor]]
