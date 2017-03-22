#!/usr/bin/python

from plotter import palette
	
def colorize(metadat,comparison='',resname=''):

	"""
	Master listing of colors for different PtdIns comparisons.
	"""
	
	if resname != '':
		colordict = {
			'DOPC':'grey',
			'DOPS':'red',
			'PIP2':'purple',}
		return palette.colors[colordict[resname]]
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
		colordict['symmetric'] = dict(colordict['ions'])
		colordict['asymmetric'] = dict(colordict['ions'])
		descriptor = tuple([metadat[i] for i in compare_to_value[comparison]])
		return palette.colors[colordict[comparison][descriptor]]

