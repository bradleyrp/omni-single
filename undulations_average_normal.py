#!/usr/bin/env python

from codes.mesh import measure_normal_deviation_from_wavy_surface

def undulations_average_normal(**kwargs):
	"""
	Compute normal fluctuations from an average surface with high excess area.
	"""
	#---parameters
	sn = kwargs['sn']
	work = kwargs['workspace']
	calc = kwargs['calc']
	upname = kwargs['upstream'].keys()
	if len(upname)!=1: 
		raise Exception('this calculation needs only a single upstream: regular grid types. received: %s'%upname)
	else: data = kwargs['upstream'][upname[0]]

	#---send the heights to the average normal calculation routine
	#---note all membrane data are stored with dimensions: leaflet, frame, X,Y even if they only have one leaflet
	surf = measure_normal_deviation_from_wavy_surface(data['mesh'].mean(axis=0),data['vecs'])

	#---pack
	result,attrs = {},{}
	result['average_normal_heights'] = surf
	return result,attrs	
