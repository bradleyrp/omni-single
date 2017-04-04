#!/usr/bin/env python

import mayavi
from mayavi import mlab

def review3d(**kwargs):
	"""
	Review things in 3D. Originally developed for coiled coils.
	!Consolidate from structural_biology.py to a central location with a standard interface. 
	Note that scaling factors assume nanometers.
	"""
	def pointplot(x,y,z,colormap='Spectral',r=1.0,color=(0,0,0)):
		mlab.points3d(x,y,z,colormap='Spectral',scale_factor=r,color=color)
	def lineplot(x,y,z,colormap='Spectral',r=0.1,color=(0,0,0)):
		mlab.plot3d(x,y,z,tube_radius=r,color=color)
	tube_thickness = kwargs.get('tube',0.1)
	sphere_radius = kwargs.get('radius',0.1)
	for lines in kwargs.get('lines',[]): lineplot(*lines.T,r=tube_thickness)
	for points in kwargs.get('points',[]): pointplot(*points.T,r=sphere_radius)
	mlab.show()
