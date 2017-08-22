#!/usr/bin/env python

import mayavi
from mayavi import mlab
import numpy as np

def review3d(**kwargs):
	"""
	Review things in 3D. Originally developed for coiled coils.
	!Consolidate from structural_biology.py to a central location with a standard interface. 
	Note that scaling factors assume nanometers.
	Note that you may get VTK/TVTK version warnings that tell you to rebuild the classes.
	Note that you cannot drop ipdb.set_trace when vTK assumes the streams because of text colors.
	"""
	def pointplot(x,y,z,colormap='Spectral',r=1.0,color=(0,0,0)):
		return mlab.points3d(x,y,z,colormap=colormap,scale_factor=r,color=color)
	def lineplot(x,y,z,colormap='Spectral',r=0.1,color=(0,0,0)):
		mlab.plot3d(x,y,z,tube_radius=r,color=color,colormap=colormap)
	def pointplot_bare(x,y,z,colormap='Spectral',r=0.1): 
		return mlab.points3d(x,y,z,colormap=colormap,scale_factor=r)
	#---get incoming settings
	tube_thickness = kwargs.get('tube',0.1)
	sphere_radius = kwargs.get('radius',0.1)
	cmap = kwargs.get('cmap','Spectral')
	#---two modes: if you get colorset then we use the weird mayavi hack to get different colors/sizes
	if 'colorset' in kwargs:
		#---via: http://stackoverflow.com/questions/22253298/mayavi-points3d-with-different-size-and-colors
		for points in kwargs.get('points',[]): 
			pts = pointplot_bare(*points.T,colormap=cmap,r=sphere_radius)
			pts.glyph.scale_mode = 'scale_by_vector'
			#---set the size with the norm of a vector
			pts.mlab_source.dataset.point_data.vectors = np.tile(
				sphere_radius*np.ones((len(points),)),(3,1)).T
			#---set the color on the unit line
			pts.mlab_source.dataset.point_data.scalars = kwargs['colorset']
	#---standard execution
	else:
		for lines in kwargs.get('lines',[]): lineplot(*lines.T,r=tube_thickness,colormap=cmap)
		for points in kwargs.get('points',[]): pts = pointplot(*points.T,r=sphere_radius,colormap=cmap)
	if not kwargs.get('noshow',False): mlab.show()

def pbcbox(vecs):
	"""
	Plot a wireframe around the periodic boundary conditions assuming that they begin at the origin.
	"""
	mv = vecs
	mlab.plot3d(np.array([0.,mv[0],mv[0],0,0.,0.,mv[0],mv[0],0,0]),
		np.array([0.,0.,mv[1],mv[1],0,0,0,mv[1],mv[1],0]),
		np.array([0.,0.,0.,0.,0.,mv[2],mv[2],mv[2],mv[2],mv[2]]))
	mlab.plot3d(np.array([0,0,mv[0],mv[0],mv[0],mv[0]]),
		np.array([mv[1],mv[1],mv[1],mv[1],0,0]),np.array([0,mv[2],mv[2],0,0,mv[2]]))

