#!/usr/bin/env python

"""
INSPECT THE CURVATURE COUPLING RESULTS
"""

import scipy
import scipy.interpolate
from omnicalc import store,load

deg2rad = lambda d: d/180.*np.pi
if is_live: 
	from ipywidgets import *
	from IPython.display import display

#---block: get the curvature investigation class
import codes.curvature_coupling.InvestigateCurvature
codes.curvature_coupling.InvestigateCurvature.plotload = plotload
codes.curvature_coupling.InvestigateCurvature.work = work
if 'ic' not in globals(): 
	ic = codes.curvature_coupling.InvestigateCurvature.InvestigateCurvature(mode='plot')

data = load('curvature_%s.dat'%ic.meta_spec['plot_cursor'],cwd=work.postdir)
sns = work.sns()

axes,fig = panelplot(figsize=(16,8),layout={'out':{'grid':[2,3],'hratios':[3,1],'hspace':0.0},
	'ins':[{'grid':[1,1]} for sn in sns]+[{'grid':[1,6]} for sn in sns]})
for snum,sn in enumerate(work.sns()):
	ax = axes[snum][0]
	field = data['%s_%s'%(sn,'cf')]
	vmax = max([np.abs(data['%s_%s'%(sn,'cf')]).max() for sn in work.sns()])
	ax.imshow(field.T,interpolation='nearest',
		origin='lower',vmin=vmax*-1,vmax=vmax,cmap=mpl.cm.__dict__['seismic'])
	for ax in axes[len(sns)+snum]:
		ax.imshow(field.T,interpolation='nearest',
			origin='lower',vmin=vmax*-1,vmax=vmax,cmap=mpl.cm.__dict__['seismic'])
picturesave('fig.curvature_optimized',
	work.plotdir,backup=False,version=True,meta={},extras=[])
