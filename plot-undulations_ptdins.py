#!/usr/bin/env python

"""
Plot undulation spectra.
Requires standard plots.
"""

from codes import undulate
from copy import deepcopy
from codes.undulate_plot import undulation_panel
import numpy as np

#---block: what to plot
routine = ['spectra','height'][:1]
sns = work.sns()

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)
	#---reset the plotname (if we are aliased)
	plotname = 'undulations'

#---block: plot the undulation spectra
if 'spectra' in routine:
	
	#---settings
	art = {'fs':{'legend':12}}
	layouts = {'4x4':{'out':{'grid':[1,1]},'ins':[{'grid':[2,2]}]}}
	layout = layouts['4x4']
	figsize = (12,12)
	wavevector_limits = [0.5,1.0,2.0]
	collections = ['symmetric','asymmetric','position','charge']

	#---sweep over plots
	plotspecs = sweeper(**{'layout':layouts,'wavevector_limit':wavevector_limits})
	for pnum,plotspec in enumerate(plotspecs):
		status('plotting layout',i=pnum,tag='plot')
		wavevector_limit = plotspec['wavevector_limit']
		axes,fig = panelplot(layout,figsize=figsize)
		for aa,ax in enumerate(axes):
			#---panel settings
			title = work.vars['composition_names'][collections[aa]]+', '+collections[aa]
			keys = work.specs['collections'][collections[aa]]
			labels = dict([(sn,work.meta[sn]['ptdins_label']+','+work.meta[sn]['ion_label']) for sn in keys])
			colors = dict([(sn,colorize(work.meta[sn],comparison=collections[aa])) for sn in keys])
			undulation_panel(ax,data,
				keys=keys,art=art,title=title,labels=labels,colors=colors,
				lims=(0,wavevector_limit))
		#---metadata and file tag
		meta = deepcopy(calc)
		meta.update(**plotspec)
		tag = 'qlim.%.1f'%wavevector_limit
		picturesave('fig.%s.%s'%(plotname,tag),work.plotdir,backup=False,version=True,meta=plotspec)
