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
routine = ['spectra','height'][:]
sns = work.sns()

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

#---block: plot the undulation spectra
if 'spectra' in routine:

	art = {'fs':{'legend':12}}
	for wavevector_limit in [0.5,1.0,2.0]:
		status('plotting undulations with limit %f'%wavevector_limit,tag='plot')
		axes,fig = panelplot(layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},figsize=(5,5))
		undulation_panel(axes[0],data,work.meta,keys=np.array(sns),
			art=art,title='undulation spectra',lims=(0,wavevector_limit),colors=colors,labels=labels)
		meta = deepcopy(calc)
		meta['wavevector_limit'] = wavevector_limit
		descriptor = 'qlim.%.1f'%wavevector_limit
		picturesave('fig.%s.%s'%(plotname,descriptor),work.plotdir,backup=False,version=True,meta=meta)

#---block: plot the average heights
if 'height' in routine:

	#---! settings direct from pipeline-uvideos
	cmap_mpl_name = 'RdBu_r'
	fs = {'axlabel':14,'title':20,'legend':14,'ticks':14,'tiny':10,'text_note':12,}
	nrows,ncols = 1,len(work.sns())
	figsize = 8,8
	panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[nrows,ncols],
		'hspace':0.4}]},figsize=figsize)

	#---need number of proteins to draw hulls
	#---! more systematic way to define a parameter across all simulations, via standard_plots?
	try: nprots_list = [work.meta[sn]['nprots'] for sn in sns]
	except: nprots_list = [1 for sn in sns]
	#---! ditto
	try: titles = [work.meta[sn]['label'] for sn in sns]
	except: titles = [sn for sn in sns]

	import render
	from render.wavevids import print_birdseye_snapshot_render

	framecounts = [data[sn]['data']['nframes'] for sn in data]
	nframes = min(framecounts)
	data_protein,calc = plotload('protein_abstractor',work)
	protein_pts_all = np.array([data_protein[sn]['data']['points_all'] for sn in sns])
	protein_pts = [[protein_pts_all[ii][fr][...,:2] for ii,sn in enumerate(sns)] for fr in range(nframes)]
	mvecs = np.array([data[sn]['data']['vecs'][:nframes] for sn in sns]).mean(axis=1)
	avgzs = [data[sn]['data']['mesh'].mean(axis=0).mean(axis=0) for sn in sns]
	avgzs = [i-i.mean() for i in avgzs]
	extrema = np.array(avgzs).min(),np.array(avgzs).max()
	extrema = [np.abs(extrema).max()*j for j in [-1,1]]

	#---use the standard birdseye renderer to render the average structure
	print_birdseye_snapshot_render(
		avgzs,[p.mean(axis=0)[...,:2] for p in protein_pts_all],mvecs,nprots_list,
		handle='OUT',pbc_expand=1.0,smooth=1.,panelspecs=panelspecs,
		extrema=extrema,fs=fs,titles=titles,cmap_mpl_name=cmap_mpl_name,
		cbar_label_specs=dict(rotation=0,labelpad=-20,y=1.1),
		zlabel=r'$\mathrm{\langle z \rangle\,(nm)}$',
		**({} if is_live else dict(outdir=work.paths['post_plot_spot'],fn='fig.average_height')))
