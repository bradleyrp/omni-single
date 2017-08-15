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
routine_all = ['spectra','height','dump_xyz','height_histograms']
routine = work.plots.get(plotname,{}).get('specs',{}).get('routine',routine_all)
sns = work.sns()

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

#---block: plot the undulation spectra
if 'spectra' in routine:

	#---settings
	art = {'fs':{'legend':12}}
	layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
	figsize=(5,5)
	wavevector_limits = [0.5,1.0,2.0]

	#---sweep over plots
	plotspecs = sweeper(**{'layout':layout,'wavevector_limit':wavevector_limits})
	for pnum,plotspec in enumerate(plotspecs):
		status('plotting layout',i=pnum,tag='plot')
		wavevector_limit = plotspec['wavevector_limit']
		axes,fig = panelplot(layout,figsize=figsize)
		for aa,ax in enumerate(axes):
			#---panel settings
			title = 'undulation spectra'
			keys = None
			#---plotter
			undulation_panel(ax,data,
				#---! note that labels,colors comes from art
				keys=keys,art=art,title=title,labels=labels,colors=colors,
				lims=(0,wavevector_limit))
		#---metadata and file tag
		meta = deepcopy(calc)
		meta.update(**plotspec)
		tag = 'qlim.%.1f'%wavevector_limit
		picturesave('fig.%s.%s'%(plotname,tag),work.plotdir,backup=False,version=True,meta=plotspec)

#---block: plot the average heights
if 'height' in routine:

	#---! settings direct from pipeline-uvideos
	cmap_mpl_name = 'RdBu_r'
	fs = {'axlabel':14,'title':20,'legend':14,'ticks':14,'tiny':10,'text_note':12,}
	fs = dict([(i,j-4) for i,j in fs.items()])
	nrows,ncols = 1,len(work.sns())
	figsize = 8,8
	panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[nrows,ncols],
		'hspace':0.4,'wspace':0.4}]},figsize=figsize)

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
	try: data_protein,_ = plotload('protein_abstractor',work)
	except: data_protein = None
	if data_protein:
		protein_pts_all = np.array([data_protein[sn]['data']['points_all'] for sn in sns])
		protein_pts = [[protein_pts_all[ii][fr][...,:2] for ii,sn in enumerate(sns)] for fr in range(nframes)]
	mvecs = np.array([data[sn]['data']['vecs'][:nframes] for sn in sns]).mean(axis=1)
	avgzs = [data[sn]['data']['mesh'].mean(axis=0).mean(axis=0) for sn in sns]
	avgzs = [i-i.mean() for i in avgzs]
	extrema = np.array(avgzs).min(),np.array(avgzs).max()
	extrema = [np.abs(extrema).max()*j for j in [-1,1]]

	#---use the standard birdseye renderer to render the average structure
	print_birdseye_snapshot_render(
		avgzs,None if not data_protein else [p.mean(axis=0)[...,:2] for p in protein_pts_all],
		mvecs,nprots_list,handle='OUT',pbc_expand=1.0,smooth=1.,panelspecs=panelspecs,
		extrema=extrema,fs=fs,titles=titles,cmap_mpl_name=cmap_mpl_name,
		cbar_label_specs=dict(rotation=0,labelpad=-20,y=1.1),
		zlabel=r'$\mathrm{\langle z \rangle\,(nm)}$',
		**({} if is_live else dict(outdir=work.paths['post_plot_spot'],fn='fig.average_height')))

#---block: export data for comparison to other labmates' methods
if 'dump_xyz' in routine:

	sn = 'simulation-v003-IRSp53-large-bilayer-50-50-DOPS-DOPC'
	outname = 'banana-v003-0-400000-2000'
	dat = data[sn]['data']
	mesh = dat['mesh']	
	nframes = min(data[sn]['data']['mesh'].shape[1],500)
	out = mesh[:,:nframes].reshape((-1,49*49))
	backin = out.reshape((2,-1,49,49))
	print('checked? %s'%np.all(backin==mesh[:,:nframes]))
	np.savetxt(
		work.postdir+'%s-%d-frames.monolayer_heights_49x49.txt'%(outname,nframes),
		out)
	sys.path.insert(0,'calcs')
	from codes.undulate_plot import calculate_undulations
	surf = backin.mean(axis=0)
	vecs = dat['vecs']	
	lims = (0,1.0) 
	uspec = calculate_undulations(surf,vecs,
		chop_last=True,perfect=True,lims=lims,raw=False)
	np.savetxt(
		work.postdir+'banana-v003-0-400000-2000-%d-frames.hqhq.txt'%nframes,
		np.array(np.transpose([uspec['x'],uspec['y']])))

	#---the above uses the undulation data but we also load some of the raw points and dump them
	#---...the following section uses lipid_abstractor, collections all, slices current_protein
	data,calc = plotload('lipid_abstractor',work)
	nframes = min(data[sn]['data']['vecs'].shape[0],nframes)
	np.savetxt(
		work.postdir+'%s-%s-frames.xyz.txt'%(outname,nframes),
		data[sn]['data']['points'][:nframes].reshape((-1,3)))
	np.savetxt(
		work.postdir+'%s-%s-frames.vecs.txt'%(outname,nframes),
		data[sn]['data']['vecs'][:nframes])
	np.savetxt(
		work.postdir+'%s.monolayers.txt'%outname,
		data[sn]['data']['monolayer_indices'])
	#---checking reshape 
	out = data[sn]['data']['points'][:nframes].reshape((-1,3))
	backin = out.reshape((nframes,-1,3))
	dat = data[sn]['data']['points'][:nframes]
	print('checked? %s'%np.all(dat[0]==backin[0]))

#---block: plot the average heights
if 'height' in routine:

	#---! settings direct from pipeline-uvideos
	cmap_mpl_name = 'RdBu_r'
	fs = {'axlabel':14,'title':20,'legend':14,'ticks':14,'tiny':10,'text_note':12,}
	fs = dict([(i,j-4) for i,j in fs.items()])
	nrows,ncols = 1,len(work.sns())
	figsize = 8,8
	panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[nrows,ncols],
		'hspace':0.4,'wspace':0.4}]},figsize=figsize)

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
	data_protein,_ = plotload('protein_abstractor',work)
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

#---block: plot the average heights
if 'height_histograms' in routine:

	postdat = {}
	zs = dict([(sn,data[sn]['data']['mesh'].mean(axis=0)) for sn in sns])
	for snum,sn in enumerate(sns):
		heights = (zs[sn]-zs[sn].mean()).reshape(-1)
		postdat[sn] = heights

	panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1],
		'hspace':0.4,'wspace':0.4}]},figsize=(8,6))
	axes,fig = panelplot(**panelspecs)
	zmax = max([np.abs(postdat[sn]).max() for sn in sns])
	halfrange = np.linspace(0,zmax*1.1,100)
	hbins = np.unique(np.concatenate((halfrange,-1*halfrange)))
	ax = axes[0]
	for sn in sns:
		counts,bins = np.histogram(postdat[sn],bins=hbins)
		ax.plot((hbins[1:]+hbins[:-1])/2.,counts,label=sn)
	legend = ax.legend()
	ax.axvline(0.0,c='k')
	picturesave('fig.undulations.height_histograms',work.plotdir,
		backup=False,version=True,meta={},extras=[legend])
