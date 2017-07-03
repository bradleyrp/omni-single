#!/usr/bin/env python

"""
Plot hydration shell statistics.
Requires standard plots.
"""

#---block: what to plot
do_strict = False
routine = ['basic','heatmap'][-1:]
sns = work.sns()
#---choosing a scanning range
hydration_scan_plot = np.sort(np.concatenate((np.arange(1.2,4.0,0.2),
	np.arange(4.0,30.0,4.0)),np.array([2.2,4.6])))
hydration_scan_plot = np.sort(np.concatenate((np.arange(1.2,4.0,0.2),
	np.arange(4.0,30.0,4.0))))

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

#---block: organize the data by cutoffs
if 'counts_by_zone_all' not in globals():

	#---warnings are errors
	if do_strict:
		import warnings
		warnings.filterwarnings('error')

	counts_by_zone_all = {}
	for sn in sns:
		dat = data[sn]['data']
		nframes = dat['nframes']
		near_lipids = dat['near_lipids']
		shell_counts = dat['shell_counts']
		frames_inds = np.arange(dat['nframes'])

		valid_frames = np.intersect1d(dat['valid_frames_shell_counts'],dat['valid_frames_near_lipids'])
		ind_map_shell_counts = np.in1d(dat['valid_frames_shell_counts'],valid_frames)
		ind_map_near_lipids = np.in1d(dat['valid_frames_near_lipids'],valid_frames)
		nframes_eff = len(valid_frames)

		nears = near_lipids[ind_map_near_lipids]
		counts = shell_counts[ind_map_shell_counts]

		#---cutoff sweep in angstroms
		lenscale = 10.0
		#---loop over pairs (wrote this hasty)
		pair_inds = (np.arange(0,2*len(hydration_scan))/2)[1:-1].reshape((-1,2))
		counts_by_zone = {}
		for ii,i in enumerate(pair_inds):
			status('correlating hydration by distance for %s'%sn,i=ii,looplen=len(pair_inds),tag='compute')
			lims = tuple([hydration_scan[j]/lenscale for j in i])
			#---loop over frames
			counts_by_frame = []
			#---! cannot vectorize this quickly. worth revisiting
			for fr in range(nframes_eff):
				counts_this = counts[fr][np.where(np.all((nears[fr]>=lims[0],nears[fr]<lims[1]),axis=0))[0]]
				counts_by_frame.extend(counts_this)
			counts_by_zone[lims] = np.array(counts_by_frame)
		counts_by_zone_all[sn] = counts_by_zone

#---block: plot the hydration versus distance
if 'basic' in routine:

	plotspecs = [{'show_errorbars':False},{'show_errorbars':True}]

	#---plot schemes
	for spec in plotspecs: 
		symbols = {}
		for sn in sns:
			if work.meta[sn]['ptdins_resname']=='P35P': symbols[sn] = 's'
			elif work.meta[sn]['cation']=="Na,Cal": symbols[sn] = 'x'
			else: symbols[sn] = 'o'
		figsize=(5,5)
		layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
		axes,fig = panelplot(layout,figsize=figsize)
		ax = axes[0]
		for sn in sns:
			counts_by_zone = counts_by_zone_all[sn]
			keys = np.array([key for key,val in counts_by_zone.items()])
			keys = keys[np.argsort(keys[:,0])]
			x_all = keys[:,1]
			y = np.array([counts_by_zone[tuple(key)].mean() for key in keys])
			y_err = np.array([counts_by_zone[tuple(key)].std() for key in keys])
			y_valid = np.where(~np.isnan(y))
			x,y,y_err = lenscale*keys[:,1][y_valid],y[y_valid],y_err[y_valid]
			x_spots = np.arange(len(x))
			color = colorize(work.meta[sn],comparison='protonation')
			if spec['show_errorbars']: ax.errorbar(x_spots+0.1*sns.index(sn),y,yerr=y_err,fmt='--o',color=color)
			else: ax.plot(x_spots,y,'-%s'%symbols[sn],ms=5,
				color=color,label='%s, %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label']))
		fn_fig = 'fig.hydration%s'%('.error_bars' if spec['show_errorbars'] else '')
		ax.set_ylabel('waters in the hydration shell')
		ax.set_xlabel('distance to lipids')
		ax.set_xticks(np.arange(keys[:,1].shape[0]))
		ax.set_xticklabels(keys[:,1],rotation=90)
		legend = ax.legend()
		picturesave(fn_fig,work.plotdir,backup=False,version=True,meta=spec)

#---block: plot heat maps of the hydration versus distance
if 'heatmap' in routine:

	#---reformulate the counts by zone
	reform = {}
	for sn in sns:
		counts_by_zone = counts_by_zone_all[sn]
		keys = np.array([key for key,val in counts_by_zone.items()])
		keys = keys[np.argsort(keys[:,0])]
		max_hydration = 10.0
		xvals = keys[:,1]
		yvals = np.arange(0,max_hydration)
		heat = np.zeros((len(keys),len(yvals)-1))
		for keynum,key in enumerate(keys):
			counts,bins = np.histogram(counts_by_zone[tuple(key)],bins=yvals)
			heat[keynum] = counts#/counts.sum().astype(float)
		heat[np.isnan(heat)] = 0
		reform[sn] = {'heat':heat,'r':xvals}
	maxz = max([r['heat'].max() for r in reform.values()])
	maxz = None

	#---filter distances
	if not all([np.all(reform[sn]['r']==reform[sns[0]]['r']) for sn in sns]):
		raise Exception('inconsistent cutoffs from hydration data')

	xvals_filter = np.array([0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,
		0.3,0.32,0.34,0.36,0.38,0.4,1.,1.4,1.8,2.2,2.6,2.8])
	distance_inds = np.where(np.in1d(xvals,xvals_filter))
	vlines = [0.4]

	figsize=(5,5)
	colormap = 'plasma'
	layout_name = 'summary2'
	aspect_equal = False
	layout = figlayout[layout_name]
	for ii,i in enumerate(layout['ins']): 
		layout['ins'][ii]['wspace'] = None
		layout['ins'][ii]['hspace'] = None

	axes,fig = panelplot(layout=layout,figsize=(12,16))
	for snum,sn in enumerate(sns):

		axrow,axcol = figplacer(sn,figplace[layout_name])
		ax = axes[axrow][axcol]
		heat,xvals = [reform[sn][i] for i in ['heat','r']]
		#---filter the distances
		heat = heat[distance_inds]
		xvals = xvals[distance_inds]
		kwargs = dict(extent=[0,len(xvals),0,yvals.max()+1],cmap=mpl.cm.__dict__[colormap])
		if maxz: kwargs.update(vmin=0,vmax=maxz)
		ax.imshow(heat.T,origin='lower',interpolation='nearest',**kwargs)
		ax.set_xticks(np.arange(len(xvals))+1)
		ax.set_xticklabels(xvals,rotation=90,fontsize=6)
		ax.set_yticks(yvals+0.5)
		ax.set_yticklabels(['%d'%i for i in yvals+1])
		if aspect_equal:
			x0,x1 = ax.get_xlim()
			y0,y1 = ax.get_ylim()
			ax.set_aspect((x1-x0)/(y1-y0))
		#---add vertical lines to denote changes in the step
		if not vlines:
			gaps = (xvals[1:]-xvals[:-1]).round(3)
			for xval in xvals[np.where((gaps[1:]-gaps[:-1])!=0.0)]: 
				x = np.where(xvals==xvals[np.where((gaps[1:]-gaps[:-1])!=0.0)])[0]+2
				ax.axvline(x,c='k',lw=1.5)
		else:
			for vline in vlines: ax.axvline(np.where(xvals==vline)[0]+1,c='k',lw=1.5)
		#---dislike tick marks
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.set_title('%s and %s%s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label'],
			' (phys)' if work.meta[sn]['composition_name']=='asymmetric' else ''))
	blank_unused_axes(axes,fig,figplace[layout_name])
	fn_fig = 'fig.hydration_review'
	picturesave(fn_fig,work.plotdir,backup=False,version=True,meta={})
