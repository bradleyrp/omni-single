#!/usr/bin/env python

"""
Plot hydration shell statistics.
Requires standard plots.
"""

#---block: what to plot
do_strict = False
routine = ['basic','heatmap'][:1]
sns = work.sns()
#---choosing a scanning range
roundlev = 4  

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)
	#---!
	atom_filter = calc['calcs']['specs']['atom_filter']
	distance_metric = calc['calcs']['specs']['distance_metric']

	#---plotspecs are formulated from key parameters in the metadata
	distance_ranges = work.plots.get('hydration',{}).get('specs',{}).get('distance_ranges_by_metric',{})
	#---previously we had multiple plotspecs but now they are set in metadata and we use one
	plotspecs = [dict([(['xmin','xmax','interval'][ii],i) 
		for ii,i in enumerate(distance_ranges[distance_metric])])]

	#---one plot per group of simulations set in the metadata
	sns_groups = work.plots.get('hydration',{}).get('specs',{}).get('sns',{})
	#---only load what we need
	sns = list(set([i for j in sns_groups.values() for i in j]))
	if not sns_groups: raise Exception('please set plots (hydration), specs, sns '
		'to set simulation groups for each plot')
	#---we discard the key in the sns groups set in the metadata
	plotspecs = [dict(sns=sns_groups[key],**spec) for spec in plotspecs for key in sns_groups]

	#---get all unique values in the ranges for post-processing
	hydration_scan = np.unique(np.concatenate([
		np.arange(plotspec['xmin'],plotspec['xmax']+plotspec['interval'],
			plotspec['interval']).round(roundlev) for plotspec in plotspecs]))
	#---use reduced keys
	key_pack = lambda *args: tuple([('%0.0'+'%d'%roundlev+'f')%float(k) for k in args])
	key_unpack = lambda arg: [float(i) for i in arg.split('_')]

#---block: remote work!
if 'postdat' in globals():

	#---unpack a transmitted postdat
	postdat_in = postdat
	postdat = {}
	for key,val in postdat_in.items():
		try: sn,lims = re.match('^(.*?)_(.*?)$',key).groups()
		except: 
			postdat[key] = val
			continue
		if sn not in postdat: postdat[sn] = {}
		#---reduced keys
		postdat[sn][key_pack(*key_unpack(lims))] = val

#---block: organize the data by cutoffs
if 'postdat' not in globals():

	#---warnings are errors
	if do_strict:
		import warnings
		warnings.filterwarnings('error')

	#---cutoff sweep in angstroms
	lenscale = 10.0
	postdat = {}
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

		#---loop over pairs (wrote this hasty)
		pair_inds = (np.arange(0,2*len(hydration_scan))/2)[1:-1].reshape((-1,2))
		dat = {}
		for ii,i in enumerate(pair_inds):
			status('correlating hydration by distance for %s'%sn,i=ii,looplen=len(pair_inds),tag='compute')
			lims = tuple([hydration_scan[j]/lenscale for j in i])
			#---loop over frames
			counts_by_frame = []
			#---! cannot vectorize this quickly. worth revisiting
			for fr in range(nframes_eff):
				counts_this = counts[fr][np.where(np.all((nears[fr]>=lims[0],nears[fr]<lims[1]),axis=0))[0]]
				counts_by_frame.extend(counts_this)
			dat[key_pack(*lims)] = np.array(counts_by_frame)
		postdat[sn] = dat

	#---make a copy for transmission
	#---! make this systematic
	reform = dict([('%s_%s_%s'%(sn,i[0],i[1]),j) 
		for sn in postdat for i,j in postdat[sn].items()])
	attrs = dict(sns=sns,meta=work.meta,vars=work.vars)
	if False: work.store(reform,'hydration.dat','~',attrs=attrs)
	
#---block: plot heat maps of the hydration versus distance
if 'heatmap' in routine:

	#---reformulate the counts by zone
	reform = {}
	for sn in sns:
		dat = postdat[sn]
		keys = np.array([key for key,val in dat.items()])
		keys = keys[np.argsort(keys[:,0])]
		max_hydration = 10.0
		xvals = keys[:,1]
		yvals = np.arange(0,max_hydration)
		heat = np.zeros((len(keys),len(yvals)-1))
		for keynum,key in enumerate(keys):
			counts,bins = np.histogram(dat[tuple(key)],bins=yvals)
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


#---block: plot the hydration versus distance
if 'basic' in routine:

	def stepper(x,y):
		"""Turn a plot into a step plot."""
		#---! decided not to use this function
		inds_y = np.arange(0,x.shape[0],0.5).astype(int)[1:-1]
		inds_x = np.arange(0,x.shape[0],0.5).astype(int)[:-2]
		return x[inds_x],y[inds_y]

	def slatplot(x,y,y_err,ax,color,w,thick_data,plus_off,basic_bar=False,**kwargs):
		"""Custom "slat" plots plot a line on top of an opaque confidence interval."""
		if basic_bar and thick_data!=None: raise Exception('cannot use basic_bar with thick_data')
		#---center or left markers
		align = kwargs.get('align','middle')
		#---independent variables are given by steps that reach out from the midpoint according to width (w)
		if align=='middle': xs = [x-w/2.,x+w/2.]
		elif align=='left': xs = [x,x+w]
		else: raise Exception('align must be middle or left')
		#---check overlapping points from other timeseries
		global point_registry
		#---rounding tolerance level for identifying overlapping points
		roundlev = kwargs.get('roundlev',4)
		#---the reduced point is used for checking for overlaps with a tolerance given by roundlev
		pts_red = [np.round(x,roundlev),np.round(y,roundlev)]
		#---if points are (even nearly) overlapping, we offset the height slightly
		if thick_data!=None and pts_red in point_registry: 
			#---note that we never offset without thick_data
			off = 2*thick_data*sum([i==pts_red for i in point_registry ])
		else: off = 0.0
		if thick_data==None:
			ax.plot(xs,[y+off,y+off],solid_capstyle='butt',lw=2,color=color)
			if basic_bar:
				ax.fill_between(xs,[0,0],[y+off,y+off],color=color,alpha=0.35,lw=0,zorder=2)
		else:
			#---the thicknesss of the step lines in the y-direction is given in data units
			yt = thick_data
			#---plot the horizontal step lines
			ax.fill_between(xs,y+off-yt,y+off+yt,color=color,zorder=4,lw=0)
		#---remember this point for checking overlaps later
		point_registry.append([np.round(x,roundlev),np.round(y,roundlev)])
		#---for hydration, there may be some points with no variance
		if y_err==0.0 and not y_err==None:
			#---! previously plotted these points via: ax.scatter(x,y,marker='+',color=color,zorder=3)
			#---make a custom plus sign with a height equal to thrice the height of the steps
			#---...and a width which is inset from the edges of the step by the plus_off flag
			#---...which we recommend setting as some proportion of the width of the bin/interval
			#---note that we only plot the special signs if there is no confidence interval i.e. if y_err is
			#---...zero because otherwise the viewer would identify the overlap by the color, which is 
			#---...darkened by layering opaque panels (note that more than two colors makes this ugly though)
			ax.fill_between(
				[xs[0]+plus_off,xs[1]-plus_off],
				[y+off-yt*3,y+off-yt*3],
				[y+off+yt*3,y+off+yt*3],
				#---to debug the plus offsets, try adding "color if not off else 'k'"
				color=color,zorder=4,lw=0)
		#---make the confidence interval plot
		elif y_err!=None:
			ys1,ys2 = [y-y_err for j in range(2)],[y+y_err for j in range(2)]
			ax.fill_between(xs,ys1,ys2,color=color,alpha=0.35,lw=0,zorder=2)
			#---lay down some white so that gridlines cannot be seen behind the opaque panels
			ax.fill_between(xs,ys1,ys2,color='w',alpha=1.0,lw=0,zorder=1)
	
	figspecs = {'fs_ticks':16,'fs_ylabel':20,'fs_xlabel':20,'fs_legend':14}

	#---plot schemes
	for spec in plotspecs:
		symbols = {}
		sns_this = [sn for sn in work.sns() if sn in spec['sns']]
		for sn in sns_this:
			if work.meta[sn]['ptdins_resname']=='P35P': symbols[sn] = 's'
			elif work.meta[sn]['cation']=="Na,Cal": symbols[sn] = 'x'
			else: symbols[sn] = 'o'
		#---track overlapping points
		point_registry = []
		legendspec = []
		figsize=(10,10)
		layout = {'out':{'grid':[1,1]},'ins':[{'grid':[2,1],'hratios':[1,4],'hspace':0.02}]}
		axes,fig = panelplot(layout,figsize=figsize)
		#---the main plot lies below
		ax,axtop = axes[1],axes[0]
		x_collects,maxcount = [],0

		do_alt_colors = len(set([work.meta[s]['cation'] for s in sns_this]))==1

		for sn in sns_this:
			dat = postdat[sn]
			hydration_scan_plot = np.arange(plotspec['xmin'],plotspec['xmax']+
				plotspec['interval'],plotspec['interval']).round(roundlev)
			#---reference x values in angstroms
			x = hydration_scan_plot[1:]
			x_pairs = np.array([(hydration_scan_plot[i:i+2]/lenscale) 
				for i in range(len(hydration_scan_plot)-1)]).round(roundlev)
			y = np.array([dat[key_pack(*xp)].mean() 
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs])
			#---! repetitive
			yn = np.array([dat[key_pack(*xp)].shape[0]/float(data[sn]['data']['nframes'])
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs])
			#---! handling alternate compositions here
			if work.meta[sn]['composition_name']=='symmetric': yn = yn/2.0
			maxcount = max([yn.max(),maxcount])
			y_err = np.array([dat[key_pack(*xp)].std() 
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs])
			valid = np.where(y!=-1.0)[0]
			#---use alternate colors if the ions are the same
			#---! currently hard-coded
			if do_alt_colors:
				color = {'membrane-v509':'r','membrane-v530':'b'}[sn]
			else: color = colorize(work.meta[sn],comparison='protonation')
			kwargs = dict(ms=5,color=color,
				label='%s, %s%s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'],
					'(%s)'%work.meta[sn]['composition_name'] if do_alt_colors else ''))
			#---for each valid point we plot a "slat" showing the confidence interval or a plus sign
			for v in valid:
				slatplot(x[v],y[v],y_err=y_err[v],ax=ax,align='left',
					#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
					w=spec['interval'],color=color,
					#---set the linewidth of the steps in data units, a (small) fraction of a water
					thick_data=0.02,
					#---set the inset trim of the zero-variance plus sign
					plus_off=spec['interval']/3.)
			x_collects.extend(x[valid])
			#---ticks follow th left alignment
			minor_xticks = np.unique(np.reshape([[i,i-spec['interval']] for i in x],-1))
			#---plot the normalization constants above ...!!!
			for v in valid:
				slatplot(x[v],yn[v],y_err=None,ax=axtop,align='left',
					#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
					w=spec['interval'],color=color,thick_data=None,plus_off=None,basic_bar=True)
			#---identical tick-mark formatting
			for ax_this in [ax,axtop]:
				ax_this.set_xticks(minor_xticks,minor=True)
				ax_this.set_axisbelow(True)
				ax_this.xaxis.grid(True,which='minor',zorder=0)
			#---major ticks are set automatically and just enlarged a bit here
			ax.set_yticklabels(['%d'%i for i in ax.get_yticks()],
				fontsize=figspecs['fs_ticks'])
			axtop.set_xticklabels([])
			ax.yaxis.grid(True,which='major',zorder=0)
			#---! should we match the opacity with slatplot?
			legendspec.append(dict(name='%s\n%s%s'%(
				work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'],
				'\n(%s)'%work.meta[sn]['composition_name'] if do_alt_colors else ''),
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color,alpha=0.5)))
		axtop.set_ylim((0,maxcount*1.1))
		#---dislike tick marks
		axtop.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		axtop.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		for ax_this in [ax,axtop]:
			ax_this.set_xlim(np.array(x_collects).min()-spec['interval'],np.array(x_collects).max())
		#---major ticks are set automatically and just enlarged a bit here
		ax.set_xticklabels(['%.1f'%i for i in ax.get_xticks()],
			fontsize=figspecs['fs_ticks'])
		fn_fig = 'fig.hydration'
		axtop.set_ylabel(r'$\mathrm{N_{ions}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_ylabel(r'$\mathrm{N_{waters}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_xlabel('cation-lipid distance ($\mathrm{\AA}$)',fontsize=figspecs['fs_xlabel'])
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		patch = mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)
		patches.extend([
			mpl.lines.Line2D([],[],color='k',lw=4),
			mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)])
		labels.extend([r'$\mathrm{\langle N_{waters} \rangle}$','zero\nvariance'])
		legend = ax.legend(patches,labels,loc='upper left',fontsize=figspecs['fs_legend'],
			#---no more than 8 items in a column
			bbox_to_anchor=(1.05,0.0,1.,1.),labelspacing=1.2,ncol=int(np.ceil(len(sns_this)/8.)),
			handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
		#---we must pass-through
		#---! note this in the lab notebook!
		picturesave(fn_fig,work.plotdir,backup=False,version=True,
			meta=dict(atom_filter=atom_filter,distance_metric=distance_metric,**spec),extras=[legend])
