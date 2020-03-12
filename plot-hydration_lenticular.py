#!/usr/bin/env python

"""
Plot hydration shell statistics.
Requires standard plots.
"""

import scipy
import scipy.interpolate

#---block: what to plot
do_strict = False
routine = ['basic','heatmap','basic_combo','ultimate'][-1:]
sns = work.sns()
#---choosing a scanning range
roundlev = 4  

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname)

	#---! updated slightly for new plotload
	atom_filter = calc['hydration']['calcs']['specs']['atom_filter']
	distance_metric = calc['hydration']['calcs']['specs']['distance_metric']

	#---plotspecs are formulated from key parameters in the metadata
	distance_ranges = work.plots.get('hydration_lenticular',{}).get(
		'specs',{}).get('distance_ranges_by_metric',{})
	#---previously we had multiple plotspecs but now they are set in metadata and we use one
	plotspecs = [dict([(['xmin','xmax','interval'][ii],i) 
		for ii,i in enumerate(distance_ranges[distance_metric])])]

	#---one plot per group of simulations set in the metadata
	sns_groups = work.plots.get('hydration_lenticular',{}).get('specs',{}).get('sns',{})
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

	def slatplot(x,y,y_err,ax,color,w,thick_data,plus_off,plus_show=True,
		basic_bar=False,bar_gradient=False,distn=None,lw_mean=2,**kwargs):
		"""Custom "slat" plots plot a line on top of an opaque confidence interval."""
		mix_black = list(np.mean([mpl.colors.to_rgb('k'),color],axis=0))
		step_color = color if not bar_gradient else mix_black
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
			ax.plot(xs,[y+off,y+off],solid_capstyle='butt',lw=lw_mean,color=step_color)
			if basic_bar:
				ax.fill_between(xs,[0,0],[y+off,y+off],color=step_color,alpha=0.35,lw=0,zorder=2)
		else:
			#---the thicknesss of the step lines in the y-direction is given in data units
			yt = thick_data
			#---plot the horizontal step lines
			ax.fill_between(xs,y+off-yt,y+off+yt,color=step_color,zorder=4,lw=0)
		#---remember this point for checking overlaps later
		point_registry.append([np.round(x,roundlev),np.round(y,roundlev)])
		#---for hydration, there may be some points with no variance
		if plus_show and y_err==0.0 and not y_err==None:
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
			#---! tuple is for box plot
			if type(y_err)==tuple:
				ys1,ys2 = [y_err[0] for j in range(2)],[y_err[1] for j in range(2)]
			else:
				ys1,ys2 = [y-y_err for j in range(2)],[y+y_err for j in range(2)]
			if not bar_gradient:
				ax.fill_between(xs,ys1,ys2,color=color,alpha=0.35,lw=0,zorder=2)
				#---lay down some white so that gridlines cannot be seen behind the opaque panels
				ax.fill_between(xs,ys1,ys2,color='w',alpha=1.0,lw=0,zorder=1)
			#---fancy gradient method
			else:
				#---lay down some white so that gridlines cannot be seen behind the opaque panels
				ax.fill_between(xs,ys1,ys2,color='w',alpha=1.0,lw=0,zorder=1)
				#---make the density gradient bars
				counts,bins = np.histogram(distn,bins=np.linspace(distn.min(),distn.max(),101))
				mids = (bins[1:]+bins[:-1])/2.
				mids[counts>0],counts[counts>0]
				try: 
					interp = scipy.interpolate.griddata(mids[counts>0],counts[counts>0],mids)
					gradient_data = np.array([interp])
				except: gradient_data = np.array([counts])
				_,cmap = colorscale(bands=[list(mpl.colors.to_rgb(color))+[0.07],color],return_cmap=True)
				ax.imshow(gradient_data.T,origin='lower',interpolation='bicubic',
					cmap=cmap,extent=(xs[0],xs[1],distn.min(),distn.max()),alpha=1,zorder=2)
	
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
		figsize=(18,10)
		layout = {'out':{'grid':[1,1]},'ins':[{'grid':[2,1],'hratios':[1,1],'hspace':0.02}]}
		axes,fig = panelplot(layout,figsize=figsize)
		#---the main plot lies below
		ax,axtop = axes[1],axes[0]
		x_collects,maxcount = [],0

		do_alt_colors = len(set([work.meta[s]['cation'] for s in sns_this]))==1

		if len(sns_this)!=2: raise Exception('lenticular only works with pairs!')
		interval_lenticular = spec['interval']/float(len(sns_this))
		#---! spec for the standard deviation vs box plots
		lenticular_bar_style = ['std','boxplot','dev'][-1]

		for snum,sn in enumerate(sns_this):

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
			distns = None
			if lenticular_bar_style=='boxplot':
				y_err = [(np.percentile(dat[key_pack(*xp)],25),np.percentile(dat[key_pack(*xp)],75))
					if len(dat[key_pack(*xp)])>0 
					else -1 for xp in x_pairs]
			elif lenticular_bar_style=='std':
				y_err = np.array([dat[key_pack(*xp)].std() 
					if len(dat[key_pack(*xp)])>0 
					else -1 for xp in x_pairs])
			elif lenticular_bar_style=='dev':
				y_err = [(dat[key_pack(*xp)].min(),dat[key_pack(*xp)].max())
					if len(dat[key_pack(*xp)])>0 
					else -1 for xp in x_pairs]
				distns = [dat[key_pack(*xp)] if len(dat[key_pack(*xp)])>0 else -1 for xp in x_pairs]
			else: raise Exception('unclear lenticular bar style %s'%lenticular_bar_style)
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
				gradient_details = {}
				if type(distns)!=type(None): gradient_details.update(distn=distns[v],bar_gradient=True)
				slatplot(x[v]+interval_lenticular*snum,y[v],y_err=y_err[v],ax=ax,align='left',
					#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
					w=interval_lenticular,color=color,
					#---set the linewidth of the steps in data units, a (small) fraction of a water
					thick_data=0.05,
					#---set the inset trim of the zero-variance plus sign
					plus_off=spec['interval']/3.,plus_show=False,**gradient_details)
			x_collects.extend(x[valid])
			#---ticks follow th left alignment
			minor_xticks = np.unique(np.reshape([[i,i-spec['interval']] for i in x],-1))
			#---plot the normalization constants above ...!!!
			for v in valid:
				slatplot(x[v]+interval_lenticular*snum,yn[v],y_err=None,ax=axtop,align='left',
					#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
					w=interval_lenticular,color=color,thick_data=None,plus_off=None,lw_mean=0,basic_bar=True)
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
		if lenticular_bar_style=='dev': 
			ax.set_aspect('auto')
			ax.set_ylim((0,8))
		axtop.set_ylim((0,maxcount*1.1))
		#---dislike tick marks
		axtop.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		axtop.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		for ax_this in [ax,axtop]:
			ax_this.set_xlim(np.array(x_collects).min()-spec['interval'],np.array(x_collects).max())
		#---major ticks are set automatically and just enlarged a bit here
		ax.set_xticklabels(['%.1f'%i for i in ax.get_xticks()],
			fontsize=figspecs['fs_ticks'])
		fn_fig = 'fig.hydration_lenticular'
		axtop.set_ylabel(r'$\mathrm{N_{ions}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_ylabel(r'$\mathrm{N_{waters}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_xlabel('cation-lipid distance ($\mathrm{\AA}$)',fontsize=figspecs['fs_xlabel'])
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		patch = mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)
		patches.extend([
			mpl.lines.Line2D([],[],color='k',lw=4),
			mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)])
		if False: labels.extend([r'$\mathrm{\langle N_{waters} \rangle}$','zero\nvariance'])
		legend = ax.legend(patches,labels,loc='upper left',fontsize=figspecs['fs_legend'],
			#---no more than 8 items in a column
			bbox_to_anchor=(1.05,0.0,1.,1.),labelspacing=1.2,ncol=int(np.ceil(len(sns_this)/8.)),
			handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
		#---we must pass-through
		#---! note this in the lab notebook!
		picturesave(fn_fig,work.plotdir,backup=False,version=True,
			meta=dict(atom_filter=atom_filter,distance_metric=distance_metric,**spec),extras=[legend])

#---block: plot the hydration versus distance
#! this section is very repetitive with the above
if 'basic_combo' in routine:

	def stepper(x,y):
		"""Turn a plot into a step plot."""
		#---! decided not to use this function
		inds_y = np.arange(0,x.shape[0],0.5).astype(int)[1:-1]
		inds_x = np.arange(0,x.shape[0],0.5).astype(int)[:-2]
		return x[inds_x],y[inds_y]

	def slatplot(x,y,y_err,ax,color,w,thick_data,plus_off,plus_show=True,
		basic_bar=False,bar_gradient=False,distn=None,alpha_basic=0.35,lw_mean=2,**kwargs):
		"""Custom "slat" plots plot a line on top of an opaque confidence interval."""
		mix_black = list(np.mean([mpl.colors.to_rgb('k'),color],axis=0))
		step_color = color if not bar_gradient else mix_black
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
			ax.plot(xs,[y+off,y+off],solid_capstyle='butt',lw=lw_mean,color=step_color)
			if basic_bar:
				ax.fill_between(xs,[0,0],[y+off,y+off],color=step_color,alpha=alpha_basic,lw=0,zorder=2)
		else:
			#---the thicknesss of the step lines in the y-direction is given in data units
			yt = thick_data
			#---plot the horizontal step lines
			ax.fill_between(xs,y+off-yt,y+off+yt,color=step_color,zorder=4,lw=2,edgecolor='w')
		#---remember this point for checking overlaps later
		point_registry.append([np.round(x,roundlev),np.round(y,roundlev)])
		#---for hydration, there may be some points with no variance
		if plus_show and y_err==0.0 and not y_err==None:
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
			#---! tuple is for box plot
			if type(y_err)==tuple:
				ys1,ys2 = [y_err[0] for j in range(2)],[y_err[1] for j in range(2)]
			else:
				ys1,ys2 = [y-y_err for j in range(2)],[y+y_err for j in range(2)]
			if not bar_gradient:
				ax.fill_between(xs,ys1,ys2,color=color,alpha=0.35,lw=0,zorder=2)
				#---lay down some white so that gridlines cannot be seen behind the opaque panels
				ax.fill_between(xs,ys1,ys2,color='w',alpha=1.0,lw=0,zorder=1)
			#---fancy gradient method
			else:
				#---lay down some white so that gridlines cannot be seen behind the opaque panels
				ax.fill_between(xs,ys1,ys2,color='w',alpha=1.0,lw=0,zorder=1)
				#---make the density gradient bars
				counts,bins = np.histogram(distn,bins=np.linspace(distn.min(),distn.max(),101))
				mids = (bins[1:]+bins[:-1])/2.
				mids[counts>0],counts[counts>0]
				try: 
					interp = scipy.interpolate.griddata(mids[counts>0],counts[counts>0],mids)
					gradient_data = np.array([interp])
				except: gradient_data = np.array([counts])
				_,cmap = colorscale(bands=[list(mpl.colors.to_rgb(color))+[0.07],color],return_cmap=True)
				ax.imshow(gradient_data.T,origin='lower',interpolation='bicubic',
					cmap=cmap,extent=(xs[0],xs[1],distn.min(),distn.max()),alpha=1,zorder=2)
	
	figspecs = {'fs_ticks':16,'fs_ylabel':20,'fs_xlabel':20,'fs_legend':14}

	#---! moved and modified for combo
	figsize=(18,10)
	layout = {'out':{'grid':[1,2],'wspace':0.15},
		'ins':[{'grid':[2,1],'hratios':[1,4],'hspace':0.02} for i in range(2)]}
	axes,fig = panelplot(layout,figsize=figsize)

	symbols = {}
	#! error on revisiting this script. changed spec to plotspec below
	sns_this = [sn for sn in work.sns() if sn in plotspec['sns']]
	for sn in sns_this:
		if work.meta[sn]['ptdins_resname']=='P35P': symbols[sn] = 's'
		elif work.meta[sn]['cation']=="Na,Cal": symbols[sn] = 'x'
		else: symbols[sn] = 'o'
	#---track overlapping points
	point_registry = []
	legendspec = []
	x_collects,maxcount = [],0

	do_alt_colors = len(set([work.meta[s]['cation'] for s in sns_this]))==1

	if len(sns_this)!=2: raise Exception('lenticular only works with pairs!')
	interval_lenticular = spec['interval']/float(len(sns_this))
	#---! spec for the standard deviation vs box plots
	lenticular_bar_style = ['std','boxplot','dev'][-1]

	axinds = {'membrane-v531':0,'membrane-v532':0,'membrane-v533':1,'membrane-v534':1}
	snum_offs = {'membrane-v531':0,'membrane-v532':1,'membrane-v533':0,'membrane-v534':1}
	sns_this = ['membrane-v531','membrane-v532','membrane-v533','membrane-v534',]

	for snum,sn in enumerate(sns_this):

		#---moved and modified for combo
		#---the main plot lies below
		ax,axtop = axes[axinds[sn]][1],axes[axinds[sn]][0]

		snum_off = snum_offs[sn]

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
		distns = None
		if lenticular_bar_style=='boxplot':
			y_err = [(np.percentile(dat[key_pack(*xp)],25),np.percentile(dat[key_pack(*xp)],75))
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs]
		elif lenticular_bar_style=='std':
			y_err = np.array([dat[key_pack(*xp)].std() 
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs])
		elif lenticular_bar_style=='dev':
			y_err = [(dat[key_pack(*xp)].min(),dat[key_pack(*xp)].max())
				if len(dat[key_pack(*xp)])>0 
				else -1 for xp in x_pairs]
			distns = [dat[key_pack(*xp)] if len(dat[key_pack(*xp)])>0 else -1 for xp in x_pairs]
		else: raise Exception('unclear lenticular bar style %s'%lenticular_bar_style)
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
			gradient_details = {}
			if type(distns)!=type(None): gradient_details.update(distn=distns[v],bar_gradient=True)
			slatplot(x[v]+interval_lenticular*snum_off,y[v],y_err=y_err[v],ax=ax,align='left',
				#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
				w=interval_lenticular,color=color,
				#---set the linewidth of the steps in data units, a (small) fraction of a water
				thick_data=0.05,
				#---set the inset trim of the zero-variance plus sign
				plus_off=spec['interval']/3.,plus_show=False,**gradient_details)
		x_collects.extend(x[valid])
		#---ticks follow th left alignment
		minor_xticks = np.unique(np.reshape([[i,i-spec['interval']] for i in x],-1))
		#---plot the normalization constants above ...!!!
		for v in valid:
			slatplot(x[v]+interval_lenticular*snum_off,yn[v],y_err=None,ax=axtop,align='left',
				#---the width of the bin is the sampling interval in the lipid-ion cutoff distance
				w=interval_lenticular,color=color,alpha_basic=0.9,
				thick_data=None,plus_off=None,lw_mean=0,basic_bar=True)
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
		if lenticular_bar_style=='dev': 
			ax.set_aspect('auto')
			ax.set_ylim((0,8))
		axtop.set_ylim((0,maxcount*1.1))
		#---dislike tick marks
		axtop.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		axtop.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		for ax_this in [ax,axtop]:
			ax_this.set_xlim(np.array(x_collects).min()-spec['interval'],np.array(x_collects).max())
		#---major ticks are set automatically and just enlarged a bit here
		ax.set_xticklabels(['%.1f'%i for i in ax.get_xticks()],
			fontsize=figspecs['fs_ticks'])
		axtop.set_ylabel(r'$\mathrm{N_{ions}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_ylabel(r'$\mathrm{N_{waters}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_xlabel('cation-lipid distance ($\mathrm{\AA}$)',fontsize=figspecs['fs_xlabel'])
		#---! the face color makes things too busy: ax.set_facecolor('black')
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		patch = mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)
		patches.extend([
			mpl.lines.Line2D([],[],color='k',lw=4),
			mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)])
		if False: labels.extend([r'$\mathrm{\langle N_{waters} \rangle}$','zero\nvariance'])
		if snum==3: legend = ax.legend(patches,labels,loc='upper left',fontsize=figspecs['fs_legend'],
			#---no more than 8 items in a column
			bbox_to_anchor=(1.05,0.0,1.,1.),labelspacing=1.2,ncol=int(np.ceil(len(sns_this)/8.)),
			handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
		if snum in [0,2]: axtop.set_title(work.meta[sn]['ptdins_label'],fontsize=16)
		axtop.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
	fn_fig = 'fig.hydration_lenticular.combo'
	picturesave(fn_fig,work.plotdir,backup=False,version=True,
		meta=dict(atom_filter=atom_filter,distance_metric=distance_metric,**spec),extras=[legend])

if 'ultimate' in routine:

	do_single = False
	do_horizontal = True
	if do_single and do_horizontal: raise Exception('incompatible')

	import matplotlib.patheffects as path_effects

	# last edited on 2018.08.06 and then made vertical on 2019.08.22 
	#   by switching the figsize and outer grid. see the if False below
	figsize=(18,8)
	if not do_horizontal: figsize=(8,18)
	if do_single: figsize = (10,10)
	figspecs = {'fs_ticks':16,'fs_ylabel':20,'fs_xlabel':20,'fs_legend':14}
	mean_marker_style = ['bar','dot'][-1]
	if do_horizontal:
		layout = {'out':{'grid':[1,2],'wspace':0.15},
			'ins':[{'grid':[2,1],'hratios':[1,1],'hspace':0.02} for i in range(2)]}
	else:
		layout = {'out':{'grid':[2,1],'hspace':0.2},
			'ins':[{'grid':[2,1],'hratios':[1,1],'hspace':0.02} for i in range(2)]}
	if do_single:
		layout = {'out':{'grid':[1,1],'hspace':0.2},
			'ins':[{'grid':[2,1],'hratios':[1,1],'hspace':0.02} for i in range(2)]}
	axes,fig = panelplot(layout,figsize=figsize)
	extras,legendspec = [],[]
	def color_by_simulation(sn):
		return colorize(work.meta[sn],comparison='asymmetric_all')
	# assign simulations to plots
	if do_single:
		panels = [
			{'ax':axes[1],'axtop':axes[0],'sns':['membrane-v531','membrane-v532']},]
	else:
		panels = [
			{'ax':axes[0][1],'axtop':axes[0][0],'sns':['membrane-v531','membrane-v532']},
			{'ax':axes[1][1],'axtop':axes[1][0],'sns':['membrane-v533','membrane-v534']}]
	sns_all = [i for j in [p['sns'] for p in panels] for i in j]
	for panelspec in panels[:(1 if do_single else None)]:
		counted,means = {},{}
		max_val = int(max([np.concatenate(postdat[sn].values()).max() for sn in panelspec['sns']]))
		ybins = np.arange(max_val+2)-0.5
		for sn in panelspec['sns']:
			means[sn] = {}
			dat = postdat[sn]
			nbins = len(dat)
			raw = np.zeros((nbins,max_val+1))
			for bnum,bin_this in enumerate(sorted(dat.keys())):
				if len(dat[bin_this])>0:
					counts,_ = np.histogram(dat[bin_this],bins=ybins,normed=True)
					raw[bnum] = counts
					means[sn][bin_this] = dat[bin_this].mean()
			raw_colors = (np.concatenate((np.tile(np.ones(raw.shape),(3,1,1)),[raw]))*np.tile(
				np.array(mpl.colors.to_rgba(color_by_simulation(sn))),
				(raw.shape[0],raw.shape[1],1)).transpose((2,0,1))).transpose((2,1,0))
			#! raw_bk = np.array(mpl.colors.to_rgba(color_by_simulation(sn)))*np.array([1.,1.,1.,raw])
			counted[sn] = raw_colors
			legendspec.append(dict(name='%s\n%s%s'%(
				work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'],
				#! removed for vertical '\n(%s)'%work.meta[sn]['composition_name'] if True else ''),
				''),
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sn),alpha=0.5)))
		if len(panelspec['sns'])!=2: raise Exception
		if len(list(set([c.shape for c in counted.values()])))!=1: raise Exception
		folded = np.concatenate([[counted[panelspec['sns'][0]][:,i],counted[panelspec['sns'][1]][:,i]] 
			for i in range(nbins)]).transpose((2,0,1))
		ax = panelspec['ax']
		xticks = np.array(sorted(set([i for j in dat.keys() for i in j]))).astype(float)
		extent = [xticks.min(),xticks.max(),-0.5,max_val+0.5]
		ax.imshow(folded.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['binary'],
			extent=extent)
		# plot the means
		for snum,sn in enumerate(panelspec['sns']):
			for k,v in means[sn].items():
				x = [float(i) for i in k]
				width = x[1]-x[0]
				offset = snum*width/2.
				if mean_marker_style=='bar':
					seg = ax.plot([x[0]+offset,x[0]+offset+width/2],[v,v],
						color=color_by_simulation(sn),lw=5,zorder=4)
					for s in seg: 
						s.set_solid_capstyle('butt')
						s.set_path_effects([path_effects.Stroke(linewidth=9,
							foreground='k',capstyle='butt'),path_effects.Normal()])
				elif mean_marker_style=='dot':
					ax.plot([x[0]+width/4+offset],[v],'o',
						markersize=3,markerfacecolor=color_by_simulation(sn),markeredgewidth=1,
						markeredgecolor='k',lw=0,zorder=4)
				else: raise Exception
		#! use auto for wider view ax.set_aspect((extent[1]-extent[0])/(extent[3]-extent[2]))
		ax.set_aspect('auto')
		ax.grid(which='minor',zorder=0)
		ax.yaxis.grid(which='major',color='k',linestyle='-',linewidth=0,alpha=0,zorder=3)
		ax.set_axisbelow(True)
		ax.set_xlim((extent[0],extent[1]))
		ax.set_ylim((0.5,max_val+0.5))
		ax.set_xticks(xticks,minor=True)
		xticks_major = np.array([1.5,2.0,2.5,3.0,3.5])/10.
		ax.set_xticks(xticks_major,minor=False)
		ax.set_xticklabels(['%.1f'%(float(i)*10) for i in xticks_major],rotation=-90,fontsize=14)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.yaxis.grid(True,which='major',zorder=3)
		ax.xaxis.grid(True,which='major',zorder=3)
		ax.set_yticks(np.arange(max_val+1))
		ax.set_yticks(ybins,minor=True)
		axtop = panelspec['axtop']
		axtop.set_xlim((extent[0],extent[1]))
		counts = np.array([[postdat[sn][k].shape[0]/data[sn]['data']['nframes'] 
			for k in sorted(postdat[sn].keys())] for sn in panelspec['sns']])
		if sorted(postdat[panelspec['sns'][0]].keys())!=sorted(postdat[panelspec['sns'][1]].keys()):
			raise Exception
		# trick where the bin edges are actually paired already
		xbins = np.array(sorted(postdat[sn].keys())).astype(float)
		xvals = np.concatenate(zip(xbins[:,0],(xbins[:,0]+xbins[:,1])/2.))
		colors = np.concatenate(zip(*np.array([[color_by_simulation(sn) 
			for sn in panelspec['sns'] ] for i in xvals]).transpose((1,0,2))))
		axtop.bar(xvals,np.concatenate(zip(*counts)),color=colors,width=xvals[1]-xvals[0],align='edge')
		axtop.set_xticks(xticks)
		axtop.set_xlim((extent[0],extent[1]))
		#! frankly I have no idea how the equal-sized axes were accomplished in the other versions
		ax.get_shared_x_axes().join(ax,axtop)
		axtop.yaxis.grid(True,which='major',zorder=0)
		axtop.xaxis.grid(True,which='major',zorder=0)
		axtop.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		axtop.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		axtop.set_xticklabels([])
		axtop.set_axisbelow(True)
		axtop.set_title(work.meta[sn]['ptdins_label'],fontsize=16)
		axtop.set_ylabel(r'$\mathrm{N_{ions}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_ylabel(r'$\mathrm{N_{waters}}$',fontsize=figspecs['fs_ylabel'])
		ax.set_xlabel('cation-lipid distance ($\mathrm{\AA}$)',fontsize=figspecs['fs_xlabel'])
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	patch = mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)
	patches.extend([
		mpl.lines.Line2D([],[],color='k',lw=4),
		mpl.lines.Line2D([],[],color='k',marker='+',markersize=30,lw=0,mew=4)])
	# changed to the first axis when I made it vertical
	if False:
		legend = ax.legend(patches,labels,loc='upper left',fontsize=figspecs['fs_legend'],
			bbox_to_anchor=(1.05,0.0,1.,1.),labelspacing=1.2,ncol=int(np.ceil(len(sns_all)/8.)),
			handleheight=2.0,markerscale=0.5,shadow=False,fancybox=False)
	else:
		ax = axes[0][0] if not do_single else axes[0]
		legend = ax.legend(patches,labels,loc='upper center',fontsize=figspecs['fs_legend'],
			bbox_to_anchor=(0.5,1.4,),labelspacing=1.2,ncol=len(sns_all),
			handleheight=2.0,markerscale=0.5,shadow=False,fancybox=False)
	frame = legend.get_frame()
	frame.set_edgecolor('w')
	frame.set_facecolor('w')
	extras.append(legend)
	fn_fig = 'fig.hydration_lenticular.imshow%s%s'%(
		('.single' if do_single else ''),('.wide' if do_horizontal else ''))
	picturesave(fn_fig,work.plotdir,backup=False,version=True,form='pdf',
		meta=dict(atom_filter=atom_filter,distance_metric=distance_metric,**spec),extras=extras)
