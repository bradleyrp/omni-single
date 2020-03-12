#!/usr/bin/env python

import copy

def msd_fit_deprecated(time,displacement,leftcut=0.1,rightcut=0.9,dims=3):

	"""
	Fit mean-squared displacement curves to determine the diffusion coefficient.
	"""
	
	if dims==3: factor = 6.
	elif dims==2: factor = 4.
	else: raise Exception('[ERROR] dimensions must be 2 or 3 for diffusion')
	xvals,yvals = time,displacement
	datbuf = [int(len(yvals)*leftcut)+2,int(len(yvals)*rightcut)+1]
	[q,r],covariance = np.polyfit((xvals[datbuf[0]:datbuf[1]]),(yvals[datbuf[0]:datbuf[1]]),1,cov=True)
	return q/factor,r/factor
	
def msd_fit(time,displacement,cuts=None,dims=3,factor=None):

	"""
	Fit mean-squared displacement curves to determine the diffusion coefficient.
	Uses cutofss in the same units as time.
	"""
	
	if factor==None:
		if dims==3: factor = 6.
		elif dims==2: factor = 4.
	xvals,yvals = time,displacement
	#---! added log recently
	ins = np.where(np.all([xvals>=np.log(cuts[0]),xvals<=np.log(cuts[1])],axis=0))[0]
	try: q,r = np.polyfit(xvals[ins],yvals[ins],1)
	except: 
		print "ERROR"
		import pdb;pdb.set_trace()
	return q/factor,r/factor

###---lifted from legacy specs/figures.py

import brewer2mpl
import matplotlib.patheffects as path_effects

def barspecs_stylized(white="w"):

	"""
	"""

	#---listing of hatch patterns	
	patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ][:len(sns)]

	colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
		for key,val in {
			'red':0,
			'blue':1,
			'green':2,
			'purple':3,
			'orange':4,
			'yellow':5,
			'brown':6,
			'pink':7,
			'grey':8,
			}.items()])
	colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")

	#---prepare the dictionaries for the bar plots
	colors_ions = {
		'NA':'green',
		'Na,Cal':'green',
		'MG':'red',
		'Cal':'blue',
		'K':'grey',
		}
	hatches_lipids = {
		'PI2P':'//',
		'P35P':'-',
		'PIPU':'xx',
		'PIPP':'++',
		'SAPI':'',
		}

	barspecs = dict([(sn,{}) for sn in sns])
	for sn in sns:
		barspecs[sn]['hatch'] = hatches_lipids[work.meta[sn]['ptdins_resname']]
		if work.meta[sn]['cation'] == 'Na,Cal':
			barspecs[sn]['hatch'] = '.'
		barspecs[sn]['edgecolor'] = colors[colors_ions[work.meta[sn]['cation']]]
		if work.meta[sn]['ptdins_resname']=='SAPI':
			barspecs[sn]['color'] = tuple(list(barspecs[sn]['edgecolor'])+[0.25])
		else: barspecs[sn]['color'] = mpl.colors.ColorConverter().to_rgb(white)
		for c in ['color','edgecolor'][:]:
			barspecs[sn][c] = tuple(list(barspecs[sn][c])+[0.65])
		#---repeat is necessary
		if work.meta[sn]['ptdins_resname']=='SAPI':
			barspecs[sn]['color'] = tuple(list(barspecs[sn]['edgecolor'])[:3]+[0.25])

	return barspecs

def legend_maker_stylized(ax,title=None,ncol=1,
	bbox=(1.05,0.0,1.,1.),loc='upper left',fs=16,extra_legends=None):

	"""
	"""

	barspecs = barspecs_stylized()
	ion_names = ['K','NA','Na,Cal','MG','Cal']
	marks_sns = [sn for sn in [m[0] for m in [[sn for sn in sns if work.meta[sn]['cation']==i] 
		for i in ion_names] if m]]
	ion_labels = [work.meta[sn]['ion_label'] for sn in marks_sns]
	rectangle_specs = dict([(sn,{}) for sn in marks_sns])
	for ss,sn in enumerate(marks_sns):
		if ion_names[ss] == 'Na,Cal': 
			rectangle_specs[sn]['hatch'] = '..'
			rectangle_specs[sn]['fc'] = 'w'
			rectangle_specs[sn]['edgecolor'] = barspecs[sn]['edgecolor']
			rectangle_specs[sn]['lw'] = 3
		else: 
			rectangle_specs[sn]['fc'] = barspecs[sn]['edgecolor']
			rectangle_specs[sn]['lw'] = 0
	patches = [mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]) for sn in marks_sns]
	for ptdins_resname in ['PI2P','P35P','SAPI'][:-1]:
		sn = [sn for sn in sns if work.meta[sn]['ptdins_resname']==ptdins_resname][0]
		ion_labels += [work.meta[sn]['ptdins_label']]
		patches += [mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=0.5,fc='w',lw=3,
			hatch=barspecs[sn]['hatch']*(1 if ptdins_resname=='PI2P' else 2))]
	if extra_legends:
		for i,j in extra_legends:
			patches.append(i)
			ion_labels.append(j)
	legend = ax.legend(patches,ion_labels,loc=loc,fontsize=fs,
		ncol=ncol,title=title,bbox_to_anchor=bbox,labelspacing=1.2,
		handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	return legend,patches

def errorbars_stylized(ax,mean,std,counter,**kwargs):

	"""
	NOT CURRENTLY SET FOR BAR GROUPS!
	"""

	width = kwargs.get('width',1.0)
	trim_std = kwargs.get('trim_std',0.5-0.15)
	color = kwargs.get('color','blue')
	shadow_color = kwargs.get('shadow_color','k')
	shadow_lw = kwargs.get('shadow_lw',6)
	lw = kwargs.get('lw',4)
	gap = kwargs.get('gap_bar',0.5)
	xoffset = kwargs.get('xoffset',0.0)

	ss = xoffset + gap/2.0 + counter*(width+gap/2.0)
	dc,dc_std = mean,std
	#---somewhat convoluted way to outline the error bars
	std_lines = [([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc-dc_std,dc-dc_std]),
		([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc+dc_std,dc+dc_std]),
		([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
	std_lines_x = [(ss+trim_std)*width,(ss+1-trim_std)*width,(ss+0.5)*width,
		(ss+0.5)*width,(ss+trim_std)*width,(ss+1-trim_std)*width]
	std_lines_y = [dc-dc_std]*3+[dc+dc_std]*3
	sub_index = [[2,3],[0,1],[4,5]]
	for i,j in sub_index:
		ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
			color=color,lw=lw,solid_capstyle='round',alpha=1,zorder=5,
			path_effects=[path_effects.withStroke(linewidth=shadow_lw+lw,foreground=shadow_color)])
	i,j = sub_index[0]
	ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
		color=color,lw=lw,solid_capstyle='round',alpha=1,zorder=5)

def ppi_bar_stylized(**kwargs):

	"""
	Simple bar plot of the physiological areas per lipid in 2D.
	"""

	figsize = (7,10)

	data = kwargs['data']
	sns = kwargs['sns']
	zero_low = kwargs.get('zero_low',False)
	ylabel = kwargs.get('ylabel',None)

	#---aesthetic settings
	fsbase = 14
	trim = 0.15
	shrink = 0.038
	shrink = 0.05
	trim_std = 0.5-0.15
	width = 1.0

	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.4")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.5,boxstyle="round,pad=0.6")

	#---prepare the figure
	if zero_low:
		#---style which indicates that the plot starts at zero
		axes,fig = panelplot(figsize=figsize,layout={'out':{'grid':[2,1],'hratios':[10,1],'hspace':0.05},
			'ins':[{'grid':[1,1]} for j in range(2)]})
		axdup = [axes[0][0],axes[1][0]]
	else:
		axes,fig = panelplot(figsize=figsize,layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]})
		axdup = [axes[0]]

	barspecs = barspecs_stylized()

	if False:
		#---listing of hatch patterns	
		patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ][:len(sns)]

		###---BEGIN AESTHETIC CHOICES

		import brewer2mpl
		colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
			for key,val in {
				'red':0,
				'blue':1,
				'green':2,
				'purple':3,
				'orange':4,
				'yellow':5,
				'brown':6,
				'pink':7,
				'grey':8,
				}.items()])
		colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")

		#---prepare the dictionaries for the bar plots
		colors_ions = {
			'NA':'green',
			'Na,Cal':'green',
			'MG':'red',
			'Cal':'blue',
			'K':'grey',
			}
		hatches_lipids = {
			'PI2P':'//',
			'P35P':'-',
			'PIPU':'xx',
			'PIPP':'++',
			'SAPI':'',
			}

		barspecs = dict([(sn,{}) for sn in sns])
		for sn in sns:
			barspecs[sn]['hatch'] = hatches_lipids[work.meta[sn]['ptdins_resname']]
			if work.meta[sn]['cation'] == 'Na,Cal':
				barspecs[sn]['hatch'] = '.'
			barspecs[sn]['edgecolor'] = colors[colors_ions[work.meta[sn]['cation']]]
			if work.meta[sn]['ptdins_resname']=='SAPI':
				barspecs[sn]['color'] = tuple(list(barspecs[sn]['edgecolor'])+[0.25])
			else: barspecs[sn]['color'] = mpl.colors.ColorConverter().to_rgb("w")
			for c in ['color','edgecolor'][:]:
				barspecs[sn][c] = tuple(list(barspecs[sn][c])+[0.65])
			#---repeat is necessary
			if work.meta[sn]['ptdins_resname']=='SAPI':
				barspecs[sn]['color'] = tuple(list(barspecs[sn]['edgecolor'])[:3]+[0.25])

		###---END AESTHETIC CHOICES

	#---looping for broken axis
	for which_axis,ax in enumerate(axdup):

		yvals = data
		xvals = np.arange(len(sns))+trim
		for yy,y in enumerate([yvals[sn]['mean'] for sn in sns]):
			x = xvals[yy]
			ax.bar(xvals[yy],y,
				width=1.0-trim*2,
				lw=3.0,
				label=work.meta[sns[yy]]['cation'],
				**barspecs[sns[yy]])

		#---moved this to errorbars_stylized for other functions but haven't reached consensus
		if True:
			#---also plot standard deviations
			for ss,sn in enumerate(sns):
				dc,dc_std = yvals[sn]['mean'],yvals[sn]['std']
				#---somewhat convoluted way to outline the error bars
				std_lines = [([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc-dc_std,dc-dc_std]),
					([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc+dc_std,dc+dc_std]),
					([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
				std_lines_x = [(ss+trim_std)*width,(ss+1-trim_std)*width,(ss+0.5)*width,
					(ss+0.5)*width,(ss+trim_std)*width,(ss+1-trim_std)*width]
				std_lines_y = [dc-dc_std]*3+[dc+dc_std]*3
				sub_index = [[2,3],[0,1],[4,5]]
				for i,j in sub_index:
					ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
						color=barspecs[sn]['edgecolor'],
						lw=4,solid_capstyle='round',alpha=1,zorder=5,
						path_effects=[path_effects.withStroke(linewidth=10,foreground='k')])
				i,j = sub_index[0]
				ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
					color=barspecs[sn]['edgecolor'],
					lw=4,solid_capstyle='round',alpha=1,zorder=5)
		else:
			#---reworking to match error bars with the regular bars
			width = width - trim
			gap = trim
			for ss,sn in enumerate(sns):
				errorbars_stylized(ax,yvals[sn]['mean'],yvals[sn]['std'],ss,gap=gap,
					width=width,trim_std=trim_std,color=barspecs[sn]['edgecolor'])

		mean_vals = [data[sn]['mean'] for sn in sns]
		ymin,ymax = np.mean(mean_vals)*0.9,np.mean(mean_vals)*1.05

		if zero_low and which_axis == 0:
			ax.set_ylim(ymin,ymax)
			ax.set_yticks(np.arange(np.round(ymin,2),np.round(ymax,2)+10**-2,10**-2))
		if not zero_low or which_axis == 1:
			ax.set_xticklabels([work.meta[sn]['ptdins_label'] for sn in sns],rotation=90,fontsize=fsbase+6)
			ax.set_xticks(xvals+width/2.0)
		elif zero_low and which_axis == 0: ax.set_xticks([])
		ax.set_xlim(0.0,len(xvals))

		loc = 'upper right'
		title = None

		#---important aesthetic choices below
		if not zero_low or which_axis==0:
			if False:
				ion_names = ['K','NA','Na,Cal','MG','Cal']
				marks_sns = [sn for sn in [m[0] for m in [[sn for sn in sns if work.meta[sn]['cation']==i] 
					for i in ion_names] if m]]
				ion_labels = [work.meta[sn]['ion_label'] for sn in marks_sns]
				rectangle_specs = dict([(sn,{}) for sn in marks_sns])
				for ss,sn in enumerate(marks_sns):
					if ion_names[ss] == 'Na,Cal': 
						rectangle_specs[sn]['hatch'] = '..'
						rectangle_specs[sn]['fc'] = 'w'
						rectangle_specs[sn]['edgecolor'] = barspecs[sn]['edgecolor']
						rectangle_specs[sn]['lw'] = 3
					else: 
						rectangle_specs[sn]['fc'] = barspecs[sn]['edgecolor']
						rectangle_specs[sn]['lw'] = 0
				patches = [mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]) for sn in marks_sns]
				for ptdins_resname in ['PI2P','P35P','SAPI'][:-1]:
					sn = [sn for sn in sns if work.meta[sn]['ptdins_resname']==ptdins_resname][0]
					ion_labels += [work.meta[sn]['ptdins_label']]
					patches += [mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=0.5,fc='w',lw=3,
						hatch=barspecs[sn]['hatch']*(1 if ptdins_resname=='PI2P' else 2))]
				legend = ax.legend(patches,ion_labels,loc='upper left',fontsize=art['fs']['legend']+6,
					ncol=1,title=title,bbox_to_anchor=(1.05,0.0,1.,1.),labelspacing=1.2,
					handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
			legend,patches = legend_maker_stylized(ax)
			#---using shadow, fancybox above instead of thicker lines
			if False:
				legend.get_frame().set_linewidth(2.0)
				for axis in ['top','bottom','left','right']:
					ax.spines[axis].set_linewidth(2.0)

		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ylabelform = kwargs.get('ylabel_format',None)
		if ylabelform:
			ax.set_yticklabels([ylabelform(i) for i in ax.get_yticks()],fontsize=fsbase+10)
		else: ax.set_yticklabels(['%.2f'%i for i in ax.get_yticks()],fontsize=fsbase+10)
		if not zero_low or which_axis == 0:
			if ylabel: ax.set_ylabel(ylabel,fontsize=fsbase+8,labelpad=20)

	if zero_low:
		ax = axdup[1]
		interval = [i[1]-i[0] for i in [axdup[0].get_yticks()]][0]
		ax.set_ylim(0,interval)
		ax.set_yticks([0])
		axdup[0].spines['bottom'].set_linewidth(0)
		axdup[1].spines['top'].set_linewidth(0)
	
		#---from http://stackoverflow.com/questions/5656798/\
		#---...python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis	
		#---size of break marks in axis coordinates
		break_pos = lambda rx,ry,d=0.008 : {
			'top_l':[(1-d*rx,1+d*rx),(-d*ry,+d*ry)],
			'bot_l':[(1-d*rx,1+d*rx),(1-d*ry,1+d*ry)],
			'top_r':[(-d*rx,d*rx),(-d*ry,+d*ry)],
			'bot_r':[(-d*rx,d*rx),(1-d*ry,1+d*ry)],
			}
		bbox = axdup[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		#---set set the aspect here for a break in the y-direction
		aspect = 1.0,(bbox.width/bbox.height)
		axspec = dict(transform=axdup[0].transAxes,color='k',clip_on=False)
		axdup[0].plot(*break_pos(*aspect)['top_l'],lw=2,**axspec)
		axdup[0].plot(*break_pos(*aspect)['top_r'],lw=2,**axspec)
		bbox = axdup[1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
		aspect = 1.0,(bbox.width/bbox.height)
		axspec.update(transform=axdup[1].transAxes)
		axdup[1].plot(*break_pos(*aspect)['bot_l'],lw=2,**axspec)
		axdup[1].plot(*break_pos(*aspect)['bot_r'],lw=2,**axspec)
		bbox = axdup[0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

	picturesave(kwargs['name'],work.plotdir,backup=False,
		version=True,meta={},pdf=True,extras=[legend])
	plt.close()

	"""
	http://stackoverflow.com/questions/5656798/python-matplotlib-is-there-a-way-to-make-a-discontinuous-axis
	"""

def ppi_bar_stylized_panel(ax,**kwargs):

	"""
	REPETITIVE WITH THE PLOT ABOVE.
	The non-panel version does a standalone.
	This one was designed to work for single axes so you could use the style on many plots.
	It would be useful to unify them all.
	Had to remove zero-low, which requires multiple axes and hence doesn't make sense here.
	"""

	outvals = {}
	data = kwargs['data']
	sns = kwargs['sns']
	ylabel = kwargs.get('ylabel',None)
	do_legend = kwargs.get('legend',True)
	altback = kwargs.get('altback',None)
	std_dev_ec = kwargs.get('std_dev_ec','k')
	errorbars_only = kwargs.get('errorbars_only',False)

	#---aesthetic settings
	fsbase = 14
	trim = 0.15
	shrink = 0.038
	shrink = 0.05
	trim_std = 0.5-0.15
	width = 1.0

	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.4")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.5,boxstyle="round,pad=0.6")

	#---alternate background fill
	if altback: barspecs = barspecs_stylized(white=altback)
	else: barspecs = barspecs_stylized()

	yvals = data
	xvals = np.arange(len(sns))#+trim
	for yy,y in enumerate([yvals[sn]['mean'] for sn in sns]):
		if not errorbars_only: ax.bar((yy+0.5)*width,y,
			width=1.0-trim*2,
			lw=3.0,
			label=work.meta[sns[yy]]['cation'],
			**barspecs[sns[yy]])

	#---! NOTE THAT INCOMING ERROR BAR SETTINGS SHOULD BE ROUTED TO errorbars_stylized
	eb_alpha = kwargs.get('errorbar_alpha',1.0)

	#---moved this to errorbars_stylized for other functions but haven't reached consensus
	if True:
		#---also plot standard deviations
		for ss,sn in enumerate(sns):
			dc,dc_std = yvals[sn]['mean'],yvals[sn]['std']
			#---somewhat convoluted way to outline the error bars
			std_lines = [([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc-dc_std,dc-dc_std]),
				([(ss+trim_std)*width,(ss+1-trim_std)*width],[dc+dc_std,dc+dc_std]),
				([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
			std_lines_x = [(ss+trim_std)*width,(ss+1-trim_std)*width,(ss+0.5)*width,
				(ss+0.5)*width,(ss+trim_std)*width,(ss+1-trim_std)*width]
			std_lines_y = [dc-dc_std]*3+[dc+dc_std]*3
			sub_index = [[2,3],[0,1],[4,5]]
			for i,j in sub_index:
				ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
					color=barspecs[sn]['edgecolor'],
					lw=4,solid_capstyle='round',alpha=eb_alpha,zorder=5,
					path_effects=[path_effects.withStroke(linewidth=10,foreground=std_dev_ec)])
			i,j = sub_index[0]
			ax.plot([std_lines_x[i],std_lines_x[j]],[std_lines_y[i],std_lines_y[j]],
				color=barspecs[sn]['edgecolor'],
				lw=4,solid_capstyle='round',alpha=eb_alpha,zorder=5)
	else:
		#---reworking to match error bars with the regular bars
		width = width - trim
		gap = trim
		for ss,sn in enumerate(sns):
			errorbars_stylized(ax,yvals[sn]['mean'],yvals[sn]['std'],ss,gap=gap,
				width=width,trim_std=trim_std,color=barspecs[sn]['edgecolor'])

	mean_vals = [data[sn]['mean'] for sn in sns]
	ymin,ymax = np.mean(mean_vals)*0.9,np.mean(mean_vals)*1.05

	#---note that the labels are perfectly flish with the beginning of the bar, for whatever reason
	ax.set_xticklabels([work.meta[sn]['ptdins_label'] for sn in sns],rotation=90,fontsize=fsbase+6)
	ax.set_xticks(xvals+width/2.0)
	ax.set_xticks([(yy+0.5)*width for yy in range(len(xvals))])
	ax.set_xlim(0,len(xvals))

	loc = 'upper right'
	title = None

	#---important aesthetic choices below
	if do_legend: 
		legend,patches = legend_maker_stylized(ax,extra_legends=kwargs.get('extra_legends',None),
			ncol=kwargs.get('ncol',1))
		frame = legend.get_frame()
		frame.set_edgecolor('black')
		frame.set_facecolor('white')
		frame.set_linewidth(2)
		outvals['legend'] = legend
		outvals['patches'] = patches

	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ylabelform = kwargs.get('ylabel_format',None)
	if ylabelform:
		ax.set_yticklabels([ylabelform(i) for i in ax.get_yticks()],fontsize=fsbase+10)
	else: ax.set_yticklabels(['%.2f'%i for i in ax.get_yticks()],fontsize=fsbase+10)
	if ylabel: ax.set_ylabel(ylabel,fontsize=fsbase+8,labelpad=20)

	return outvals

def barmaker(ax,yvals,**kwargs):

	"""
	standardized bar plotter with groupings
	bars are plotted as a list of lists
	"""

	#---defaults have a space between bar groups which are otherwise flush
	width = kwargs.pop('width',1.0)
	gap_out = kwargs.pop('gap',0.1)
	gap_in = kwargs.pop('gap_in',0.0)
	zero_threshold = kwargs.pop('zero_threshold',0.0)
	lw_zero = kwargs.pop('lw_zero',3.0)
	barspecs_all = kwargs.get('barspecs_all',{})
	barspecs = kwargs.get('barspecs',{})
	barspecs_stack = kwargs.get('barspecs_stack',{})
	xoffset = kwargs.get('xoffset',0)

	#---if there are no bar groups and/or stacks we reformat the objects for the triple loop
	if all([type(i)!=list for i in yvals]): yvals = [[[i]] for i in yvals]
	elif all([all([type(j)!=list for j in i]) for i in yvals]):
		yvals = [[[j] for j in i] for i in yvals]
	assert all([all([type(j)==list for j in i]) for i in yvals]),'yvals have weird dimensions'

	ymax,xpos,xvals = 0,xoffset,[]
	#---loop over groups of bars
	for ii,i in enumerate(yvals):
		xpos += gap_out
		xvals_grp = []
		#---loop over bars within the group
		for jj,j in enumerate(i):
			xvals_grp.append(xpos)
			bottom = 0
			#---loop over stacks of bars
			for kk,k in enumerate(j):
				#---assemble the right specs for the bar (by group) and stack position
				this_barspecs = dict(barspecs_all)
				if barspecs: this_barspecs.update(**barspecs[ii][jj])
				if barspecs_stack: this_barspecs.update(**barspecs_stack[ii][jj][kk])
				#---zero thresholds get a nice line to remind you they are present
				if k <= zero_threshold and bottom == 0:
					ax.plot([xpos,xpos+width],[k,k],color=this_barspecs.get('color','k'),
						lw=2,zorder=2,clip_on=False)
					ax.plot([xpos,xpos+width],[k,k],color='k',lw=4,zorder=1,clip_on=False)
				ax.bar(xpos,k,width=width,bottom=bottom,**this_barspecs)
				ymax = ymax if ymax>k+bottom else k+bottom
				bottom += k
			xpos += width
			xpos += gap_in
		xvals.append(xvals_grp)

	return {'ymax':ymax,'xvals':xvals,'width':width}

###---MAIN

#---settings
withchar = ','
cutoff_time = [[1,20],[8,72]][0]
fit_type = 'linear'
routine = [
	'average_diffusion',
	'average_diffusion_simple_custom',
	'all_lipids_msds_fitline',
	'lipid_diffusion_histograms_by_lipid',
	'lipid_diffusion_histograms',
	'diffusion_scaling',
	'average_diffusion_simple_custom2',
	][:]

if 'data' not in globals(): 
	data,calc = plotload(plotname,work)
	next_script = sys.argv[0]
if 'post' not in globals():
	post = dict([(sn,{}) for sn in work.sns()])
	post_log = dict([(sn,{}) for sn in work.sns()])
	for sn in [s for s in work.sns() if s in data.keys()]:
		dts = data[sn]['data']['dt']/1000.
		msd = data[sn]['data']['msd']
		ndts,nmol = np.shape(msd)
		names = np.unique(data[sn]['data']['resnames'])
		post[sn] = dict([(n,[]) for n in names])
		post_log[sn] = dict([(n,[]) for n in names])
		for ind in range(nmol):
			#---note that the msd only holds displacements so we square it before fitting
			if 0: q,r = msd_fit_deprecated(dts,msd[:,ind]**2,leftcut=0.1,rightcut=0.5,dims=2)
			else: q,r = msd_fit(dts,msd[:,ind]**2,cuts=cutoff_time,dims=2)
			qlog,rlog = msd_fit(np.log(dts),np.log(msd[:,ind]**2),cuts=cutoff_time,dims=2,factor=1)
			post[sn][data[sn]['data']['resnames'][ind]].append(q)
			post_log[sn][data[sn]['data']['resnames'][ind]].append(qlog)
	sns_sub = [s for s in work.sns() if s in data.keys()]

if 'average_diffusion' in routine:

	"""
	Plot (dots of) the average diffusion rate for each lipid type in each simulation.
	"""

	allres = work.vars['selectors']['resnames_lipid_chol']
	resnames_pip2 = work.vars['selectors']['resnames_PIP2']
	axes,fig = panelplot(figsize=(8,8),layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]})
	ax = axes[0]
	counter,centers = 0,[]
	for ss,sn in enumerate(sns_sub):
		color = colorize(work.meta[sn],comparison='protonation')
		dat = data[sn]['data']
		reslist = [r for r in allres if r in np.unique(dat['resnames'])]
		#---scale by 10**6 to convert from nm2/ps to um2/s
		if 0: ax.scatter([ss for r in reslist],[np.mean(post[sn][r])*10**6/1000. for r in reslist],lw=0,s=100,
			color=[colorize(work.meta,resname=r) for r in reslist])
		dc = np.mean(np.concatenate([post[sn][r] for r in reslist])*10**6/1000.)
		ax.bar(range(counter,counter+len(reslist)),[np.mean(post[sn][r])*10**6/1000.-dc for r in reslist],
			width=1.0,lw=0,color=[colorize(work.meta,resname=r) for r in reslist],
			alpha=1,zorder=2,bottom=dc)
		ax.plot(range(counter,counter+len(reslist)+1),np.ones(len(reslist)+1)*dc,
			'k',lw=3,solid_capstyle='round',alpha=1,zorder=3)
		centers.append(np.mean(np.arange(counter,counter+len(reslist))+1))
		counter += len(reslist)+1
	ax.set_xticks(centers)
	ax.set_xticklabels([work.meta[sn]['ion_label']+withchar+work.meta[sn]['ptdins_label'] 
		for sn in sns_sub],rotation=90)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ax.set_xlim(-1,counter)
	ax.set_title('lipid diffusion rates')
	ax.set_ylabel('$\mathrm{D\,{\mu m}^{2} {s}^{-1}}$')
	patches = [(mpl.patches.Patch(color=colorize(work.meta,resname=(r if r!='PtdIns' else 'PI2P')),label=r),r) 
		for r in [a for a in allres if a not in resnames_pip2]+['PtdIns']]
	legend = ax.legend(*zip(*patches),bbox_to_anchor=(1,0.5),loc='center left')
	picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={},extras=[legend])
	plt.close()

if 'average_diffusion_simple_custom' in routine:

	"""
	Custom plot of physiological bilayer diffusion rates, simplified.
	"""

	#---! problem with meta causes many figures to be created. possibly an issue with omnicalc 
	#---! actually everything looks cool but check it
	fsbase = 18
	ymin,ymax = 0,0
	shrink = 0.045 #---correct lines overhanging bars
	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.4")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.4")
	allres = work.vars['selectors']['resnames_lipid_chol']
	resnames_pip2 = work.vars['selectors']['resnames_PIP2']
	axes,fig = panelplot(figsize=(7,10),layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]})
	ax = axes[0]
	counter,centers,scs = 0,[],[]
	#---re-order by cation
	cations_order = ['NA','K','MG','Cal']
	sns = [sn for cat in cations_order for sn in work.specs['collections']['all']
		if work.meta[sn]['composition_name'] == 'asymmetric' and work.meta[sn]['cation']==cat]
	for ss,sn in enumerate(sns):
		if fit_type == 'linear': post_this = post[sn]
		elif fit_type == 'log': post_this = post_log[sn]
		color = colorize(work.meta[sn],comparison='protonation')
		dat = data[sn]['data']
		reslist = [r for r in allres if r in np.unique(dat['resnames'])]
		width = len(reslist)
		#---scale by 10**6 to convert from nm2/ps to um2/s
		#---hack to show the PS dot for sn=='membrane-v538
		offsets = [{('membrane-v538','DOPS'):0.75,('membrane-v538','PI2P'):-0.75}.\
			get((sn,r),0) for r in reslist]
		offsets = [{}.get((sn,r),0) for r in reslist]
		scs.append(ax.scatter([width*ss+width/2.0+offsets[rr] for rr,r in enumerate(reslist)],
			[np.mean(post_this[r])*10**6/1000. for r in reslist],lw=1.5,s=175,edgecolor='w',
			color=[colorize(work.meta,resname=r) for r in reslist],zorder=4))
		#---tiny standard deviation plots (repetitive)
		if False:
			for deviate in [-1,1]:
				scs.append(ax.scatter([width*ss+width/2.0+offsets[rr] for rr,r in enumerate(reslist)],
				[np.mean(post_this[r])*10**6/1000.+deviate*np.std(post_this[r])*10**6/1000. 
				for r in reslist],lw=0,s=100,edgecolor='w',
				color=[colorize(work.meta,resname=r) for r in reslist],zorder=4))
		dc = np.mean(np.concatenate([post_this[r] for r in reslist])*10**6/1000.)
		dc_std = np.std(np.concatenate([post_this[r] for r in reslist])*10**6/1000.)
		ymax = dc+dc_std if dc+dc_std > ymax else ymax
		if 0: ax.bar(range(counter,counter+width),[np.mean(post_this[r])*10**6/1000.-dc for r in reslist],
			width=1.0,lw=0,color=[colorize(work.meta,resname=r) for r in reslist],
			alpha=1,zorder=2,bottom=dc)
		ax.plot([counter+shrink*width,counter+width-shrink*width],np.ones(2)*dc,
			color=color,lw=5,solid_capstyle='round',alpha=1,zorder=3)
		ax.bar(ss*width,dc,width=width,alpha=0.35,zorder=2,lw=0,color=color,edgecolor=color)
		#---also plot standard deviations
		std_lines = [([(ss+0.25)*width,(ss+0.75)*width],[dc-dc_std,dc-dc_std]),
			([(ss+0.25)*width,(ss+0.75)*width],[dc+dc_std,dc+dc_std]),
			([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
		for l in std_lines: 
			ax.plot(l[0],l[1],color=color,lw=4,solid_capstyle='round',alpha=1,zorder=3)
			ymax = l[1] if dc+dc_std > ymax else ymax
		centers.append(np.mean(np.arange(counter,counter+width)+1))
		counter += width+0
	ax.set_xticks(centers)
	if False:
		ax.set_xticklabels([work.meta[sn]['ion_label']+withchar+work.meta[sn]['ptdins_label'] 
			for sn in sns],rotation=45,ha='right',fontsize=fsbase)
	outer_labels = []
	ax.set_xticklabels([])
	for snum,sn in enumerate(sns):
		tb = ax.text(width*(snum+0.5),(ymax-ymin)*-0.02,work.meta[sn]['ion_label'],
			bbox=tagbox_ion,ha='center',va='top',rotation=0,
			color='k',fontsize=fsbase-2)
		outer_labels.append(tb)
		text = work.meta[sn]['ptdins_label']
		tb = ax.text(width*(snum+0.5),(ymax-ymin)*0.02,text,
			bbox=tagbox_ptdins,rotation=90,ha="center",va="bottom",color='k',fontsize=fsbase)
	ax.set_ylim(ymin,ymax*1.1)
	ax.set_yticklabels(['%.f'%i for i in ax.get_yticks()],fontsize=fsbase+2)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ax.set_xlim(-1,counter+1)
	ax.set_title('lipid diffusion rates',fontsize=fsbase+4)
	ax.set_ylabel('$\mathrm{D\,{\mu m}^{2} {s}^{-1}}$',fontsize=fsbase+2)
	patches = [(mpl.lines.Line2D([],[],marker='o',markersize=15,markeredgecolor='w',mew=2,lw=0,
		markerfacecolor=colorize(work.meta,resname=(r if r!='PtdIns' else 'PI2P'),),label=r),r) 
		for r in [a for a in allres if a not in resnames_pip2]+['PtdIns']]
	legend = ax.legend(*zip(*patches),loc='upper right',numpoints=1)
	picturesave('fig.%s.simple'%plotname,work.plotdir,backup=False,
		version=True,meta={'cutoff_time':cutoff_time,'fit_type':fit_type},extras=[legend]+outer_labels)
	plt.close()

if 'all_lipids_msds_fitline' in routine:

	"""
	Demonstrate the lipid diffusion fits in the linear space.
	These fits assume the behavior is diffusive.
	"""
	
	#dt = work.slice(sn)[calc['slice_name']]['all']['skip']
	dt = work.slices[sn]['slices']['current']['skip']
	axes,fig = panelplot(layout=figlayout['summary1'],figsize=(12,12))
	for pnum,sn in enumerate(sns_sub):
		axrow,axcol = figplacer(sn,figplace['summary1'])
		ax = axes[axrow][axcol]
		ax.set_ylabel('MSD $\mathrm{({nm}^{2})}$')
		ax.set_xlabel('time (ns)')
		ax.set_title(' '.join([
			work.meta[sn]['composition_name'],
			work.meta[sn]['ptdins_label'],
			work.meta[sn]['ion_label'],
			]))
		dts = data[sn]['data']['dt']/1000.
		msd = data[sn]['data']['msd']**2
		for m in msd.T[::]:
			ax.plot(dts,m,alpha=0.01,zorder=3)
			q,r = msd_fit(dts,m,dims=2,cuts=cutoff_time)
			ins = np.where(np.all([dts>=cutoff_time[0],dts<=cutoff_time[1]],axis=0))[0]
			#---msd_fit comes with the diffusion factor so we rescale here for the plot
			ax.plot(dts[ins],4*q*dts[ins]+4*r,c='r',lw=1,zorder=2,alpha=0.1)
	for rr,cc in [(ii,jj) for ii,i in enumerate(figplace['summary1']) 
		for jj,j in enumerate(i) if j==None]:
		fig.delaxes(axes[rr][cc])
	picturesave('fig.%s.fits_linear'%plotname,work.plotdir,backup=False,version=True,meta={})
	plt.close()

if 'lipid_diffusion_histograms_by_lipid' in routine:

	"""
	Histograms of ion diffusion coefficients for each lipid type.
	"""

	#---convert nm2/ns to um2/s
	factor = 1000.
	edges = np.arange(-2,20+1,1)
	fl = copy.deepcopy(figlayout['summary1'])
	fl['out']['hspace'] = 0.6
	axes,fig = panelplot(layout=fl,figsize=(12,10))
	for pnum,sn in enumerate(sns_sub):
		axrow,axcol = figplacer(sn,figplace['summary1'])
		ax = axes[axrow][axcol]
		for lipid_type in post[sn].keys():
			counts,edges = np.histogram(factor*np.array(post[sn][lipid_type]),bins=edges,normed=True)
			mids = (edges[:-1]+edges[1:])/2.
			color = colorize(work.meta[sn],resname=lipid_type)
			label = lipid_type
			ax.plot(mids,counts,color=color,lw=2,label=label)
		ax.set_title(work.meta[sn]['composition_name']+' with '+work.meta[sn]['ptdins_label'])
		ax.set_xlabel('D $\mathrm{({\mu m}^{2}{s}^{-1})}$')
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	for rr,cc in [(ii,jj) for ii,i in enumerate(figplace['summary1']) 
		for jj,j in enumerate(i) if j==None]: fig.delaxes(axes[rr][cc])
	for pnum,sn in enumerate(sns_sub):
		axrow,axcol = figplacer(sn,figplace['summary1'])
		ax = axes[axrow][axcol]
		legend = ax.legend(loc='upper right',fontsize=8)
	picturesave('fig.%s.distribution_by_lipid'%plotname,work.plotdir,backup=False,version=True,meta={})
	plt.close()

if 'lipid_diffusion_histograms' in routine:

	"""
	Histograms of ion diffusion coefficients.
	"""

	#---convert nm2/ns to um2/s
	factor = 1000.
	edges = np.arange(-2,20+1,1)
	#axes,fig = panelplot(layout={'out':{'grid':[1,1]},'ins':[{'grid':[2,2],'hspace':0.5}]},figsize=(10,6))
	axes,fig = panelplot(layout=figlayout['summary1'],figsize=(12,12))
	for pnum,sn in enumerate(sns_sub):
		axrow,axcol = figplacer(sn,figplace['summary1'])
		ax = axes[axrow][axcol]
		counts,edges = np.histogram(factor*np.concatenate(post[sn].values()),bins=edges)
		mids = (edges[:-1]+edges[1:])/2.
		color = colorize(work.meta[sn],comparison='protonation')
		label = lipid_type
		ax.set_title(work.meta[sn]['composition_name']+' with '+work.meta[sn]['ptdins_label'])
		ax.plot(mids,counts,color=color,lw=2,label=label)
		ax.fill_between(mids,0,counts,color=color,alpha=0.35,zorder=2,label=label)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.set_xlabel('D $\mathrm{({\mu m}^{2}{s}^{-1})}$')
	picturesave('fig.%s.distribution'%plotname,work.plotdir,backup=False,version=True,meta={})
	plt.close()

if 'diffusion_scaling' in routine:

	"""
	Check diffusion type.
	"""
	
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	edges = np.arange(-1.0,2+0.25,0.25)
	#dt = work.slice(sn)[calc['slice_name']]['all']['skip']
	dt = work.slices[sn]['slices']['current']['skip']
	axes,fig = panelplot(layout=figlayout['summary1'],figsize=(12,12))
	for pnum,sn in enumerate(sns_sub):
		axrow,axcol = figplacer(sn,figplace['summary1'])
		ax = axes[axrow][axcol]
		ax.set_xlabel('scaling coefficient')
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.set_title(' '.join([
			work.meta[sn]['composition_name'],
			work.meta[sn]['ptdins_label'],
			work.meta[sn]['ion_label'],
			]))
		axins = inset_axes(ax,width="30%",height="30%",loc=1,
			bbox_to_anchor=(0.,0.,1.,1.),
            bbox_transform=ax.transAxes,
            borderpad=0)
		dts = data[sn]['data']['dt']/1000.
		msd = data[sn]['data']['msd']**2
		for m in msd.T[::10]:
			axins.plot(np.log(dts),np.log(m),alpha=0.1,zorder=3)
			ins = np.where(np.all([dts>=cutoff_time[0],dts<=cutoff_time[1]],axis=0))[0]
		axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
		axins.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='off')
		ax.axvline(0.,c='k')
		ax.axvline(1.,c='k')
		for lipid_type in post_log[sn].keys():
			counts,edges = np.histogram(post_log[sn][lipid_type],normed=True,bins=edges)
			mids = (edges[:-1]+edges[1:])/2.
			color = colorize(work.meta[sn],resname=lipid_type)
			label = lipid_type
			ax.set_title(work.meta[sn]['composition_name']+' with '+work.meta[sn]['ptdins_label'])
			ax.plot(mids,counts,color=color,lw=2,label=label)
	for rr,cc in [(ii,jj) for ii,i in enumerate(figplace['summary1']) 
		for jj,j in enumerate(i) if j==None]:
		fig.delaxes(axes[rr][cc])
	plt.suptitle('diffusion scaling')
	picturesave('fig.%s.scaling'%plotname,work.plotdir,backup=False,version=True,meta={})
	plt.close()
	
if 'average_diffusion_simple_custom2' in routine:

	#---use canonical orderings
	sns = work.vars['orders']['canon']['asymmetric']
	allres = work.vars['selectors']['resnames_lipid_chol']
	diffusion_consts = dict([(sn,{}) for sn in sns])
	for ss,sn in enumerate(sns):
		dat = data[sn]['data']
		reslist = [r for r in allres if r in np.unique(dat['resnames'])]
		if fit_type == 'linear': post_this = post[sn]
		elif fit_type == 'log': post_this = post_log[sn]
		catted = np.concatenate([post_this[r] for r in reslist])
		diffusion_consts[sn]['mean'] = np.mean(catted)*10**6/1000.
		diffusion_consts[sn]['std'] = np.std(catted)*10**6/1000.

	ppi_bar_stylized(
		name='fig.diffusion_lipids.simple2',
		data=diffusion_consts,
		sns=sns,
		zero_low=False,
		ylabel_format=lambda x:'%d'%x,
		ylabel='$\mathrm{D\,{\mu m}^{2} {s}^{-1}}$',)	
