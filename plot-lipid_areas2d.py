#!/usr/bin/python -i

#! updated from legacy factory!!! everything works except the plot that shows the delta between areas
def unique(k): return list(np.unique(k))

# from calcs/specs/figures.py

if True:

	import brewer2mpl
	import matplotlib.patheffects as path_effects

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
					lw=3.0,align='edge',
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

# end stuff from calcs/specs/figures.py

if False:
	#---plot prep
	if 'plotload' not in globals(): execfile('/etc/pythonstart')
	execfile('./omni/base/header.py')
	from plotter import *
	from base.store import plotload
	execfile('./calcs/specs/figures.py')
	execfile('./calcs/specs/colors.py')
	import numpy as np
	import matplotlib.patheffects as path_effects

#---settings
plotname = 'lipid_areas2d'
withchar = ','
routine = ['detailed_view','simple_view','simple_view_newbars'][:]
print_version = 'v20160923'

#---load everything
if 'data' not in globals(): 
	data,calcs = plotload(plotname)
	#! pulled from legacy and modified here
	keys = dict([(delve(v,
		'calcs','specs','upstream','lipid_mesh','upstream','lipid_abstractor','selector'),k)
		for k,v in calcs.items() if k!='extras'])
	data = [data[keys[k]] for k in ['lipid_com','lipid_chol_com']]
	ns = next_script = sys.argv[0]

if 'detailed_view' in routine:

	axes,fig = panelplot(figsize=(12,10),
		layout={'out':{'grid':[1,1],'hspace':1.0},
		'ins':[{'grid':[1,3],'wspace':0.5}]})
	allres = work.vars['selectors']['resnames_lipid']
	sns = [i for i in work.sns() if i in data[0]]
	counters = [1,1,1]
	for axnum,composition,has_chol in [
		(0,'symmetric',False),
		(1,'asymmetric',False),
		(2,'asymmetric',True)]:
		sns_sub = [sn for sn in sns if work.meta[sn]['composition_name']==composition]
		#---better ordering
		sns_sub = [i for i in work.vars['orders']['areas'][composition] if i in sns_sub]
		if has_chol: allres = work.vars['selectors']['resnames_lipid_chol']
		else: allres = work.vars['selectors']['resnames_lipid']
		tick_positions = []
		for sn in sns_sub:
			ax = axes[axnum]
			dat = data[1 if has_chol else 0][sn]['data']
			resnames = dat['resnames']
			reslist = [i for i in allres if i in unique(resnames)]
			ntypes = len(reslist)
			mono_areas = [dat['areas%d'%mn] for mn in range(2)]
			imono = dat['monolayer_indices']
			monolayer_residues = [[np.concatenate([np.where(np.where(imono==mn)[0]==i)[0] 
				for i in np.where(resnames==r)[0]]) for r in reslist] for mn in range(2)]
			for mn in range(2):
				#---plot base
				ax.plot(range(counters[axnum],counters[axnum]+ntypes+1),np.ones(ntypes+1)*np.mean(mono_areas[mn]),
					'k',lw=3,solid_capstyle='round',alpha=1,zorder=3)
				#---plot deviations
				for rr,r in enumerate(reslist):
					actual_r = reslist.index(r)
					obs = mono_areas[mn][:,monolayer_residues[mn][actual_r]]
					if len(obs)>0:
						alone = True if len(monolayer_residues[1-mn][actual_r])==0 else False
						if not alone: ax.bar(counters[axnum]+rr+0.5*mn,np.mean(obs)-np.mean(mono_areas[mn]),width=0.5,align='edge',
							bottom=np.mean(mono_areas[mn]),lw=0,color=colorize(work.meta,resname=r),alpha=1,zorder=2)
						else: ax.bar(counters[axnum]+rr,np.mean(obs)-np.mean(mono_areas[mn]),width=1,
							bottom=np.mean(mono_areas[mn]),lw=0,color=colorize(work.meta,resname=r),alpha=1,zorder=2,align='edge')
			tick_positions.append((counters[axnum]+(ntypes+1)/2.))
			counters[axnum] += ntypes+1
			ax.set_title(composition+('\nwith cholesterol' if has_chol else ''))
			ax.set_ylabel(r'area per lipid $\mathrm{({nm}^{2})}$')
		labels = [work.meta[sn]['ion_label']+withchar+work.meta[sn]['ptdins_label'] for sn in sns_sub]
		ax.set_xticklabels(labels,rotation=90,ha='center')
		ax.set_xticks(tick_positions)
		for aa,ax in enumerate(axes[:2]): ax.set_xlim(0,counters[aa]+1)	
	picturesave('fig.%s'%plotname,work.plotdir,backup=False,
		version=True,meta={'print_version':print_version})
	plt.close()

if 'simple_view' in routine:

	"""
	Simple bar plot of the physiological areas per lipid in 2D.
	"""

	fsbase = 14
	has_chol = False
	shrink = 0.038 #---correct lines overhanging bars
	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.4")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.5,boxstyle="round,pad=0.6")
	allres = work.vars['selectors']['resnames_lipid']
	#---re-order by cation
	cations_order = ['NA','K','MG','Cal']
	#! updated for new omnicalc here
	sns = [sn for cat in cations_order for sn in work.metadata.collections['all']
		if work.meta[sn]['composition_name'] == 'asymmetric' and work.meta[sn]['cation']==cat]
	#---! overriding with canonical orderings
	sns = work.vars['orders']['canon']['asymmetric']
	if 'areadat' not in globals():
		areadat = {}
		for sn in sns:
			dat = data[1 if has_chol else 0][sn]['data']
			resnames = dat['resnames']
			reslist = [i for i in allres if i in unique(resnames)]
			ntypes = len(reslist)
			mono_areas = [dat['areas%d'%mn] for mn in range(2)]
			areadat[sn] = {}
			areadat[sn]['mean'] = np.mean(mono_areas[mn])
			#---! note that simulations with less than e.g. 800 frames will have different errorbar stats
			nsegs = int(mono_areas[mn].shape[0]/100.0)
			areadat[sn]['std'] = np.std([np.mean(mono_areas[mn][i*100:(i+1)*100]) for i in range(nsegs)])
	axes,fig = panelplot(figsize=(7,10),layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]})
	ax = axes[0]
	yvals = areadat
	colors = [colorize(work.meta[sn],comparison='protonation') for sn in sns]
	ax.bar(np.arange(len(sns)),[yvals[sn]['mean'] for sn in sns],color=colors,width=1.0,alpha=0.65,lw=0,align='edge')
	#---also plot standard deviations
	for ss,sn in enumerate(sns):
		width = 1.0
		dc,dc_std = yvals[sn]['mean'],yvals[sn]['std']
		std_lines = [([(ss+0.25)*width,(ss+0.75)*width],[dc-dc_std,dc-dc_std]),
			([(ss+0.25)*width,(ss+0.75)*width],[dc+dc_std,dc+dc_std]),
			([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
		for l in std_lines: ax.plot(l[0],l[1],color=colors[ss],lw=4,solid_capstyle='round',alpha=1,zorder=3)
	for snum,sn in enumerate(sns):
		ax.plot([snum+shrink,snum+1-shrink],[yvals[sn]['mean'],yvals[sn]['mean']],
			color=colors[snum],lw=5,solid_capstyle='round',alpha=1,zorder=3)
	mean_vals = [areadat[sn]['mean'] for sn in sns]
	ymin,ymax = np.mean(mean_vals)*0.9,np.mean(mean_vals)*1.05
	ax.set_ylim(ymin,ymax)
	ax.set_yticks(np.arange(np.round(ymin,2),np.round(ymax,2)+10**-2,10**-2))
	for snum,sn in enumerate(sns):
		tb = ax.text(snum+0.5,yvals[sn]['mean']+0.006,work.meta[sn]['ion_label'],
			bbox=tagbox_ion,ha='center',va='bottom',rotation=0,
			color='k',fontsize=fsbase+2)
		text = work.meta[sn]['ptdins_label']
		#text = re.sub('(3|4)',r"\large{\1}",text)
		tb = ax.text(snum+0.5,(ymax-ymin)*0.05+ymin,text,
			bbox=tagbox_ptdins,rotation=90,ha="center",va="bottom",color='k',fontsize=fsbase+4)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ax.set_yticklabels(['%.2f'%i for i in ax.get_yticks()],fontsize=fsbase+6)
	ax.set_ylabel(r'(projected) area per lipid $\mathrm{({nm}^{2})}$',fontsize=fsbase+8,labelpad=20)
	ax.set_xticklabels([])
	picturesave('fig.lipid_areas2d.simple.phys',work.plotdir,backup=False,version=True,meta={})
	plt.close()

if 'simple_view_newbars_dev' in routine:

	"""
	Simple bar plot of the physiological areas per lipid in 2D.
	MOVED TO figures.ppi_bar_stylized.
	"""

	#---aesthetic settings
	fsbase = 14
	trim = 0.15
	shrink = 0.038
	shrink = 0.05
	trim_std = 0.5-0.15

	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.4")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.5,boxstyle="round,pad=0.6")

	#---composition settings
	has_chol = False
	allres = work.vars['selectors']['resnames_lipid']
	#---use canonical orderings
	sns = work.vars['orders']['canon']['asymmetric']
	#---collect data
	if 'areadat' not in globals():
		areadat = {}
		for sn in sns:
			dat = data[1 if has_chol else 0][sn]['data']
			resnames = dat['resnames']
			reslist = [i for i in allres if i in unique(resnames)]
			ntypes = len(reslist)
			mono_areas = [dat['areas%d'%mn] for mn in range(2)]
			areadat[sn] = {}
			areadat[sn]['mean'] = np.mean(mono_areas[mn])
			#---! note that simulations with less than e.g. 800 frames will have different errorbar stats
			nsegs = int(mono_areas[mn].shape[0]/100.0)
			areadat[sn]['std'] = np.std([np.mean(mono_areas[mn][i*100:(i+1)*100]) for i in range(nsegs)])

	#---prepare the figure
	axes,fig = panelplot(figsize=(7,10),layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]})
	ax = axes[0]

	patterns = [ "/" , "\\" , "|" , "-" , "+" , "x", "o", "O", ".", "*" ][:len(sns)]

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
		'K':'beige',
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
		else: barspecs[sn]['color'] = 'w'
		for c in ['color','edgecolor'][1:]:
			barspecs[sn][c] = tuple(list(barspecs[sn][c])+[0.8])

	yvals = areadat
	xvals = np.arange(len(sns))+trim
	for yy,y in enumerate([yvals[sn]['mean'] for sn in sns]):
		x = xvals[yy]
		ax.bar(xvals[yy],y,
			width=1.0-trim*2,
			lw=3.0,align='edge',
			label=work.meta[sns[yy]]['cation'],
			**barspecs[sns[yy]])

	#---also plot standard deviations
	for ss,sn in enumerate(sns):
		width = 1.0
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
	#---hide the extra mean bars for now
	if False:
		for snum,sn in enumerate(sns):
			ax.plot([snum+shrink+trim,snum+1-shrink-trim],[yvals[sn]['mean'],yvals[sn]['mean']],
				color=barspecs[sn]['edgecolor'],lw=5,solid_capstyle='round',alpha=1,zorder=4,
				path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])

	mean_vals = [areadat[sn]['mean'] for sn in sns]
	ymin,ymax = np.mean(mean_vals)*0.9,np.mean(mean_vals)*1.05
	ax.set_ylim(ymin,ymax)
	ax.set_yticks(np.arange(np.round(ymin,2),np.round(ymax,2)+10**-2,10**-2))
	ax.set_xticklabels([work.meta[sn]['ptdins_label'] for sn in sns],rotation=90,fontsize=fsbase+6)
	ax.set_xticks(xvals+width/2.0)
	
	loc = 'upper right'
	title = None

	if False:
		loc = 'upper right'
		title = ['ions',None][-1]
		lw = 0
		leginds = [0,1,4,6]
		h,l = ax.get_legend_handles_labels()
		h = [mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,fc=barspecs[sn]['edgecolor']) for ss,sn in enumerate(sns)]
		legend = ax.legend([h[i] for i in leginds],[l[i] for i in leginds],
			loc=loc,fontsize=art['fs']['legend']+6,ncol=1,title=title)
		for legobj in legend.legendHandles: legobj.set_linewidth(lw)

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
		handleheight=2.0,markerscale=0.5)

	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ax.set_yticklabels(['%.2f'%i for i in ax.get_yticks()],fontsize=fsbase+10)
	ax.set_ylabel('(projected) area per lipid $\mathrm{({nm}^{2})}$',fontsize=fsbase+8,labelpad=20)

	picturesave('fig.lipid_areas2d.simple2.phys',work.plotdir,backup=False,
		version=True,meta={},pdf=True,extras=[legend])
	plt.close()

if 'simple_view_newbars' in routine:

	has_chol = False
	allres = work.vars['selectors']['resnames_lipid']
	#---use canonical orderings
	sns = work.vars['orders']['canon']['asymmetric']
	#---collect data
	if 'areadat' not in globals():
		areadat = {}
		for sn in sns:
			dat = data[1 if has_chol else 0][sn]['data']
			resnames = dat['resnames']
			reslist = [i for i in allres if i in unique(resnames)]
			ntypes = len(reslist)
			mono_areas = [dat['areas%d'%mn] for mn in range(2)]
			areadat[sn] = {}
			areadat[sn]['mean'] = np.mean(mono_areas[mn])
			#---! note that simulations with less than e.g. 800 frames will have different errorbar stats
			nsegs = int(mono_areas[mn].shape[0]/100.0)
			areadat[sn]['std'] = np.std([np.mean(mono_areas[mn][i*100:(i+1)*100]) for i in range(nsegs)])

	ppi_bar_stylized(
		name='fig.lipid_areas2d.simple2.phys',
		data=areadat,
		sns=sns,
		zero_low=True,
		ylabel='(projected) area per lipid $\mathrm{({nm}^{2})}$',
		)
