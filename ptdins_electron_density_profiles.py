#!/usr/bin/env python

# normalize to bulk which starts 10nm from the midplane (conservative)
bilayer_width = 10.

@autoload(plotrun)
def load():
	data = plotload('ptdins_electron_density_profiles')
	sns = work.sns()
	kwargs = dict(sns=sns,data=data,bilayer_width=bilayer_width)
	post = post_process_edps(**kwargs)

def post_process_edps(sns,data,**kwargs):
	bilayer_width = kwargs['bilayer_width']
	post = {}
	for sn in sns:
		dat = data.this[sn]
		try: bin_size = dat['bin_size']
		except:
			dat = data.this[sn]['data'] #! no idea why we had to add this
			bin_size = dat['bin_size']
		nbins = dat['nbins']
		vecs = dat['vecs']
		half_box_height = vecs.mean(axis=0)[2]/2.
		bin_size_actual = half_box_height/nbins
		groups = list(dat['resnames'])+list(dat['group_regexes'])
		post[sn] = {}
		for gnum,group in enumerate(groups):
			xvals = np.linspace(-1.*half_box_height,half_box_height,nbins)
			yvals = dat['tabulated'][:,gnum].mean(axis=0)
			name = {u'.+':'all',u'^(OW)|(HW(1|2))$':'water', u'^C[0-9]+':'carbon'}.get(group,group)
			bulk_inds = np.where(np.abs(xvals)>bilayer_width)[0]
			norm_factor = np.product(vecs.mean(axis=0)[:2])*bin_size
			post[sn][gnum] = {'xvals':xvals,'yvals':yvals,
				'n':len(dat['group_%d'%gnum]),'name':name,'norm_factor':norm_factor}
	return post

@autoplot(plotrun)
def plot_summary():
	figsize = (12,12)
	axes,fig = square_tiles(ntiles=len(sns),figsize=figsize,hspace=0.35,wspace=0.35)
	for snum,sn in enumerate(sns):
		ax = axes[snum]
		for key in sorted(post[sn].keys()):
			xvals,yvals = [post[sn][key][k] for k in ['xvals','yvals']]
			ax.plot(xvals,yvals/post[sn][key]['norm_factor'],label=post[sn][key]['name'])
		ax.set_ylabel('electron density $(e\,{nm}^{-3})$')
		ax.set_xlabel('z (nm)')
		ax.tick_params(axis='both',which='both',length=0)
		ax.set_title('%s and %s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label']))
	picturesave('fig.electron_density_profiles.summary',work.plotdir,
		backup=False,version=True,meta={},extras=[])
	plt.close()

@autoplot(plotrun)
def plot_interesting(**kwargs):
	"""
	Plot a subset of the EDPs. Standalone or called by plot-ptdins_ion_binding.py
	"""
	global layout_style
	extras = []
	layout_style = kwargs.get('layout_style',None)
	layout_styles = kwargs.get('layout_styles',['by_sn','by_ion'])
	do_external = kwargs.get('do_external',False)
	mark_peaks = kwargs.get('mark_peaks',True)
	xlim = kwargs.get('xlim',False)
	early_color = kwargs.get('early_color','k')
	ion_names = ['cation','anion']
	fs_labels = kwargs.get('fs_labels',10)
	fs_legend = kwargs.get('fs_legend',10)
	show_title = kwargs.get('show_title',True)
	nsegs = kwargs.get('nsegs',3)
	do_interpolation = False
	lw = 2.0

	def get_panel(sn,ion_name):
		global layout_style
		if layout_style=='by_sn': return sns_this.index(sn)
		elif layout_style=='by_ion': return ion_names.index(ion_name)
		else: raise Exception
	def get_color(sn,ion_name):
		global layout_style
		if layout_style=='by_sn': 
			return {'cation':'r','anion':'b'}[ion_name]
		elif layout_style=='by_ion': 
			return colorize(work.meta[sn],comparison='protonation')
		else: raise Exception
	def get_label(sn,ion_name):
		global layout_style
		if layout_style=='by_sn': return (
			work.meta[sn]['ion_label'] if ion_name=='cation' else '${Cl}^{-}$')
		elif layout_style=='by_ion': 
			return '%s and %s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label'])
		else: raise Exception
	get_color = kwargs.get('get_color',get_color)

	def plot_edp(axes,get_panel,sns_this):
		for snum,sn in enumerate(sns_this):
			legend_labels,legend_patches = [],[]
			for ion_num,ion_name in enumerate(ion_names):
				dat = data.this[sn]
				try: vecs = dat['vecs']
				except:
					dat = data.this[sn]['data'] #! had to add this
					vecs = dat['vecs']
				bin_size = dat['bin_size']
				norm_factor = np.product(vecs.mean(axis=0)[:2])*bin_size
				nbins = dat['nbins']
				vecs = dat['vecs']
				half_box_height = vecs.mean(axis=0)[2]/2.
				bin_size_actual = half_box_height/nbins
				groups = list(dat['resnames'])+list(dat['group_regexes'])
				gnum = groups.index(work.meta[sn][ion_name])
				xvals = np.linspace(-1.*half_box_height,half_box_height,nbins)
				nframes = dat['tabulated'][:,gnum].shape[0]
				yvals = np.array_split(dat['tabulated'][:,gnum],nsegs)
				half_box_round = float(np.floor(half_box_height/bin_size)*bin_size)
				ax = axes[get_panel(sn=sn,ion_name=ion_name)]
				ax.set_ylabel('electron density $(e\,{nm}^{-3})$',fontsize=fs_labels)
				ax.set_xlabel('z (nm)',fontsize=fs_labels)
				ax.tick_params(axis='both',which='both',length=0,labelsize=fs_labels)
				if xlim: ax.set_xlim(xlim)
				if layout_style=='by_sn' and show_title:
					ax.set_title('%s and %s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label']))
				elif layout_style=='by_ion' and show_title: ax.set_title(ion_name)
				elif not show_title: pass
				else: raise Exception
				cmap = mpl.colors.LinearSegmentedColormap.from_list(
					'calcium_shift',[early_color,get_color(sn=sn,ion_name=ion_name)])
				colors = [cmap(float(i)/(nsegs-1)) for i in range(nsegs)]
				for ynum,y in enumerate(yvals): 
					if do_interpolation:
						interp_f = scipy.interpolate.interp1d(xvals,y.mean(axis=0))
						y_interp = np.array([interp_f(i) 
							for i in np.arange(-half_box_round,half_box_round+bin_size,bin_size)])
						yvals_this = y_interp/norm_factor
					else: yvals_this = y.mean(axis=0)/norm_factor
					label = get_label(sn,ion_name)
					ax.plot(xvals,yvals_this,color=colors[ynum],lw=lw,
						label=label if ynum==len(yvals)-1 else None,zorder=10 if ion_name=='cation' else 9)
					if ynum==len(yvals)-1 and mark_peaks:
						peak_pos = xvals[xvals>0][yvals_this[xvals>0].argmax()]
						ax.axvline(peak_pos,lw=1,c='k')
					# sequential lines are too busy so we just code in the gradient
					if ynum==len(yvals)-1:
						legend_labels.append(label if ynum==len(yvals)-1 else None)
						legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,
							facecolor=colors[ynum],edgecolor='w',lw=2))
				#! note legend might be wrong for the EDP singleton plots because it ignores layout_style?
				if ion_num==len(ion_names)-1: 
					legend_labels.append('early')
					legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,
						facecolor=early_color,edgecolor='w',lw=2))
				legend = ax.legend(legend_patches,legend_labels,fontsize=fs_legend,
					labelspacing=0.,handleheight=1.5,markerscale=0.5)
				frame = legend.get_frame()
				frame.set_edgecolor('w')
				frame.set_facecolor('w')
				extras.append(legend)

	# famous original main loop
	if not do_external:
		for layout_style in layout_styles:
			if do_interpolation: import scipy
			sns_this = ['membrane-v531','membrane-v532']
			figsize = (6,8)
			extras = []
			axes,fig = square_tiles(ntiles=len(sns_this),figsize=figsize,hspace=0.4)
			plot_edp(axes,get_panel)
			picturesave('fig.electron_density_profiles.interesting.%s'%layout_style,work.plotdir,
				backup=False,version=True,meta={},extras=[])
			plt.close()
	else: 
		plot_edp(kwargs['axes'],kwargs['get_panel'],kwargs['sns_this'])
	return extras

if __name__=='__main__': pass
