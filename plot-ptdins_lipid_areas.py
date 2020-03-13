#!/usr/bin/env python

"""
Plot bilayer areas for the ptdins project.
"""

#---autoplot settings
plotrun.routine = None
global seepspace

@autoload(plotrun)
def loader():
	"""Load data."""
	#---only load once
	if 'data' not in globals():
		data,calc = plotload(plotname)
		sns = work.sns()
		#---autoplot exports to globals
		global seepspace
		seepspace = 'data,sns'.split(',')
		for key in seepspace: globals()[key] = locals()[key]

@autoplot(plotrun)
def make_plots():
	"""
	Plot summary and comprehensive plots.
	"""
	if True: 

		if 0:
			# added v599 in sensible places now that it is ready
			plot_areas(out_fn='lipid_areas',lims2d=(0.65,0.7), #! (195,210), lims now for area per lipid
				sns2d=['membrane-v%3d'%i for i in [534,532,599,533,531,536,530,538]],
				sns3d=['membrane-v%3d'%i for i in [599,536,534,532,533,531,530,538]])

		# on 2019.09.11 replotting
		else:
			#! (195,210), lims now for area per lipid
			
			if 0:
				plot_areas2d(out_fn='lipid_areas_2d',lims2d=(0.65,0.7), 
					sns2d=['membrane-v%3d'%i for i in [534,532,599,533,531,536,530,538]],)
				plot_areas3d(out_fn='lipid_areas_3d',lims2d=(0.65,0.7), 
					sns3d=['membrane-v%3d'%i for i in [599,536,534,532,533,531,530,538]],)
			#! revised on 2020.03.12
			else:
				plot_areas2d(out_fn='lipid_areas_2d',lims2d=(0.65,0.7), 
					sns2d=['membrane-v%3d'%i for i in [534,532,599,533,531,538]],)
				plot_areas3d(out_fn='lipid_areas_3d',lims2d=(0.65,0.7), 
					sns3d=['membrane-v%3d'%i for i in [599,534,532,533,531,538]],)
	else:
		sns_symmetric = ['membrane-v%3d'%i for i in [509,514,515,543,510,511]]
		sns_symmetric3d = ['membrane-v%3d'%i for i in [509,514,515,543,510,511]]
		plot_areas(out_fn='lipid_areas.symmetric',labels3d=True,
			sns2d=sns_symmetric,sns3d=sns_symmetric3d,lims2d=(240,265))

def plot_areas(out_fn,sns2d,sns3d,lims2d=None,labels3d=False):
	"""
	Plot 2D and 3D lipid areas.
	"""
	plotspec = ptdins_manuscript_settings()
	patches = []
	#---include broken axes in the layout
	axes,fig = panelplot(figsize=(12,10),
		layout={'out':{'grid':[1,2],'wspace':0.4},
		'ins':[{'grid':[3,1],'hratios':[10,1,4],'hspace':0.1},{'grid':[1,1]}],})
	sns_this = sns2d
	#---repeated plotting for broken axis
	for dup_num,ax in enumerate(axes[0][:2]):
		for snum,sn in enumerate(sns_this):
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			dat = data['lipid_areas2d'][sn]['data']
			#! hacking to per-lipid areas
			norm_factor = 300.0 #! for cholesterol systems
			#---! note that 2D areas should be equal between leaflets up to minor approximation errors
			area2d = dat['areas%d'%top_mono].sum(axis=1).mean()/norm_factor
			area2d_err = dat['areas%d'%top_mono].sum(axis=1).std()/norm_factor
			ax.bar([snum],[area2d],width=1.0,lw=0,edgecolor='w',
				color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
				hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
			ax.errorbar([snum],[area2d],yerr=[area2d_err],alpha=1.0,lw=4.0,c='k')
		if dup_num==0 and lims2d!=None: ax.set_ylim(lims2d)
	#---broken axis requires duplicates
	axdup = [axes[0][0],axes[0][1]]
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
		'bot_r':[(-d*rx,d*rx),(1-d*ry,1+d*ry)],}
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
	axdup[0].set_xticks([])
	axdup[1].set_xticks([])
	for ax in axdup: 
		ax.set_yticklabels(['%.02f'%i for i in ax.get_yticks()],fontsize=16)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	axdup[0].set_ylabel('area per lipid (2D) $\mathrm{({nm}^{2})}$',fontsize=20)
	#---end broken axis modifications
	#---plot 3D areas
	ax = axes[1][0]
	sns_this = sns3d
	for snum,sn in enumerate(sns_this):
		top_mono = work.meta[sn].get('index_top_monolayer',0)
		bot_mono = {0:1,1:0}[top_mono]
		dat = data['lipid_areas3d'][sn]['data']
		#! hacking to per-lipid areas
		norm_factor = 300.0 #! for cholesterol systems
		top_area = dat['areas%d'%top_mono].sum(axis=1).mean()/norm_factor
		bot_area = dat['areas%d'%bot_mono].sum(axis=1).mean()/norm_factor
		ax.bar([snum],[top_area-bot_area],bottom=[bot_area],width=1.0,lw=0,edgecolor='w',
			color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
			hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
		yerr_top = dat['areas%d'%top_mono].sum(axis=1).std()/norm_factor
		yerr_bot = dat['areas%d'%bot_mono].sum(axis=1).std()/norm_factor
		ax.errorbar([snum],[top_area],yerr=[yerr_top],alpha=1.0,lw=4.0,c='k')
		ax.errorbar([snum],[bot_area],yerr=[yerr_bot],alpha=1.0,lw=4.0,c='k')
	#---aesthetics
	ax.set_xticks([])
	ax.set_yticklabels(['%.02f'%i for i in ax.get_yticks()],fontsize=16)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.set_ylabel('area per lipid by leaflet (3D) $\mathrm{({nm}^{2})}$',fontsize=20)
	ax.set_xlim((-0.5-0.25,len(sns_this)+3.0)) # increased padding from 2.5 to 3 for labels
	#---symmetric has bars that are too small so we add extra labels
	if labels3d:
		tagbox_ptdins = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
		ax.set_xticks([])
		ylevel = 270.0
		xcenters = np.arange(len(sns3d))
		for snum,sn in enumerate(sns3d):
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			bot_mono = {0:1,1:0}[top_mono]
			dat = data['lipid_areas3d'][sn]['data']
			bot_area = dat['areas%d'%bot_mono].sum(axis=1).mean()
			bot_std = dat['areas%d'%bot_mono].sum(axis=1).std()
			text = '%s %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'])
			#---custom offset below
			tb = ax.text(1.0*(xcenters[snum]),266.5,text,
				bbox=tagbox_ptdins,rotation=-90,ha="center",va="top",color='k',fontsize=12)
			patches.append(tb)
			#---custom ylim
			ax.set_ylim((265,290))
			ax.axvline(snum,lw=1,c='k',alpha=0.5,zorder=1)
	#---annotations after the last bar (adapted from new charging curve plots)
	el = mpl.patches.Ellipse((2,-1),0.5,0.5)
	for yy in range(2):
		color = '#bdbdbd'
		ann = ax.annotate(['inner\nleaflet','outer\nleaflet'][yy],
			xy=(snum+0.5,[top_area,bot_area][yy]),xycoords='data',
			xytext=(35,0),textcoords='offset points',
			size=12,va="center",bbox=dict(boxstyle="round",fc=color,ec="none"),
			arrowprops=dict(arrowstyle="wedge,tail_width=1.",
				fc=color,ec="none",patchA=None,patchB=el,relpos=(0.2,0.5),))
		patches.append(ann)
	#---fancy legend sequence (repeated in a few places)
	bar_formats = make_bar_formats(sns_this,work=work)
	comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] for sn in sns_this])),
		ion_names=list(set([work.meta[sn]['cation'] for sn in sns_this])))
	#---legend below
	kwargs = dict(bbox=(0.5,-0.1),loc='upper center',ncol=2)
	legend,patches = legend_maker_stylized(axes[0][1],work=work,
		sns_this=sns_this,bar_formats=bar_formats,comparison_spec=comparison_spec,fancy=False,**kwargs)
	patches.append(legend)
	fig.delaxes(axes[0][2])
	meta = {'sns':sns3d}
	picturesave('fig.%s'%out_fn,
		work.plotdir,backup=False,version=True,meta=meta,extras=patches,form='pdf')

## HACKING VIA COPY on 2019.09.11 for expedience

def plot_areas2d(out_fn,sns2d,lims2d=None,labels3d=False):
	"""
	Plot 2D and 3D lipid areas.
	"""
	plotspec = ptdins_manuscript_settings()
	patches = []
	#---include broken axes in the layout
	axes,fig = panelplot(figsize=(5,10),
		layout={'out':{'grid':[1,1],'wspace':0.4},
		'ins':[{'grid':[3,1],'hratios':[10,1,4],'hspace':0.1}],})
	axes = [axes]
	sns_this = sns2d
	#---repeated plotting for broken axis
	for dup_num,ax in enumerate(axes[0][:2]):
		for snum,sn in enumerate(sns_this):
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			dat = data['lipid_areas2d'][sn]['data']
			#! hacking to per-lipid areas
			norm_factor = 300.0 #! for cholesterol systems
			#---! note that 2D areas should be equal between leaflets up to minor approximation errors
			area2d = dat['areas%d'%top_mono].sum(axis=1).mean()/norm_factor
			area2d_err = dat['areas%d'%top_mono].sum(axis=1).std()/norm_factor
			ax.bar([snum],[area2d],width=1.0,lw=0,edgecolor='w',
				color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
				hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
			ax.errorbar([snum],[area2d],yerr=[area2d_err],alpha=1.0,lw=4.0,c='k')
		if dup_num==0 and lims2d!=None: ax.set_ylim(lims2d)
	#---broken axis requires duplicates
	axdup = [axes[0][0],axes[0][1]]
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
		'bot_r':[(-d*rx,d*rx),(1-d*ry,1+d*ry)],}
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
	axdup[0].set_xticks([])
	axdup[1].set_xticks([])
	for ax in axdup: 
		ax.set_yticklabels(['%.02f'%i for i in ax.get_yticks()],fontsize=16)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	axdup[0].set_ylabel('area per lipid (2D) $\mathrm{({nm}^{2})}$',fontsize=20)
	#---end broken axis modifications
	#---fancy legend sequence (repeated in a few places)
	bar_formats = make_bar_formats(sns_this,work=work)
	comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] for sn in sns_this])),
		ion_names=list(set([work.meta[sn]['cation'] for sn in sns_this])))
	#---legend below
	kwargs = dict(bbox=(0.5,-0.1),loc='upper center',ncol=2)
	legend,patches = legend_maker_stylized(axes[0][1],work=work,
		sns_this=sns_this,bar_formats=bar_formats,comparison_spec=comparison_spec,fancy=False,**kwargs)
	patches.append(legend)
	fig.delaxes(axes[0][2])
	meta = {'sns':sns2d}
	picturesave('fig.%s'%out_fn,
		work.plotdir,backup=False,version=True,meta=meta,extras=patches,form='pdf')

def plot_areas3d(out_fn,sns3d,lims2d=None,labels3d=False):
	"""
	Plot 2D and 3D lipid areas.
	"""
	plotspec = ptdins_manuscript_settings()
	patches = []
	#---include broken axes in the layout
	axes,fig = panelplot(figsize=(6,10),
		layout={'out':{'grid':[1,1],'wspace':0.4},
		'ins':[{'grid':[1,1]}],})
	axes = [axes]
	ax = axes[0][0]
	sns_this = sns3d
	for snum,sn in enumerate(sns_this):
		top_mono = work.meta[sn].get('index_top_monolayer',0)
		bot_mono = {0:1,1:0}[top_mono]
		dat = data['lipid_areas3d'][sn]['data']
		#! hacking to per-lipid areas
		norm_factor = 300.0 #! for cholesterol systems
		top_area = dat['areas%d'%top_mono].sum(axis=1).mean()/norm_factor
		bot_area = dat['areas%d'%bot_mono].sum(axis=1).mean()/norm_factor
		ax.bar([snum],[top_area-bot_area],bottom=[bot_area],width=1.0,lw=0,edgecolor='w',
			color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
			hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
		yerr_top = dat['areas%d'%top_mono].sum(axis=1).std()/norm_factor
		yerr_bot = dat['areas%d'%bot_mono].sum(axis=1).std()/norm_factor
		ax.errorbar([snum],[top_area],yerr=[yerr_top],alpha=1.0,lw=4.0,c='k')
		ax.errorbar([snum],[bot_area],yerr=[yerr_bot],alpha=1.0,lw=4.0,c='k')
	#---aesthetics
	ax.set_xticks([])
	ax.set_yticklabels(['%.02f'%i for i in ax.get_yticks()],fontsize=16)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.set_ylabel('area per lipid by leaflet (3D) $\mathrm{({nm}^{2})}$',fontsize=20)
	ax.set_xlim((-0.5-0.25,len(sns_this)+2.0)) # increased padding from 2.5 to 3 for labels
	#---symmetric has bars that are too small so we add extra labels
	if labels3d:
		tagbox_ptdins = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
		ax.set_xticks([])
		ylevel = 270.0
		xcenters = np.arange(len(sns3d))
		for snum,sn in enumerate(sns3d):
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			bot_mono = {0:1,1:0}[top_mono]
			dat = data['lipid_areas3d'][sn]['data']
			bot_area = dat['areas%d'%bot_mono].sum(axis=1).mean()
			bot_std = dat['areas%d'%bot_mono].sum(axis=1).std()
			text = '%s %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'])
			#---custom offset below
			tb = ax.text(1.0*(xcenters[snum]),266.5,text,
				bbox=tagbox_ptdins,rotation=-90,ha="center",va="top",color='k',fontsize=12)
			patches.append(tb)
			#---custom ylim
			ax.set_ylim((265,290))
			ax.axvline(snum,lw=1,c='k',alpha=0.5,zorder=1)
	#---annotations after the last bar (adapted from new charging curve plots)
	el = mpl.patches.Ellipse((2,-1),0.5,0.5)
	for yy in range(2):
		color = '#bdbdbd'
		ann = ax.annotate(['inner\nleaflet','outer\nleaflet'][yy],
			xy=(snum+0.5,[top_area,bot_area][yy]),xycoords='data',
			xytext=(35,0),textcoords='offset points',
			size=12,va="center",bbox=dict(boxstyle="round",fc=color,ec="none"),
			arrowprops=dict(arrowstyle="wedge,tail_width=1.",
				fc=color,ec="none",patchA=None,patchB=el,relpos=(0.2,0.5),))
		patches.append(ann)
	#---fancy legend sequence (repeated in a few places)
	bar_formats = make_bar_formats(sns_this,work=work)
	comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] for sn in sns_this])),
		ion_names=list(set([work.meta[sn]['cation'] for sn in sns_this])))
	#---legend below
	kwargs = dict(bbox=(0.5,-0.1),loc='upper center',ncol=4)
	legend,patches = legend_maker_stylized(axes[0][0],work=work,
		sns_this=sns_this,bar_formats=bar_formats,comparison_spec=comparison_spec,fancy=False,**kwargs)
	patches.append(legend)
	meta = {'sns3d':sns3d}
	picturesave('fig.%s'%out_fn,
		work.plotdir,backup=False,version=True,meta=meta,extras=patches,form='pdf')

plotrun.routine = None

if __name__=='__main__':

	if 1:

		art = {'fs':{'label':16,'title':20}}
		mpl.rcParams['hatch.linewidth'] = 1.5
		kb = 1.38064852*10**-23 # m2 kg s-2 K-1 or J/K or N-m/K
		factor = 2*kb*310*(10**-9/(10**-9)**2)**2*10**3 # nm to meters in two terms, then from N/m to mN/m
		unit_factor = 1000.
		#! note that values are 2.5x10^-7 mN/m
		"""
		see Waheed and Edholm: K_A = 2 a kBT / (sig_a)^2 / N
		a is the area per lipid, or 2A/N
		sig_a is the mean-squared fluctuations in the area per lipid
		use factor above for N/m or just 2.0 for kBT
		NOTE. you get way higher values when you use lipid_com or lipid_chol_com
		instead, it is best to use lipid_chol_com_headless, which gives experimental values
		this makes sense since cholesterol is a molecule and the headgroups can really screw up the fluctuations
		I would imagine the increase in compressibility is because the headgroup decreases the variance in lipid 
		areas, which is in the denominator
		recall the design decisions: cholesterol is a lipid, use tails (i.e. headless), 
		then average over time, then lipids
		"""
		# use a formal for total area and area fluctuation i.e. without including N_lipids. !! this is way off
		ka_func = lambda x: (factor*x.sum(axis=1).mean()/area2d.sum(axis=1).std()**2*unit_factor).mean()
		# several ways to do this. average time and lipid together
		ka_func = lambda x: (factor*x.mean()/area2d.std()**2/N_lipids*unit_factor)
		# several ways to do this. average time then by lipid
		#   this means we are using the timewise fluctuations, and then collecting statistics over lipids
		#   furthermore I think this more closely matches 
		ka_func = lambda x: (factor*x.mean(axis=0)/x.std(axis=0)**2*unit_factor).mean()
		# trying time on the oustide, lipids on the inside since MSD is typically over lipids, instantaneous t
		ka_func = lambda x: (factor*x.mean(axis=1)/x.std(axis=1)**2/N_lipids*unit_factor).mean()
		# final decision is the same as above, for obvious reasons. hence order doesn't matter and we 
		#   treat time and lipid-wise deviations distinctly. note that this formula is 
		#   best explained in the intro to Waheed and Edholm 2009
		ka_func = lambda x: (factor*x.mean(axis=0)/x.std(axis=0)**2/N_lipids*unit_factor).mean()

		#! either cholesterol is included or not
		area_compressibility = {}
		for sn in sns:
			dat = data['lipid_areas2d'][sn]['data']
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			area2d = dat['areas%d'%top_mono]
			if sn=='membrane-v532':
				# [STATUS] variance for v532 is 0.0663 # when lipid_chol_com
				# [STATUS] variance for v532 is 0.0931 # when headless. pattern is the same
				status('variance for v532 is %.4f'%area2d.std(axis=0).mean())
			N_lipids = area2d.shape[1]*2 # multiply by 2 for number-symmetric bilayers for total lipid count
			area_compressibility[sn] = dict(mean=ka_func(area2d))
			# easy way to get error bars
			nsegs = 10.
			subd = ((np.arange(len(area2d))/(len(area2d)/nsegs)).astype(int))
			if True:
				area_compressibility[sn]['std'] = np.std([
					ka_func(area2d[np.where(subd==i)]) for i in np.unique(subd).astype(int)])
			else: #! probably incorrect
				x = area2d
				area_compressibility[sn]['std'] = \
					(factor*x.std(axis=0)/area2d.std(axis=0)**2/N_lipids*unit_factor).mean()

		if 0: sns_collect = [['membrane-v%d'%i for i in j] for j in [
			[509,514,515,510,511],
			[538,536,530,531,533,599,532,534],]]
		sns_collect = [['membrane-v%d'%i for i in j] for j in [
			[510,511],
			[538,531,533,599,532,534],]]
		axes,fig = panelplot(figsize=(12,10),
			layout={'out':{'grid':[1,1]},
			'ins':[{'grid':[1,2],'wratios':[len(i) for i in sns_collect],'wspace':0.3}]})
		for pnum,sns_this in enumerate(sns_collect):
			ax = axes[pnum]
			plotspec = ptdins_manuscript_settings()
			for snum,sn in enumerate(sns_this):
				ax.bar([snum],[area_compressibility[sn]['mean']],width=1.0,lw=0,edgecolor='w',
					color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
					hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
				ax.errorbar([snum],[area_compressibility[sn]['mean']],yerr=[area_compressibility[sn]['std']],
					alpha=1.0,lw=4.0,c='k')
			ax.set_xticks([])
			ax.set_ylabel(r'$K_A\,(dyn/cm)$',fontsize=art['fs']['label'])
			ax.tick_params(axis='y',which='both',left='off',right='off',
				labelleft='on',labelsize=art['fs']['label'])	
			if pnum==0: ax.set_title('symmetric',fontsize=art['fs']['title'])
			elif pnum==1: ax.set_title('physiological',fontsize=art['fs']['title'])
			else: raise Exception
			if pnum==0:
				sns_all = list(set([i for j in sns_collect for i in j]))
				bar_formats = make_bar_formats(sns_all,work=work)
				# extra hatching for legend to make the lines smaller because too thick
				for sn in bar_formats: 
					if bar_formats[sn].get('hatch'):
						# only the square one looks bad so we fix it
						if '+' in bar_formats[sn]['hatch']: 
							bar_formats[sn]['hatch'] = bar_formats[sn]['hatch'][:-1]
				# legend_maker_stylized takes 'bbox' for 'bbox_to_anchor'
				# custom position for legend so it is beneath both axes
				kwargs = dict(loc='center left',ncol=5,bbox=(-0.15,-0.1))
				comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] 
					for sn in sns_all])),
					ion_names=list(set([work.meta[sn]['cation'] for sn in sns_all])))
				legend,patches = legend_maker_stylized(ax,work=work,
					sns_this=sns_all,bar_formats=bar_formats,
					comparison_spec=comparison_spec,fancy=False,**kwargs)
		picturesave('fig.area_compressibility',
			work.plotdir,backup=False,version=True,meta={},extras=[legend],form='pdf')
