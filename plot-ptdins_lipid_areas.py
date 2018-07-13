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
		plot_areas(out_fn='lipid_areas',lims2d=(0.65,0.7), #! (195,210), lims now for area per lipid
			sns2d=['membrane-v%3d'%i for i in [534,532,533,531,536,530,538]],
			sns3d=['membrane-v%3d'%i for i in [536,534,532,533,531,530,538]])
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
	ax.set_xlim((-0.5-0.25,len(sns_this)+2.5))
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
	picturesave('fig.%s'%out_fn,
		work.plotdir,backup=False,version=True,meta={},extras=patches)

plotrun.routine = None #[]
if __name__=='__main__':

	kb = 1.38064852*10**-23 # m2 kg s-2 K-1 or J/K or N-m/K
	factor = 2*kb*300*(10**-9/(10**-9)**2)**2*10**3 # nm to meters in two terms, then from N/m to mN/m
	#! note that values are 2.5x10^-7 mN/m

	area_compressibility = {}
	for sn in sns:
		dat = data['lipid_areas2d'][sn]['data']
		top_mono = work.meta[sn].get('index_top_monolayer',0)
		area2d = dat['areas%d'%top_mono]
		area_compressibility[sn] = factor*area2d.mean()/area2d.std()**2

	axes,fig = panelplot(figsize=(12,10),
		layout={'out':{'grid':[1,1],'wspace':0.4},
		'ins':[{'grid':[1,2]}]})

	sns_collect = [['membrane-v%d'%i for i in j] for j in [
		[536,538,530,531,532,533,534],
		[509,514,515,510,511]]]
	for pnum,sns_this in enumerate(sns_collect):
		ax = axes[pnum]
		plotspec = ptdins_manuscript_settings()
		for snum,sn in enumerate(sns_this):
			ax.bar([snum],[area_compressibility[sn]],width=1.0,lw=0,edgecolor='w',
				color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
				hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
		ax.set_xticks([])
		ax.set_ylabel('$K_A\,(k_B T)??$')
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')	
		if pnum==0: ax.set_title('asymmetric')
		elif pnum==1: ax.set_title('symmetric')
		else: raise Exception
		if pnum==0:
			sns_all = list(set([i for j in sns_collect for i in j]))
			bar_formats = make_bar_formats(sns_all,work=work)
			kwargs = dict(bbox=(0.5,-0.1),loc='upper center',ncol=3)
			comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] 
				for sn in sns_all])),
				ion_names=list(set([work.meta[sn]['cation'] for sn in sns_all])))
			legend,patches = legend_maker_stylized(ax,work=work,
				sns_this=sns_all,bar_formats=bar_formats,
				comparison_spec=comparison_spec,fancy=False,**kwargs)
	picturesave('fig.area_compressibility',
		work.plotdir,backup=False,version=True,meta={},extras=[legend])
