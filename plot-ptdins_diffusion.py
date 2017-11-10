#!/usr/bin/env python

"""
Lipid diffusion plots.
Only change sections with "user additions here".
"""

import brewer2mpl
from codes.looptools import basic_compute_loop

#---settings
cutoff_time = [[1,20],[8,72]][0]
fit_type = 'linear'

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
	q,r = np.polyfit(xvals[ins],yvals[ins],1)
	return q/factor,r/factor

@autoplot(plotrun)
def diffusion_plots():
	"""Plot different versions of the diffusion plot."""
	#---simple version for main text
	diffusion_plot(out_fn='diffusion_lipids',hatch_lw=2,
		sns_this=[i for i in work.specs['collections']['asymmetric_all'] if work.meta[i]['cation']!="Na,Cal"])
	#---comprehensive version for supplement with a space between groups
	sns_l = ['membrane-v%3d'%i for i in [509,514,515,543,510,511]]
	sns_r = ['membrane-v%3d'%i for i in [542,536,538,530,531,533,532,534]]
	xtick_details = dict(xticks=np.cumsum([np.arange(len(i)).mean()+o 
		for i,o in zip([sns_l,sns_r],[0,len(sns_l)-1.0])]),
		xtick_labels=['symmetric','asymmetric'],
		bar_x=np.concatenate((np.arange(len(sns_l)),np.arange(len(sns_r))+1.0+len(sns_l))))
	diffusion_plot(out_fn='diffusion_lipids.comprehensive',
		sns_this=sns_l+sns_r,hatch_lw=1,xtick_details=xtick_details)

def diffusion_plot(sns_this,out_fn,spacer=None,hatch_lw=1.0,xtick_details=None):
	"""
	Plot diffusion rates.
	"""
	if xtick_details==None: xtick_details = {}
	#---need thicker hatchess
	#---! cannot change hatch frequency except via figure size and DPI
	mpl.rcParams['hatch.linewidth'] = hatch_lw
	plotspec = ptdins_manuscript_settings()
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
	axes,fig = panelplot(figsize=(8,10),
		layout={'out':{'grid':[1,1],'wspace':0.0,'hspace':0.3},
		'ins':[{'grid':[1,1]} for i in range(1)],})
	ax = axes[0]
	for snum,sn in enumerate(sns_this):
		xpos = xtick_details.get('bar_x',np.arange(len(sns_this)))[snum]
		ax.bar([xpos],diffusion_consts[sn]['mean'],width=1.0,lw=0,edgecolor='w',
			color=plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]],
			hatch=plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']],**(
				{'label':work.meta[sn]['ion_label']} if True else {}))
		ax.errorbar([xpos],[diffusion_consts[sn]['mean']],yerr=[diffusion_consts[sn]['std']],
			alpha=1.0,lw=4.0,c='k')
	#---labels
	ax.set_title('lipid diffusion',fontsize=18)
	ax.set_ylabel('$\mathrm{D\,({\mu m}^{2} {s}^{-1}})$',fontsize=18)
	ax.set_yticklabels(ax.get_yticks(),fontsize=18)
	if xtick_details:
		ax.set_xticks(xtick_details['xticks'])
		ax.set_xticklabels(xtick_details['xtick_labels'],fontsize=18)
	else: ax.set_xticks([])
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labeltop='off')
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	#---fancy legend sequence (repeated in a few places)
	bar_formats = make_bar_formats(sns_this,work=work)
	comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] for sn in sns_this])),
		ion_names=list(set([work.meta[sn]['cation'] for sn in sns_this])))
	#---legend below
	kwargs = dict(bbox=(0.5,-0.05),loc='upper center',ncol=4)
	legend,patches = legend_maker_stylized(ax,work=work,
		sns_this=sns_this,bar_formats=bar_formats,comparison_spec=comparison_spec,
		**kwargs)
	picturesave('fig.%s'%out_fn,
		work.plotdir,backup=False,version=True,meta={},extras=[legend])

@autoload(plotrun)
def reload():
	"""Load everything for the plot only once."""
	#---user additions here. include a list of all of your globals and a once-only load sequence
	global data,sns,post
	data,calc = plotload(plotname)
	sns = work.sns()
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
