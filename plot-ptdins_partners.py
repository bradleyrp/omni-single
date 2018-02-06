#!/usr/bin/env python

"""
Plot probabilities of adjacent lipids on the mesh.
"""

import scipy
import scipy.stats
import scipy.interpolate
import itertools,time
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import brewer2mpl
from codes.looptools import basic_compute_loop

### STANDARD

def get_upstream_mesh_partners(abstractor_key):
	"""Filter the upstream sweep for presence/absence of cholesterol."""
	# historical note: previosuly we had to check the calc specs for the right abstractor_key here
	# ... because the data structure was prepared by a function called collect_upstream_calculations_over_loop
	# ... however the new omnicalc lets you loop over calculations automatically. it's your job to unpack them
	try: candidates = [key for key in data_collect.keys() if 
		#! literally the longest nested dictionary call I have ever written I'm sorry -RPB
		calc_collect[key]['calcs']['specs']['upstream']['lipid_mesh'][
		'upstream']['lipid_abstractor']['selector']==abstractor_key]
	except: candidates = []
	if len(candidates)!=1: 
		raise Exception('failed to get a key for abstractor "%s" from data from plotload'%abstractor_key)
	else: 
		key = candidates[0]
		return dict([(sn,data_collect[key][sn]['data']) for sn in data_collect[key].keys()])

def prep_postdat(nn,sns,abstractor='lipid_com',monolayer_indexer=None):
	data = get_upstream_mesh_partners(abstractor)
	postdat = dict()
	for sn in sns:
		if monolayer_indexer:
			top_mono = monolayer_indexer(sn=sn,abstractor=abstractor)
		else: top_mono = work.meta[sn].get('index_top_monolayer',0)
		combonames = data[sn]['combonames_%d'%nn]
		#---! deprecated method before using more basic method
		if False:
			counts_random = data[sn]['counts_random_%d'%nn].mean(axis=0)[top_mono].mean(axis=0)
			counts_obs = data[sn]['counts_observed_%d'%nn][top_mono].mean(axis=0)
			counts_random_err = data[sn]['counts_random_%d'%nn].std(axis=0)[top_mono].mean(axis=0)
			inds = np.where(counts_random>0)[0]
			ratios = counts_obs[inds]/counts_random[inds]
			ratios_err = np.array([data[sn]['counts_random_%d'%nn][i][top_mono].mean(axis=0) 
				for i in range(len(data[sn]['counts_random_%d'%nn]))]).std(axis=0)
			rename_ptdins = lambda x: 'PtdIns' if x==work.meta[sn]['ptdins_resname'] else x
			combonames_std = np.array([[rename_ptdins(i) for i in j] for j in combonames[inds]])
			postdat[sn] = dict(ratios=ratios,ratios_err=ratios_err,combos=combonames_std)
		combi = scipy.misc.comb
		combos = data[sn]['combonames_%d'%nn]
		reslist = data[sn]['reslist']
		comps = dict([(reslist[int(i)],j) 
			for i,j in zip(*np.unique(data[sn]['monolayer_residues_%d'%top_mono],return_counts=True))])
		inds = np.where([np.in1d(c,comps.keys()).all() for c in combos])[0]
		ratios_random = np.array([np.product([combi(comps[k],v) 
			for k,v in zip(*np.unique(c,return_counts=True))])/combi(sum(comps.values()),nn) 
			for c in combos[inds]])
		counts_obs = data[sn]['counts_observed_%d'%nn][top_mono].mean(axis=0)
		counts_obs_std = data[sn]['counts_observed_%d'%nn][top_mono].std(axis=0)
		nmols = data[sn]['monolayer_residues_%d'%top_mono].shape[0]
		ratios = counts_obs[inds]/(counts_obs[inds].sum())/ratios_random
		ratios_err = counts_obs_std[inds]/(counts_obs[inds].sum())/ratios_random
		rename_ptdins = lambda x: 'PtdIns' if x==work.meta[sn]['ptdins_resname'] else x
		combonames_std = np.array([[rename_ptdins(i) for i in j] for j in combonames[inds]])
		postdat[sn] = dict(ratios=ratios,ratios_err=ratios_err,combos=combonames_std)
	if not all([len(postdat[sns[0]]['combos'])==len(postdat[sn]['combos']) 
		and np.all(postdat[sns[0]]['combos']==postdat[sn]['combos']) for sn in sns]):
		for sn in sns:
			status('simulation %s has combinations:\n%s'%(sn,postdat[sn]['combos']),tag='debug')
		sns_wrong = [sn for sn in sns if not np.all(postdat[sns[0]]['combos']==postdat[sn]['combos'])]
		status(tag='debug',string='the following simulations have uneven combination names '
			'(compared to %s): %s'%(sns[0],sns_wrong))
		raise Exception('uneven combination names. see debugging above. check your monolayer indices')
	postdat['combos'] = postdat[sns[0]]['combos']
	return postdat

@autoplot(plotrun)
def plot_partners_ptdins():
	"""Create several plots of the partner data."""
	specs = {
		'summary':{
			'sns':work.metadata.collections['position'],
			'extras':{'special':True,'summary':True,'legend_below':True},
			'specs':{
				0:dict(nn=3,abstractor='lipid_com',combos=[['Ptdins','Ptdins','Ptdins']]),
				1:dict(nn=2,abstractor='lipid_com'),},
			'panelspec':dict(figsize=(12,8),
				layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,2],
				'wratios':[1,6],'wspace':0.1},]}),},
		'comprehensive_core':{
			'sns':work.metadata.collections['position'],
			'extras':{'small_labels_ax3':True,'all_y_labels':True,'legend_everywhere':True},
			'specs':{
				0:dict(nn=2,abstractor='lipid_com'),
				1:dict(nn=3,abstractor='lipid_com'),
				2:dict(nn=2,abstractor='lipid_chol_com'),
				3:dict(nn=3,abstractor='lipid_chol_com'),},
			'panelspec':dict(figsize=(18,30),
				layout={'out':{'grid':[1,1]},
				'ins':[{'grid':[4,1],'wspace':0.1} for i in range(1)]}),},
		'comprehensive_core_wide':{
			'sns':work.metadata.collections['position'],
			'extras':{'four_plots':True},
			'specs':{
				0:dict(nn=2,abstractor='lipid_com'),
				2:dict(nn=3,abstractor='lipid_com'),
				1:dict(nn=2,abstractor='lipid_chol_com'),
				3:dict(nn=3,abstractor='lipid_chol_com'),},
			'panelspec':dict(figsize=(36,18),
				layout={'out':{'grid':[1,2],'wratios':[1,2],'wspace':0.1},
				'ins':[{'grid':[2,1]} for i in range(2)]}),},
		'comprehensive':{
			#---! what's wrong with v536 has only POPC-POPC
			'sns':[i for i in work.metadata.collections['asymmetric_all'] if i!='membrane-v536'],
			'extras':{'small_labels_ax3':True,'error_bars':False,
				'all_y_labels':True,'legend_everywhere':True},
			'specs':{
				0:dict(nn=2,abstractor='lipid_com'),
				2:dict(nn=3,abstractor='lipid_com'),
				1:dict(nn=2,abstractor='lipid_chol_com'),
				3:dict(nn=3,abstractor='lipid_chol_com'),},
			'panelspec':dict(figsize=(18,30),
				layout={'out':{'grid':[1,1]},
				'ins':[{'grid':[4,1],'wspace':0.1} for i in range(1)]}),},}
	for figname,spec in specs.items(): 
		plot_partners_basic(figname=figname,**spec)

@autoplot(plotrun)
def plot_partners_actinlink():
	"""Create several plots of the partner data."""
	sns_mdia2_ordering = ['mdia2bilayer10','mdia2bilayer10_2',
		'mdia2bilayer_nochl2','mdia2bilayer_nochl3',
		'mdia2bilayerphys','mdia2bilayerphys2','mdia2bilayer30','mdia2bilayer30_2']
	sns_mdia2_ordering_chol = ['mdia2bilayer10','mdia2bilayer10_2',
		'mdia2bilayerphys','mdia2bilayerphys2','mdia2bilayer30','mdia2bilayer30_2']
	global sns
	sns = sns_mdia2_ordering
	#! replicate mapping from actinlink_bonds_analysis
	global color_by_simulation,replicate_mapping,extra_labels
	extra_labels = {
		'pip2_20_no_chol':r'mDia2, 20% $PIP_2$, no CHOL ($\times2$)',
		'pip2_20':r'mDia2, 20% $PIP_2$ ($\times2$)',
		'pip2_30':r'mDia2, 30% $PIP_2$ ($\times2$)',
		'pip2_10':r'mDia2, 10% $PIP_2$ ($\times2$)',}
	replicate_mapping = [
		('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
		('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
		('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
		('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])]
	#! actinlink settings go in the art file please !!! this is verbatim from actinlink_bonds_analysis
	def color_by_simulation(sn):
		# tetradic colors via https://www.sessions.edu/color-calculator/
		colors = ['#ff6c28','#28c6ff','#a928ff','#ffd828']
		# tetradic colors via http://paletton.com/#uid=7030Z0kqbujggFvlnx6qcY+wDkl
		colors = ['#f2b52c','#22b93c','#2e47a4','#f23a2c']
		refs = ['^mdia2bilayer_nochl','^mdia2bilayer(?!(_nochl|[0-9]))','^mdia2bilayer10','^mdia2bilayer30']
		ref_to_color = dict(zip(refs,colors))
		matches = [color for ref,color in ref_to_color.items() if re.match(ref,sn)]
		if len(matches)>1: raise Exception
		elif len(matches)==0:
			#! handling combos here with a minor hack see global replicate_mapping
			return dict(zip(zip(*replicate_mapping)[0],colors))[sn]
		else: return matches[0]
	#! note that monolayer indexer is different for cholesterol vs no cholesterol so we are hardcoding the
	#! ... right indices here instead of the metadata
	def monolayer_indexer(sn,abstractor):
		if abstractor=='lipid_chol_com':
			if sn in ['mdia2bilayer10_2','mdia2bilayer30_2']: return 1
			else: return 0
		elif abstractor=='lipid_com':
			if sn in ['mdia2bilayer10','mdia2bilayerphys2','mdia2bilayer30_2']: return 1
			else: return 0
		else: raise Exception
	#! the following flags divert the usual flow to plot for actinlink: legend_scheme, hatch, color_scheme
	specs = {
		'mesh_with_chol':{
			'sns':sns_mdia2_ordering_chol,
			'extras':{'special':True,'summary':True,'legend_below':True,
				'color_scheme':'actinlink','legend_scheme':'actinlink_nochol','hatch':False,
				'monolayer_indexer':monolayer_indexer},
			'specs':{
				0:dict(nn=3,abstractor='lipid_chol_com'),
				1:dict(nn=2,abstractor='lipid_chol_com'),},
			'panelspec':dict(figsize=(20,8),
				layout={'out':{'grid':[1,1]},'ins':[{'grid':[2,1],
				'wspace':0.1},]}),},
		'mesh_without_chol':{
			'sns':sns_mdia2_ordering[:],
			'extras':{'special':True,'summary':True,'legend_below':True,
				'color_scheme':'actinlink','legend_scheme':'actinlink','hatch':False,
				'monolayer_indexer':monolayer_indexer},
			'specs':{
				0:dict(nn=3,abstractor='lipid_com'),
				1:dict(nn=2,abstractor='lipid_com'),},
			'panelspec':dict(figsize=(20,8),
				layout={'out':{'grid':[1,1]},'ins':[{'grid':[2,1],
				'wspace':0.1},]}),},}
	for figname,spec in specs.items(): 
		status('plotting %s'%figname,tag='plot')
		plot_partners_basic(figname=figname,**spec)

def plot_partners_basic(sns,figname,specs,panelspec,extras):
	"""Summarize the lipid mesh partners."""
	baseline = 1.0
	sns_reorder = sns
	plotspec = ptdins_manuscript_settings()
	axes,fig = panelplot(**panelspec)
	try: axes = [i for j in axes for i in j]
	except: pass
	for axnum,details in specs.items():
		nn = details['nn']
		kwargs_to_postdat = {}
		if extras.get('monolayer_indexer',False): 
			kwargs_to_postdat['monolayer_indexer'] = extras['monolayer_indexer']
		if extras.get('postdat',False): 
			postdat = dict([(sn,extras['postdat'][(sn,nn)]) for sn in sns])
			postdat['combos'] = extras['postdat']['combos_%d'%nn]
		# the original script computes the postdata
		else: postdat = prep_postdat(nn=nn,sns=sns,abstractor=details['abstractor'],**kwargs_to_postdat)
		combos = details.get('combos',postdat['combos'])
		ax = axes[axnum]
		max_y,min_y = 0,10
		combo_spacer,half_width = 0.5,0.5
		for snum,sn in enumerate(sns_reorder):
			# already ensured each simulation has the same combinations
			for cnum,combo in enumerate(combos):
				xpos = (cnum)*len(sns)+snum+cnum*combo_spacer
				ypos = postdat[sn]['ratios'][cnum]-baseline
				color_scheme = extras.get('color_scheme',False)
				if color_scheme=='actinlink':
					#! temporary hack
					color = color_by_simulation(sn)
				else: color = plotspec['colors'][plotspec['colors_ions'][work.meta[sn]['cation']]]
				hatch = extras.get('hatch',plotspec['hatches_lipids'][work.meta[sn]['ptdins_resname']])
				ax.bar([xpos],[ypos],
					width=2*half_width,bottom=baseline,
					color=color,
					hatch=hatch,**({'label':work.meta[sn]['ion_label']} if cnum==len(combos)-1 else {}))
				y_err = postdat[sn]['ratios_err'][cnum]
				if extras.get('error_bars',True): ax.errorbar(xpos,ypos+baseline,yerr=y_err,
					alpha=1.0,lw=2.0,c='k')
				max_y,min_y = max((max_y,ypos+baseline+y_err)),min((min_y,ypos+baseline-y_err))
		ax.set_xlim(-half_width-half_width/2.,(len(sns))*len(combos)-
			half_width+(len(combos)-1)*combo_spacer+half_width/2.)
		ax.axhline(1.0,c='k',lw=2)
		dividers = np.array([len(sns)*(cnum+1)+cnum*half_width-half_width/2.
			for cnum in range(-1,len(combos))])
		for divider in dividers: ax.axvline(divider,c='k',lw=1)
		mids = (dividers[1:]+dividers[:-1])/2.
		ax.set_xticks(mids)
		ax.set_xticklabels(['\n'.join(c) for c in combos],fontsize=18)
		ax.xaxis.tick_top()
		labelsize = 14
		if extras.get('small_labels_ax3',False) and axnum==3: labelsize = 10
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on',labelsize=14)
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labeltop='on',labelsize=labelsize)
		#! ax.set_ylim((min_y*0.9,max_y*1.1))
		if (axnum==len(axes)-1 or extras.get('legend_everywhere',False) 
			and extras.get('legend_scheme',False) in ['actinlink_nochol','actinlink']):
			#! custom actinlink legends
			legendspec = []
			if extras['legend_scheme']=='actinlink_nochol': 
				match_list = [replicate_mapping[i] for i in [0,2,3]]
			else: match_list = replicate_mapping
			for sn,matches in match_list:
				legendspec.append(dict(name=extra_labels[sn],
					patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(matches[0]))))
			patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
			legend = ax.legend(patches,labels,loc='lower left',
				bbox_to_anchor=(1.0,0.0,1.,1.),ncol=1,fontsize=14)
			frame = legend.get_frame()
			frame.set_edgecolor('k')
			frame.set_facecolor('white')
		elif axnum==len(axes)-1 or extras.get('legend_everywhere',False):
			bar_formats = make_bar_formats(sns_reorder,work=work)
			comparison_spec = dict(ptdins_names=list(set([work.meta[sn]['ptdins_resname'] for sn in sns])),
				ion_names=list(set([work.meta[sn]['cation'] for sn in sns])))
			kwargs = dict(bbox=(1.05,0.0,1.,1.),loc='upper left')
			if extras.get('legend_below',False):
				kwargs = dict(bbox=(0.5,-0.05),loc='upper center',ncol=4)
			legend,patches = legend_maker_stylized(ax,work=work,
				sns_this=sns_reorder,bar_formats=bar_formats,comparison_spec=comparison_spec,
				**kwargs)
	# make the unity line even between subplots (used algrebra in real life here)
	# ... decided to shift down the right axis y plot since it has a smaller range, to make things even
	#! note that this is not perfect but we should just pick our battles
	if extras.get('special',False):
		spread_prop = np.array([[i for i in ax.get_ylim()] for ax in axes[:2]])
		s0,s1,s2,s3 = spread_prop.reshape(-1)
		shift_down = ((1.-s2)-(1.-s0)*(s3-1.)/(s1-1.))
		axes[1].set_ylim((axes[1].get_ylim()[0],axes[1].get_ylim()[1]+shift_down))
	axes[0].set_ylabel('observations relative to chance',fontsize=16)
	if extras.get('four_plots',False):
		axes[1].set_ylabel('observations relative to chance',fontsize=16)
	if extras.get('all_y_labels',False):
		for ax in axes: ax.set_ylabel('observations relative to chance',fontsize=16)
	picturesave('fig.lipid_mesh_partners.%s'%figname,work.plotdir,
		backup=False,version=True,meta={},extras=[legend])#,form='pdf')

@autoload(plotrun)
def load():
	"""Load everything for the plot only once."""
	data_collect,calc_collect = work.plotload('lipid_mesh_partners')
	sns = work.sns()

# one plot script for multiple projects
plotrun.routine = ['plot_partners_%s'%os.path.basename(os.getcwd())]
