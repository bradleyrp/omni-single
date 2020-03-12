#!/usr/bin/env python

"""
incoming from script-partners.py
"""

import scipy
import scipy.stats
import scipy.interpolate
import itertools,time
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---settings
#---we must declare variables here. this should match the top of the reload function
variables = 'sns,data,calc,postdat'.split(',')

from codes.looptools import basic_compute_loop

###---STANDARD

#---declare standard variables
required_variables = 'printers,routine'.split(',')
for v in variables+required_variables:
	if v not in globals(): globals()[v] = None

def register_printer(func):
	"""Add decorated functions to a list of "printers" which are the default routine."""
	global printers
	if printers is None: printers = []
	printers.append(func.__name__)
	return func

###---PLOTS

def partner_counter(sn,fr,mn):
	"""..."""
	#---globals from parent function
	global data,reslist,rxmlook,combolookup_str,combolookup3_str,tally,tally_random,tally3,tally3_random
	dat = data[sn]['data']
	i2s = lambda *k: '.'.join([str(j) for j in k])

	simplices = dat[i2s(mn,fr,'simplices')]
	gids = dat[i2s(mn,fr,'ghost_ids')]
	links = np.array(list(set([tuple(np.sort(j)) for j in 
		np.concatenate([np.transpose((gids[simplices][:,i],gids[simplices][:,(i+1)%3])) 
		for i in range(3)])])))
	#---uniquely identify each link type
	links_type = np.transpose([np.sum(rxmlook[mn][links]==i,axis=1) for i in range(len(reslist))])
	links_type_str = [''.join(['%s'%s for s in i]) for i in links_type]
	counts = dict(scipy.stats.itemfreq(links_type_str))
	tally[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

	#---sidestepping the probability questions and randomizing the identifier
	#---...which will effectively randomize the identities of the vertices
	#---...and then later we will take the average of these and then use it to see if observed
	#---...links are more or less likely to appear in our graph than in a random one
	rxmlook_rand = [np.random.permutation(r) for r in rxmlook]
	links_type_wrong = np.transpose([np.sum(rxmlook_rand[mn][links]==i,axis=1) 
		for i in range(len(reslist))])
	links_type_str = [''.join(['%s'%s for s in i]) for i in links_type_wrong]
	counts = dict(scipy.stats.itemfreq(links_type_str))
	tally_random[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

	#---identify each triple type
	triples = gids[simplices]
	triple_type = np.transpose([np.sum(rxmlook[mn][triples]==i,axis=1) for i in range(len(reslist))])
	triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type]
	counts = dict(scipy.stats.itemfreq(triple_type_str))
	tally3[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup3_str])

	#---randomize the triples
	triple_type_wrong = np.transpose([np.sum(rxmlook_rand[mn][triples]==i,axis=1) 
		for i in range(len(reslist))])
	triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type_wrong]
	counts = dict(scipy.stats.itemfreq(triple_type_str))
	tally3_random[mn][fr] = np.array([int(counts[i]) 
		if i in counts else 0 for i in combolookup3_str])

def partner_watcher(sn,**kwargs):

	#dat = dats[0]
	global data
	dat = data[sn]['data']
	i2s = lambda *k: '.'.join([str(j) for j in k])

	#---from binding_combinator
	resnames = np.array(dat['resnames'])
	nframes = dat['nframes']
	lipids = np.array(list(resnames[np.sort(np.unique(resnames,return_index=True)[1])]))

	#---replicate some of the binding_combinator and compute_bridging_norms codes for pairs
	nn = 2-1
	combos = np.array([''.join(j) for j in 
		itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) 
		if sum([int(k) for k in j])==nn+1])
	combonames = [tuple(v) for v in [np.concatenate([[lipids[ww]]*int(w) 
		for ww,w in enumerate(l)]) for l in combos]]
	nn = 3-1
	combos3 = np.array([''.join(j) for j in 
		itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) 
		if sum([int(k) for k in j])==nn+1])
	combonames3 = [tuple(v) for v in [np.concatenate([[lipids[ww]]*int(w) 
		for ww,w in enumerate(l)]) for l in combos3]]

	#---globals for parallel
	global reslist,rxmlook,combolookup_str,combolookup3_str,tally,tally_random,tally3,tally3_random

	#---determine monolayer-specific residue indices
	imono = dat['monolayer_indices']
	nmols = [np.sum(dat['monolayer_indices']==i) for i in range(2)]
	resnames = np.array(dat['resnames'])
	reslist = list(np.array(resnames)[np.sort(np.unique(resnames,return_index=True)[1])])
	rxm = [[
		np.array([np.where(np.where(imono==mn)[0]==i)[0][0] 
		for i in np.where(np.all((imono==mn,resnames==rn),axis=0))[0]])
		for rn in reslist] for mn in range(2)]
	rxmlook = [np.zeros(n) for n in nmols]
	for mn in range(2):
		for ri,r in enumerate(rxm[mn]):
			if r != []: rxmlook[mn][r] = ri
		
	#---record the link count trajectories
	combolookup = np.sum([np.array(combonames)==r for r in reslist],axis=2).T
	combolookup_str = [''.join(['%s'%s for s in i]) for i in combolookup]
	combolookup3 = np.sum([np.array(combonames3)==r for r in reslist],axis=2).T
	combolookup3_str = [''.join(['%s'%s for s in i]) for i in combolookup3]
	tally = [np.zeros((nframes,len(combolookup))) for mn in range(2)]
	tally_random = [np.zeros((nframes,len(combolookup))) for mn in range(2)]
	tally3 = [np.zeros((nframes,len(combolookup3))) for mn in range(2)]
	tally3_random = [np.zeros((nframes,len(combolookup3))) for mn in range(2)]

	#looper = [dict(sn=sn,mn=mn,fr=fr) for mn in range(2) for fr in range(nframes)]
	#incoming = basic_compute_loop(partner_counter,looper)

	for mn in range(2):
		for fr in range(nframes):
			status('%s monolayer %d'%(sn,mn),i=fr,looplen=nframes,tag='compute')
			simplices = dat[i2s(mn,fr,'simplices')]
			gids = dat[i2s(mn,fr,'ghost_ids')]
			links = np.array(list(set([tuple(np.sort(j)) for j in 
				np.concatenate([np.transpose((gids[simplices][:,i],gids[simplices][:,(i+1)%3])) 
				for i in range(3)])])))
			#---uniquely identify each link type
			links_type = np.transpose([np.sum(rxmlook[mn][links]==i,axis=1) for i in range(len(reslist))])
			links_type_str = [''.join(['%s'%s for s in i]) for i in links_type]
			counts = dict(scipy.stats.itemfreq(links_type_str))
			tally[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

			#---sidestepping the probability questions and randomizing the identifier
			#---...which will effectively randomize the identities of the vertices
			#---...and then later we will take the average of these and then use it to see if observed
			#---...links are more or less likely to appear in our graph than in a random one
			rxmlook_rand = [np.random.permutation(r) for r in rxmlook]
			links_type_wrong = np.transpose([np.sum(rxmlook_rand[mn][links]==i,axis=1) 
				for i in range(len(reslist))])
			links_type_str = [''.join(['%s'%s for s in i]) for i in links_type_wrong]
			counts = dict(scipy.stats.itemfreq(links_type_str))
			tally_random[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

			#---identify each triple type
			triples = gids[simplices]
			triple_type = np.transpose([np.sum(rxmlook[mn][triples]==i,axis=1) for i in range(len(reslist))])
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type]
			counts = dict(scipy.stats.itemfreq(triple_type_str))
			tally3[mn][fr] = np.array([int(counts[i]) if i in counts else 0 for i in combolookup3_str])

			#---randomize the triples
			triple_type_wrong = np.transpose([np.sum(rxmlook_rand[mn][triples]==i,axis=1) 
				for i in range(len(reslist))])
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type_wrong]
			counts = dict(scipy.stats.itemfreq(triple_type_str))
			tally3_random[mn][fr] = np.array([int(counts[i]) 
				if i in counts else 0 for i in combolookup3_str])

	results = {}
	results ['nframes'] = np.array(nframes)
	results ['resnames'] = np.array(resnames).astype(str)
	for mn in range(2): results['pairs%d'%mn] = np.array(tally[mn])
	for mn in range(2): results['pairs_random%d'%mn] = np.array(tally_random[mn])
	for mn in range(2): results['triples%d'%mn] = np.array(tally3[mn])
	for mn in range(2): results['triples_random%d'%mn] = np.array(tally3_random[mn])
	results ['combonames'] = np.array(combonames).astype(str)
	results ['combonames3'] = np.array(combonames3).astype(str)
	return results

def reload():
	"""Load everything for the plot only once."""
	#---canonical globals list
	#---!? can this be made programmatic?
	global sns,data,calc,postdat
	#---reload sequence goes here
	data,calc = plotload(plotname)
	sns = work.sns()
	postdat = dict([(sn,partner_watcher(sn)) for sn in sns])

###---STANDARD

def printer():
	"""Load once per plot session."""
	global variables,routine,printers
	#---reload if not all of the globals in the variables
	if any([v not in globals() or globals()[v] is None for v in variables]): reload()
	#---after loading we run the printers
	printers = list(set(printers if printers else []))
	if routine is None: routine = list(printers)	
	#---routine items are function names
	for key in routine: 
		status('running routine %s'%key,tag='printer')
		globals()[key]()

if __name__=='__main__': printer()

if True:

	#---! monolayer top! one loop! etc

	import brewer2mpl
	colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
		for key,val in {
		'red':0,'blue':1,'green':2,'purple':3,'orange':4,
		'yellow':5,'brown':6,'pink':7,'grey':8,}.items()])
	colors['white'] = mpl.colors.ColorConverter().to_rgb("#000000")
	colors['pink'] = mpl.colors.ColorConverter().to_rgb("#f1948a")
	colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")
	bgw = 0.3
	colors['bluegreen'] = np.average((colors['blue'],colors['green']),weights=[1-bgw,bgw],axis=0)
	bgw2 = 0.6
	colors['bluegreen2'] = np.average((colors['blue'],colors['green']),weights=[1-bgw2,bgw2],axis=0)
	scale = 0.5
	colors['light_blue'] = np.average((colors['blue'],colors['white']),weights=[1-scale,scale],axis=0)

	colors_ions = {'NA':'green','Na,Cal':'bluegreen','MG':'pink','Cal':'blue','K':'grey',}
	hatches_lipids = {'PI2P':'//','P35P':'-','PIPU':'xx','PIPP':'++','SAPI':''}

	ratios = dict()
	for style in ['pairs','triples']:
		vals = [[] for sn in sns]
		for snum,sn in enumerate(sns):
			mn = 0
			post = postdat[sn]
			combonames = post['combonames%s'%{'triples':'3','pairs':''}[style]]
			for cnum,combo in enumerate(combonames):
				ratio = post['%s0'%style].mean(axis=0)[cnum]/post['%s_random0'%style].mean(axis=0)[cnum]
				ratio_errs = post['%s0'%style].mean(axis=0)[cnum]/post['%s_random0'%style].mean(axis=0)[cnum]
				if not np.isnan(ratio): vals[snum].append(ratio)
				else: vals[snum].append(0.0)
		ratios[style] = vals

	fig = plt.figure(figsize=(14,12))
	for style_num,style in enumerate(['pairs','triples']):
		ax = plt.subplot(212 if style_num==1 else 211)
		counter = 0
		ticks = []
		vals = ratios[style]
		combonames = post['combonames%s'%{'triples':'3','pairs':''}[style]]
		for cnum,combo in enumerate(combonames):
			vals_sub = np.array([vals[snum][cnum] for snum,sn in enumerate(sns)])
			if (vals_sub==0).all(): continue
			xs = np.arange(len(vals_sub))+counter
			for snum,(x,y) in enumerate(zip(xs,vals_sub-1.0)):
				ax.bar([x],[y],width=1.0,bottom=1.0,
					hatch=hatches_lipids[work.meta[sns[snum]]['ptdins_resname']],
					color=colors[colors_ions[work.meta[sns[snum]]['cation']]])
			counter += len(vals_sub)+0.5
			ax.axvline(counter-1.0+0.25,c='k',lw=1)
			ticks.append((combo,xs.mean()))
		ax.set_xticks(zip(*ticks)[1])
		ax.set_xticklabels(['\n'.join(i) for i in zip(*ticks)[0]])
		ax.axhline(1.0,c='k',lw=2)
		ax.set_ylabel('links compared to random')
	picturesave('fig.TEST',work.plotdir,backup=False,version=True,meta={},extras=[])

#---! FAILURE
if False:

	print(66)
	figsize = (12,12)
	bardat,barspec,extras = [],[],[]
	bardat = [[[(('a','b'),(i,1.0)) for i in [3,4,6,7]],[(('a','b'),(i,1.0)) for i in [3,4,6,7]]]]
	barspec = [[[{'color':'k'} for i in [3,4,6,7]],[{'color':'k'} for i in [3,4,6,7]]]]
	bar_groups = dict(bardat=bardat,barspec=barspec,extras=extras)
	#axes,fig = panelplot(figsize=figsize,layout={'out':{'grid':[1,1],},
	#	'ins':[{'grid':[1,1],} for i in range(1)]})
	#ax = axes[0]
	fig = plt.figure()
	ax = plt.subplot(111)
	bars = BarGrouper(ax=ax,fig=fig,figsize=figsize,show_xticks=False,**bar_groups)
	picturesave('fig.TEST2',work.plotdir,backup=False,version=True,meta={},extras=[])

