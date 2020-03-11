#!/usr/bin/env python

"""
Plot lipid-lipid g(r) results.
"""

import scipy
import scipy.integrate
import re

class Post:
	"""Data packaged for plotting over upstream calculation loop."""
	def __init__(self): self.dataset,self.calcset = {},{}
	def new(self,name,**kwargs):
		if name in self.dataset: 
			raise Exception('key already exists: %s'%name)
		self.calcset[name] = kwargs.pop('calc',{})
		self.dataset[name] = kwargs
		self.select(name)
	def select(self,name):
		self.data = self.dataset[name]
		self.calc = self.calcset[name]
		self.cursor = name

class Same:
	"""Make sure items are the same between replicates."""
	def __init__(self): self.this = []
	def add(self,item): self.this.append(item)
	def get(self):
		if not all([self.this[0]==i for i in self.this[1:]]): raise Exception('not same')
		else: return self.this[0]

def make_legend(ax,fig,**kwargs):
	# legend at last tile
	legendspec = []
	keys = kwargs.pop('keys',None)
	legend_below = kwargs.pop('legend_below',False)
	legend_ncol = kwargs.pop('legend_ncol',False)
	if kwargs: raise Exception
	for sn_general,sns in replicate_mapping:
		if not keys or sn_general in keys:
			legendspec.append(dict(name=extra_labels[sn_general],
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	if legend_below:
		legend = fig.legend(patches,labels,
			loc='upper center',ncol=legend_ncol if legend_ncol else len(legendspec),fontsize=14)
	else: legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor('w')
	frame.set_facecolor('w')
	return legend

def attach_selector_key(meta):
	"""Include the abstractor_selector in the metadata."""
	try: 
		abstractor_selector = work.plotspec.specs['calculation'][
			'lipid_rdfs']['upstream']['lipid_abstractor']['selector']
		meta['abstractor_selector'] = abstractor_selector
	except: pass

if False:
	def reducer(pairings,pair_residues):
		"""Final reduction in the data."""
		post = dict([])
		for pnum,(pairname,pair) in enumerate(zip(pairings,pair_residues)):
			post[pairname] = {}
			for sn_group,sns_this in replicate_mapping:
				obs,nmols,total_areas = [],Same(),[]
				for sn in sns_this:
					mn = actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')
					dat = data[sn]['data']
					key = 'counts_mn%d_%s:%s'%(mn,pairname[0],pairname[1])
					if key not in dat: 
						status('reducing %s, %s SKIP'%(sn,key))
						continue
					else: status('reducing %s, %s'%(sn,key))
					obs.append(dat[key])
					nmols.add([np.in1d(dat['resnames'][dat['monolayer_indices']==mn],p).sum() 
						for p in pair])
					total_areas.append(dat['total_area'])
				if len(obs) not in [0,2]: raise Exception('something has gone wrong')
				elif len(obs)==0: continue
				# average replicates
				counts = np.array(obs).mean(axis=0)
				nmols = nmols.get()
				total_area = np.mean(total_areas) #! total area is averaged between simulations for norming
				density = nmols[0]*nmols[1]/total_area
				post[pairname][sn_group] = dict(counts=counts,total_area=total_area,density=density)

@autoload(plotrun)
def load():
	"""Load single or looped upstream RDF data."""
	data,calc = plotload('lipid_rdfs')
	sns = work.sns()
	dpi,img_form = 600,'pdf'
	# handle loops with a Post object
	post = Post()
	# if the calculation is in the top level of calc then we have no loop
	if 'lipid_rdfs' in calc: 
		post.new('lipid_rdfs',data=data['lipid_rdfs'],calc=calc['lipid_rdfs'])
	# upstream loop in the data requires extra selection
	else:
		keys = [key for key,val in calc.items() if key!='extras']
		# name the data in the loop with a path in the specs
		path = 'calcs','specs','upstream','lipid_abstractor','selector'
		for key in keys:
			name = delve(calc[key],*path)
			post.new(name,data=data[key],calc=calc[key])
	# generics dumped to globals
	#! port to ptdins color_by_simulation = actinlink_color_by_simulation
	#! all of this is coming from ptdins
	color_by_simulation = lambda sn: ptdins_color_by_simulation(meta=work.meta[sn])
	xmax_cumulative = 2.0
	# loop over upstream calculations
	for upstream_name in post.dataset:
		post.select(upstream_name)
		cutoff,binsize = [post.calc['calcs']['specs'][i] for i in ['cutoff','binsize']]
		scanrange = np.arange(0,cutoff,binsize)
		middles = (scanrange[1:]+scanrange[:-1])/2
		areas = np.array([np.pi*(binsize)*middles[i]*2 for i in range(len(middles))])
		data_this = post.data['data']
		# collect all possible pairings
		pairings = list(set([tuple(i) for j in [
			data_this[sn]['data']['pairnames_mn%d'%
			actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')] 
			for sn in sns] for i in j]))
		pairname_to_pair = dict([(tuple(i),j) for i,j in zip(*[data_this[sn]['data']['%s_mn%d'%(p,
			actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com'))] for p in ['pairnames','pairs']])])
		pair_residues = [pairname_to_pair[p] for p in pairings]
		#! porting to ptdins replicate_mapping = actinlink_replicate_mapping
		replicate_mapping = [(i,[i,]) for i in sns]
		#! extra_labels = actinlink_extra_labels
		extra_labels = ptdins_extra_labels
		#! it would be nice to have this calculation in a separate function but then it will not have globals.
		#! ... this issue needs documenting along with the load magic
		post_this = post.data['reduced'] = dict([])
		for pnum,(pairname,pair) in enumerate(zip(pairings,pair_residues)):
			post_this[pairname] = {}
			for sn_group,sns_this in replicate_mapping:
				obs,nmols,total_areas = [],Same(),[]
				for sn in sns_this:
					mn = actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')
					dat = data_this[sn]['data']
					key = 'counts_mn%d_%s:%s'%(mn,pairname[0],pairname[1])
					if key not in dat: 
						status('reducing %s, %s SKIP'%(sn,key))
						continue
					else: status('reducing %s, %s'%(sn,key))
					obs.append(dat[key])
					nmols.add([np.in1d(dat['resnames'][dat['monolayer_indices']==mn],p).sum() 
						for p in pair])
					total_areas.append(dat['total_area'])
				if len(obs)==1:
					print('[WARNING] if this is supposed to be actinlink we only have one replicate!')
				elif len(obs) not in [0,2]: raise Exception('something has gone wrong')
				# porting to ptdins
				elif len(obs)==0: continue
				# average replicates
				counts = np.array(obs).mean(axis=0)
				nmols = nmols.get()
				total_area = np.mean(total_areas) #! total area is averaged between simulations for norming
				density = nmols[0]*nmols[1]/total_area
				post_this[pairname][sn_group] = dict(counts=counts,
					total_area=total_area,density=density,nmols=nmols)
		# load necessary globals into the post object for later
		#! note that quantities which vary with upstream data have to be stored in post otherwise
		#! ... they go directly to globals and may linger there during later loops over upstream data
		post.data['pairings'] = pairings
		post.data['replicate_mapping'] = replicate_mapping
		post.data['pair_residues'] = pair_residues
		post.data['middles'] = middles
		post.data['areas'] = areas

def plot_radial_distributions(figname,mode='summary',**kwargs):
	"""Plot RDFs or CDFs of the RDFs."""
	# plot lipid-lipid RDFs on one panel
	reps_subset = kwargs.get('groups',None)
	pairings_this = kwargs.get('pairings',pairings)
	figsize = kwargs.pop('figsize',(12,12))
	legend_spec = kwargs.get('legend_spec',{})
	if reps_subset: replicate_mapping_this = [(i,j) for i,j in replicate_mapping if i in reps_subset]
	else: replicate_mapping_this = replicate_mapping
	pair_residues_this = [pair_residues[ii] for ii,i in enumerate(pairings) if i in pairings_this]
	axes,fig = square_tiles(len(pairings_this),figsize=figsize,hspace=0.4,wspace=0.4)
	for pnum,(pairname,pair) in enumerate(zip(pairings_this,pair_residues_this)):
		for sn_group,sns_this in replicate_mapping_this:
			if sn_group not in reduced[pairname]: continue
			counts,nmols,total_area,density = [reduced[pairname][sn_group][k] 
				for k in ['counts','nmols','total_area','density']]
			ax = axes[pnum]
			#! show 90% of the maximum range
			xmax = 0.9*np.sqrt(total_area)
			valid = middles<xmax
			if mode=='cdf': 
				valid = middles<xmax_cumulative
				vals = np.array([scipy.integrate.simps(counts[:i]/areas[:i]/density) 
					for i in range(2,len(counts))])[valid]
				ax.plot(middles[valid],vals[valid],color=color_by_simulation(sns_this[0]),lw=2)
				ax.set_xlim((0,xmax_cumulative))
				ax.set_ylabel(r'$\int{g(r)}$')
			elif mode=='summary': 
				valid = middles<xmax
				ax.plot(middles[valid],(counts/areas/density)[valid],
					color=color_by_simulation(sns_this[0]),lw=2)
				ax.axhline(1.0,lw=1,c='k')
				ax.set_xlim((0,xmax))
				ax.set_ylabel('$g(r)$')
				ax.axhline(1.0,lw=0.5,c='k')
			else: raise Exception
			ax.set_xlim((0,xmax_cumulative))
			ax.set_title('%s-%s'%tuple([work.vars['names']['short'].get(p,p) for p in pairname]))
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_xlabel(r'$r\,(nm)$')
	legend = make_legend(axes[-1],fig=fig,keys=[i for i,j in replicate_mapping_this],**legend_spec)
	meta = {}
	attach_selector_key(meta)
	picturesave('fig.lipid_rdfs.%s'%figname,directory=work.plotdir,
		version=True,meta=meta,extras=[legend],dpi=dpi,form=img_form)

@autoplot(plotrun)
def plot_distributions():
	figspec = {
		'summary':{'mode':'summary'},
		'summary.20':{'mode':'summary','groups':['pip2_20_no_chol','pip2_20']},
		'summary.simple':{'mode':'summary','pairings':[
			('all lipids','all lipids'),('PI2P','PI2P'),('DOPS','PI2P'),('DOPS','DOPS'),('CHL1','DOPS')],
			'legend_spec':{'legend_below':False,'legend_ncol':1},'figsize':(6,10)},
		'cumulative':{'mode':'cdf'},}
	figspec = {'summary.simple':figspec['summary.simple']}
	# loop over figures
	for figname,spec in figspec.items(): 
		# loop over upstream data, which is propagated to meta
		for upstream_name in post.dataset:
			post.select(upstream_name)
			# dump the post for each upstream into globals before plotting
			for key in post.data: globals()[key] = post.data[key]
			plot_radial_distributions(figname=figname,**spec)

@autoplot(plotrun)
def plot_cumulative_compare():
	"""Simple bar plot of the cumulative RDF at a cutoff."""
	figname = 'cumulative_compare'
	axes,fig = square_tiles(len(pairings),figsize=(12,12),hspace=0.4,wspace=0.4)
	heights = {}
	for pnum,(pairname,pair) in enumerate(zip(pairings,pair_residues)):
		ax = axes[pnum]
		heights[pairname] = {}
		for snum,(sn_group,sns_this) in enumerate(replicate_mapping):
			if sn_group not in post[pairname]: continue
			counts,nmols,total_area,density = [post[pairname][sn_group][k] 
				for k in ['counts','nmols','total_area','density']]
			# mimic the cdf routine
			valid = middles<xmax_cumulative
			vals = np.array([scipy.integrate.simps(counts[:i]/areas[:i]/density) 
				for i in range(2,len(counts))])[valid]
			# the last value occurs at the xmax_cumulative
			#! check that this is precise?
			heights[pairname][sn_group] = vals[-1]
			ax.bar(snum,vals[-1],color=color_by_simulation(sns_this[0]),width=1.0)
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.set_xticks([])
		ax.set_title('%s-%s'%tuple([work.vars['names']['short'].get(p,p) for p in pairname]))
	ymin = min([min(i.values()) for i in heights.values()])
	ymax = max([max(i.values()) for i in heights.values()])
	for ax in axes:
		ax.set_ylim((0.9*ymin,1.1*ymax))
	legend = make_legend(axes[-1],fig)
	meta = {}
	attach_selector_key(meta)
	picturesave('fig.lipid_rdfs.%s'%figname,directory=work.plotdir,meta=meta,
		version=True,extras=[legend],dpi=dpi,form=img_form)

@autoplot(plotrun)
def plot_ptdins_convergence(): pass
 
plotrun.routine = ['plot_distributions'] # None
plotrun.routine = [] # see development below
if __name__=='__main__': 

	# customized to ptdins. has no replicate mapping
	sns_this = ['membrane-v565']
	this_dataset = 'lipid_chol_com_headless'
	pairname = [[u'DOPE', u'DOPS', u'PI2P'], [u'DOPE', u'DOPS', u'PI2P']]
	for sn in sns_this:
		#!!! is the indexer correct?
		mn = actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')
		dat = data_this[sn]['data']
		key = 'counts_mn%d_%s:%s'%(mn,pairname[0],pairname[1])
		print(dat[key])
