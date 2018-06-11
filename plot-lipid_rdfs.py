#!/usr/bin/env python

"""
Plot lipid-lipid g(r) results.
"""

import scipy
import scipy.integrate
import re

class Same:
	"""Make sure items are the same between replicates."""
	def __init__(self): self.this = []
	def add(self,item): self.this.append(item)
	def get(self):
		if not all([self.this[0]==i for i in self.this[1:]]): raise Exception('not same')
		else: return self.this[0]

def make_legend(ax,keys=None):
	# legend at last tile
	legendspec = []
	for sn_general,sns in replicate_mapping:
		if not keys or sn_general in keys:
			legendspec.append(dict(name=extra_labels[sn_general],
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor('k')
	frame.set_facecolor('w')
	return legend

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
	data,calc = plotload('lipid_rdfs')
	sns = work.sns()
	# generics
	mode = ['summary','cdf'][-1]
	xmax_cumulative = 2.0
	cutoff,binsize = [calc['lipid_rdfs']['calcs']['specs'][i] for i in ['cutoff','binsize']]
	scanrange = np.arange(0,cutoff,binsize)
	middles = (scanrange[1:]+scanrange[:-1])/2
	areas = np.array([np.pi*(binsize)*middles[i]*2 for i in range(len(middles))])
	# collect all possible pairings
	pairings = list(set([tuple(i) for j in [
		data[sn]['data']['pairnames_mn%d'%actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')] 
		for sn in sns] for i in j]))
	pairname_to_pair = dict([(tuple(i),j) for i,j in zip(*[data[sn]['data']['%s_mn%d'%(p,actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com'))] for p in ['pairnames','pairs']])])
	pair_residues = [pairname_to_pair[p] for p in pairings]
	color_by_simulation = actinlink_color_by_simulation
	replicate_mapping = actinlink_replicate_mapping
	extra_labels = actinlink_extra_labels
	#! it would be nice to have this calculation in a separate function but then it will not have globals.
	#! ... this issue needs documenting along with the load magic
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
			post[pairname][sn_group] = dict(counts=counts,total_area=total_area,density=density,nmols=nmols)

def plot_radial_distributions(figname,mode='summary',**kwargs):
	"""Plot RDFs or CDFs of the RDFs."""
	# plot lipid-lipid RDFs on one panel
	reps_subset = kwargs.get('groups',None)
	if reps_subset: replicate_mapping_this = [(i,j) for i,j in replicate_mapping if i in reps_subset]
	else: replicate_mapping_this = replicate_mapping
	axes,fig = square_tiles(len(pairings),figsize=(12,12),hspace=0.4,wspace=0.4)
	for pnum,(pairname,pair) in enumerate(zip(pairings,pair_residues)):
		for sn_group,sns_this in replicate_mapping_this:
			if sn_group not in post[pairname]: continue
			counts,nmols,total_area,density = [post[pairname][sn_group][k] 
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
				ax.set_ylabel('$\int{g(r)}$')
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
			ax.set_xlabel('$r\,(nm)$')
	legend = make_legend(axes[-1],keys=[i for i,j in replicate_mapping_this])
	picturesave('fig.lipid_rdfs.%s'%figname,directory=work.plotdir,meta={},extras=[legend])

@autoplot(plotrun)
def plot_distributions():
	figspec = {
		'summary':{'mode':'summary'},
		'summary.20':{'mode':'summary','groups':['pip2_20_no_chol','pip2_20']},
		'cumulative':{'mode':'cdf'},}
	for figname,spec in figspec.items(): plot_radial_distributions(figname=figname,**spec)

@autoplot(plotrun)
def plot_cumulative_compare():
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
	legend = make_legend(axes[-1])
	picturesave('fig.lipid_rdfs.%s'%figname,directory=work.plotdir,meta={},extras=[legend])

plotrun.routine = None
if __name__=='__main__': pass

