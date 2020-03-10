#!/usr/bin/env python

"""
Plot lipid-lipid g(r) results.
This was forked from plot-lipid_rdfs.py to highlight convergence.
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
		if not all([self.this[0]==i for i in self.this[1:]]): 
			raise Exception('not same')
		else: 
			try: return self.this[0]
			except: return None

def color_by_simulation(*args,**kwargs): return 'k'

def attach_selector_key(meta):
	"""Include the abstractor_selector in the metadata."""
	try: 
		abstractor_selector = work.plotspec.specs['calculation'][
			'lipid_rdfs']['upstream']['lipid_abstractor']['selector']
		meta['abstractor_selector'] = abstractor_selector
	except: pass

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
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,
					fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	if legend_below:
		legend = fig.legend(patches,labels,
			loc='upper center',ncol=legend_ncol 
			if legend_ncol else len(legendspec),fontsize=14)
	else: legend = ax.legend(patches,labels,loc='upper left',
		bbox_to_anchor=(1.0,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor('w')
	frame.set_facecolor('w')
	return legend

def plot_radial_distributions(figname,mode='summary',**kwargs):
	"""Plot RDFs or CDFs of the RDFs."""
	# plot lipid-lipid RDFs on one panel
	reps_subset = kwargs.get('groups',None)
	pairings_this = kwargs.get('pairings',pairings)
	figsize = kwargs.pop('figsize',(12,12))
	legend_spec = kwargs.get('legend_spec',{})
	if reps_subset: replicate_mapping_this = [(i,j) 
		for i,j in replicate_mapping if i in reps_subset]
	else: replicate_mapping_this = replicate_mapping
	pair_residues_this = [pair_residues[ii] 
		for ii,i in enumerate(pairings) if i in pairings_this]
	axes,fig = square_tiles(len(pairings_this),
		figsize=figsize,hspace=0.4,wspace=0.4)
	for pnum,(pairname,pair) in enumerate(
		zip(pairings_this,pair_residues_this)):
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
				vals = np.array([
					scipy.integrate.simps(counts[:i]/areas[:i]/density) 
					for i in range(2,len(counts))])[valid]
				ax.plot(middles[valid],vals[valid],
					color=color_by_simulation(sns_this[0]),lw=2)
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

@autoload(plotrun)
def load():

	# SETTINGS
	sns = work.sns()
	#! assume zero for these monolayers
	monolayer_indexer = lambda sn,abstractor='lipid_chol_com': 0
	#! identity replicate mapping
	replicate_mapping = [(i,[i,]) for i in sns]
	# timescale for convergence (ns)
	#! set this to 10 for the 100ns simulations or 50 for the 500ns
	ts_converge = 10
	# plot settings
	xmax_cumulative = 4.0
	extra_labels = {'membrane-v565':'Ca2+','membrane-v563':'Mg2+',}
	dpi,img_form = 600,'pdf'

	if 'post' not in globals():
		# load the data
		data,calc = plotload('lipid_rdfs_detail')
		# handle loops with a Post object
		post = Post()
		# if the calculation is in the top level of calc then we have no loop
		post.new('lipid_rdfs',data=data,calc=calc['lipid_rdfs_detail'])
		# loop over upstream calculations
		for upstream_name in post.dataset:
			post.select(upstream_name)
			cutoff,binsize = [post.calc['calcs']['specs'][i] 
				for i in ['cutoff','binsize']]
			scanrange = np.arange(0,cutoff,binsize)
			middles = (scanrange[1:]+scanrange[:-1])/2
			areas = np.array([np.pi*(binsize)*middles[i]*2 
				for i in range(len(middles))])
			data_this = post.data['data']
			# collect all possible pairings
			pairings = list(set([tuple(i) for j in [
				data_this[sn]['data']['pairnames_mn%d'%
				actinlink_monolayer_indexer(sn,abstractor='lipid_chol_com')] 
				for sn in sns] for i in j]))
			# trawl all simulations for the pairs
			pairname_to_pair = {}
			for sn in sns:
				pairname_to_pair.update(**dict([(tuple(i),j) 
					for i,j in zip(*[data_this[sn]['data']['%s_mn%d'%(p,
					monolayer_indexer(sn,abstractor='lipid_chol_com'))] 
					for p in ['pairnames','pairs']])]))
			pair_residues = [pairname_to_pair[p] for p in pairings]
			#! see original plot-lipid_rdfs.py for a note on some magic
			post_this = post.data['reduced'] = dict([])
			for pnum,(pairname,pair) in enumerate(zip(pairings,pair_residues)):
				post_this[pairname] = {}
				for sn_group,sns_this in replicate_mapping:
					obs,total_areas,series = [],[],[]
					nmols,nsegments = Same(),Same()
					for sn in sns_this:
						mn = monolayer_indexer(sn,abstractor='lipid_chol_com')
						dat = data_this[sn]['data']
						key = 'counts_mn%d_%s:%s'%(mn,pairname[0],pairname[1])
						if key not in dat: 
							status('reducing %s, %s SKIP'%(sn,key))
							continue
						else: status('reducing %s, %s'%(sn,key))
						# get the number of frames
						nframes = dat[key].shape[0]
						duration = (calc['extras'][sn]['end'] -
							calc['extras'][sn]['start'])/(nframes-1)
						times = np.linspace(0,duration,nframes)
						# subselect the right indices
						inds = [m for m in [np.where(np.all([times<=j,times>=i],axis=0))[0] for i,j in [(ts_converge*(i-1),ts_converge*i) for i in range(1,duration/ts_converge+1)]] if len(m)>0]
						obs.append(dat[key].mean(axis=0))
						series_this = np.zeros((len(inds),middles.shape[0]))
						for ww,when in enumerate(inds):
							series_this[ww] = dat[key][when].mean(axis=0)
						series.append(series_this)
						nmols.add([np.in1d(
							dat['resnames'][dat['monolayer_indices']==mn],p).sum() 
							for p in pair])
						total_areas.append(dat['total_area'])
					# average replicates
					counts = np.array(obs).mean(axis=0)
					nmols = nmols.get()
					# not every simulation has every pairname
					if not nmols: continue
					total_area = np.mean(total_areas) #! total area is averaged between simulations for norming
					density = nmols[0]*nmols[1]/total_area
					post_this[pairname][sn_group] = dict(counts=counts,
						total_area=total_area,density=density,
						nmols=nmols,series=series)
			# load necessary globals into the post object for later
			#! note that quantities which vary with upstream data have to be stored in post otherwise
			#! ... they go directly to globals and may linger there during later loops over upstream data
			post.data['pairings'] = pairings
			post.data['replicate_mapping'] = replicate_mapping
			post.data['pair_residues'] = pair_residues
			post.data['middles'] = middles
			post.data['areas'] = areas
	
plotrun.routine = [] 

if __name__=='__main__': 

	# original set
	if 0:

		# use the "long" collection with two simulations
		assert len(calc['extras'].keys())==2
		assert ts_converge == 50
		sns_order = ['membrane-v563','membrane-v565']
		#! previous plot plot_distributions()
		pairings_this = post.data['pairings']
		# the following plot is based on `plot_distributions`
		pair_residues_this = [pair_residues[ii] 
			for ii,i in enumerate(pairings) if i in pairings_this]
		meta = {}
		if 0:
			axes,fig = panelplot(figsize=(20,10),
				layout={'out':{'grid':[1,len(pairings)],'wspace':0.8,},
				'ins':[{'grid':[len(sns),1],'hspace':0.3} for i in pairings]})
			ax_cb = axes[0][-1]
			ax_finder = lambda pnum,snum: axes[pnum][snum]
		else:
			axes,fig = square_tiles(len(pairings)*len(sns),
				figsize=(20,20),hspace=0.4,wspace=0.4)
			ax_finder = lambda pnum,snum: axes[pnum+len(pairings)*snum]
			ax_cb = axes[-1]
		for snum,sn in enumerate(sns_order):
			#!!! not set up for replicates. written in haste
			sn_group = sn
			for pnum,(pairname,pair) in enumerate(
				zip(pairings_this,pair_residues_this)):
				ax = ax_finder(pnum,snum)
				reduced = post.data['reduced']
				counts,nmols,total_area,density = [reduced[pairname][sn_group][k] 
					for k in ['counts','nmols','total_area','density']]
				#! show 90% of the maximum range
				xmax = 0.9*np.sqrt(total_area)
				valid = middles<xmax
				valid = middles<xmax
				this_series = reduced[pairname][sn_group]['series']
				if len(this_series)>1: raise Exception('no replicates for series')
				this_series = this_series[0]
				rainbow = [mpl.cm.__dict__['RdBu'](i) 
					for i in np.linspace(0,1,this_series.shape[0])]
				for series_num,series in enumerate(this_series):
					ax.plot(middles[valid],(series/areas/density)[valid],
						color=rainbow[series_num],lw=1)
				ax.axhline(1.0,lw=1,c='k')
				ax.set_xlim((0,xmax))
				ax.set_ylabel('$g(r)$')
				ax.axhline(1.0,lw=0.5,c='k')
				ax.set_xlim((0,xmax_cumulative))
				ax.set_title('%s %s - %s'%tuple([work.meta[sn]['ion_label']]+[
					work.vars['names']['short'].get(p,p) 
					for p in pairname]))
				ax.tick_params(axis='y',which='both',left='off',
					right='off',labelleft='on')
				ax.tick_params(axis='x',which='both',top='off',
					bottom='off',labelbottom='on')
				ax.set_xlabel(r'$r\,(nm)$')
		cmap = mpl.cm.__dict__['RdBu']
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		#! https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
		axins = inset_axes(ax_cb,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.15,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
		bounds = np.linspace(0,1,10+1)
		norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
		cb = mpl.colorbar.ColorbarBase(axins,
			cmap=cmap,norm=norm,
		    spacing='proportional',ticks=bounds,
		    boundaries=bounds,format='%1i')
		axins.set_title('time (ns)')
		axins.set_yticklabels(['%d '%i 
			for i in np.linspace(0,duration,len(series_this)+1)])
		legends = [axins]
		picturesave('fig.lipid_rdfs.convergence',directory=work.plotdir,
			version=True,meta=meta,extras=legends,dpi=dpi,form=img_form)

	# getting the all-lipids for all simulations
	# switch to collection "all" from collection "long" to get this to work
	if 1:

		assert len(calc['extras'].keys())>2
		assert ts_converge == 10
		sns_order = work.sns()
		pairings_this = [(u'all lipids', u'all lipids')]
		# the following plot is based on `plot_distributions`
		pair_residues_this = [pair_residues[ii] 
			for ii,i in enumerate(pairings) if i in pairings_this]
		meta = {}
		axes,fig = square_tiles(len(pairings_this)*len(sns),
			figsize=(20,20),hspace=0.4,wspace=0.4)
		ax_finder = lambda pnum,snum: axes[pnum+len(pairings_this)*snum]
		ax_cb = axes[-1]
		for snum,sn in enumerate(sns_order):
			#!!! not set up for replicates. written in haste
			sn_group = sn
			for pnum,(pairname,pair) in enumerate(
				zip(pairings_this,pair_residues_this)):
				ax = ax_finder(pnum,snum)
				reduced = post.data['reduced']
				counts,nmols,total_area,density = [reduced[pairname][sn_group][k] 
					for k in ['counts','nmols','total_area','density']]
				#! show 90% of the maximum range
				xmax = 0.9*np.sqrt(total_area)
				valid = middles<xmax
				valid = middles<xmax
				this_series = reduced[pairname][sn_group]['series']
				if len(this_series)>1: raise Exception('no replicates for series')
				this_series = this_series[0]
				rainbow = [mpl.cm.__dict__['RdBu'](i) 
					for i in np.linspace(0,1,this_series.shape[0])]
				for series_num,series in enumerate(this_series):
					ax.plot(middles[valid],(series/areas/density)[valid],
						color=rainbow[series_num],lw=1)
				ax.axhline(1.0,lw=1,c='k')
				ax.set_xlim((0,xmax))
				ax.set_ylabel('$g(r)$')
				ax.axhline(1.0,lw=0.5,c='k')
				ax.set_xlim((0,xmax_cumulative))
				ax.set_title('%s %s - %s'%tuple([work.meta[sn]['ion_label']]+[
					work.vars['names']['short'].get(p,p) 
					for p in pairname]))
				ax.tick_params(axis='y',which='both',left='off',
					right='off',labelleft='on')
				ax.tick_params(axis='x',which='both',top='off',
					bottom='off',labelbottom='on')
				ax.set_xlabel(r'$r\,(nm)$')
		cmap = mpl.cm.__dict__['RdBu']
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		#! https://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
		axins = inset_axes(ax_cb,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.15,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
		bounds = np.linspace(0,1,10+1)
		norm = mpl.colors.BoundaryNorm(bounds,cmap.N)
		cb = mpl.colorbar.ColorbarBase(axins,
			cmap=cmap,norm=norm,
		    spacing='proportional',ticks=bounds,
		    boundaries=bounds,format='%1i')
		axins.set_title('time (ns)')
		axins.set_yticklabels(['%d '%i 
			for i in np.linspace(0,duration,len(series_this)+1)])
		legends = [axins]
		picturesave('fig.lipid_rdfs.convergence.all',directory=work.plotdir,
			version=True,meta=meta,extras=legends,dpi=dpi,form=img_form)
