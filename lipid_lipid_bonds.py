#!/usr/bin/python env

"""
LIPID-LIPID hydrogen bond and salt-bridge analysis
supercedes lipid_lipid_bonds_analysis
"""

import seaborn as sb
import pandas as pd
import itertools,collections,copy

### LOAD

def load_generic(**kwargs):
	"""
	Generic load function which can be customized for your project 
	by writing a "load_<project>" function corresponding to the parent folder.
	"""
	data = kwargs['data']
	# simulations listing from the (ordered) collections in the metadata plots entry for this plot
	sns = work.sns()
	# each simulation gets a color
	colors = [mpl.cm.__dict__['jet'](i) for i in np.linspace(0,1.,len(sns))]
	def color_by_simulation(sn):
		return colors[sns.index(sn)]
	# multiple replicates can be merged later
	replicate_mapping = [(sn,[sn]) for sn in sns]
	# cleaner labels in case your "sn" simulation names are clumsy
	extra_labels = dict([(sn,sn) for sn in sns])
	# counting lipids is necessary to normalize them
	nmol_counts = dict([(sn,dict(zip(*[data['hydrogen_bonding'][sn]['data'][k] 
		for k in ['resnames','nmols']]))) for sn in sns])
	# elegant names for residues e.g. PtdIns
	residue_renamer = lambda x:x
	# exclude certain residue names
	resnames_exclude = []
	# whether to run the calculation with merged simulations
	merged_opts = [False,True][:1]
	# control the residue ordering on the grid-style summary plots
	key_residues = []
	# we can report bonds per lipid on both axes if compositions are uniform
	bonds_per_lipid = True
	# print the super title above the plots
	title_style = ['suptitle',None][0]
	# legend can be above the plots or on the right, for longer lists
	legend_position = ['upper center','upper right'][0]
	# return variables which are exported to globals in the calling function
	outgoing = dict(sns=sns,color_by_simulation=color_by_simulation,
		replicate_mapping=replicate_mapping,extra_labels=extra_labels,nmol_counts=nmol_counts,
		residue_renamer=residue_renamer,resnames_exclude=resnames_exclude,key_residues=key_residues,
		merged_opts=merged_opts,bonds_per_lipid=bonds_per_lipid,
		title_style=title_style,legend_position=legend_position)
	return outgoing

### PROJECT CUSTOMIZATIONS

def load_ptdins(**kwargs):
	"""Customizations for the PtdIns project."""
	outgoing = load_generic(**kwargs)
	# custom definitions
	sns = work.sns()
	def residue_renamer(resname):
		if resname in ['PI2P','P35P','PIPP','PIPU','SAPI']: return 'PtdIns'
		else: return resname
	def color_by_simulation(sn):
		return colorize(work.meta[sn],comparison='asymmetric_all')
	extra_labels = dict([(sn,'%s\n%s'%(
		work.meta.get(sn,{}).get('ion_label','ion'),
		work.meta.get(sn,{}).get('ptdins_label','ptdins'),
		)) for sn in sns])
	# we exclude POPC and halve the number of apparent cholesterol in order to account for leaflets
	# note that an automatic accounting of leaflet differences would require more code and leaflet indices
	resnames_exclude = ['POPC']
	nmol_counts = outgoing['nmol_counts']
	for sn in sns: 
		nmol_counts[sn]['PtdIns'] = nmol_counts[sn][work.meta[sn]['ptdins_resname']]
		if 'CHL1' in nmol_counts[sn]: nmol_counts[sn]['CHL1'] = nmol_counts[sn]['CHL1']/2.
	merged_opts = [False,True][:1]
	key_residues = ['PtdIns','DOPS','DOPE','CHL1'][::-1]
	# repack
	outgoing.update(residue_renamer=residue_renamer,color_by_simulation=color_by_simulation,
		extra_labels=extra_labels,resnames_exclude=resnames_exclude,key_residues=key_residues,
		merged_opts=merged_opts)
	return outgoing

def load_actinlink(**kwargs):
	"""Customizations for the actinlink project."""
	outgoing = load_generic(**kwargs)
	# normalize by both lipid species when composition varies
	outgoing['bonds_per_lipid'] = False
	resnames_exclude = ['POPC']
	merged_opts = [False,True][:]
	key_residues = ['PI2P','DOPS','DOPE','CHL1']
	replicate_mapping = actinlink_replicate_mapping
	color_by_simulation = actinlink_color_by_simulation
	extra_labels = actinlink_extra_labels
	title_style = None
	legend_position = 'upper right'
	outgoing.update(resnames_exclude=resnames_exclude,merged_opts=merged_opts,
		key_residues=key_residues,replicate_mapping=replicate_mapping,
		color_by_simulation=color_by_simulation,extra_labels=extra_labels,
		title_style=title_style,legend_position=legend_position)
	return outgoing

### PLOT

class PlotBonds:
	"""
	Handle bond counting and plotting.
	"""
	def __init__(self,**kwargs):
		figsize = kwargs.pop('figsize',16.)
		wspace = kwargs.pop('wspace',0.5)
		hspace = kwargs.pop('hspace',0.4)
		dims = kwargs.pop('dims',(0,1))
		self.name = kwargs.pop('name')
		self.dataspec = kwargs.pop('dataspec',{})
		self.style = self.dataspec.pop('style','violin')
		self.legend_pad = self.dataspec.pop('legend_pad',3.0)
		self.key_residues = self.dataspec.pop('key_residues',[])
		# per the hydrogen_bonding.py, the hydrogen bond donor is listed first
		self.donor_acceptor = self.dataspec.pop('donor_acceptor',True)
		self.title_style = self.dataspec.pop('title_style','suptitle')
		self.normed = self.dataspec.pop('normed',None)
		self.resnames_exclude = self.dataspec.pop('resnames_exclude',[])
		self.symmetric = self.dataspec.pop('symmetric',False)
		self.legend_position = self.dataspec.pop('legend_position','upper center')
		if dims not in [(0,1),(1,0)]: raise Exception
		# the "special" plot method uses a grid to plot the lipid types in the pair
		if self.dataspec.pop('special',False):
			self.style = 'violin_special'
			self.symmetric = self.dataspec.pop('special_symmetric',True)
		self.raw = self.count_bonds(**self.dataspec)
		if kwargs: raise Exception(kwargs)
		self.extras,self.dump = [],[]
		# unpack into tiers by dims
		self.data = collections.OrderedDict()
		self.index = collections.OrderedDict()
		for k1,v1 in self.raw.items():
			for k2,v2 in v1.items():
				keys = tuple([k1,k2][i] for i in dims)
				if keys[0] not in self.index:
					self.index[keys[0]] = []
				self.index[keys[0]].append(keys[1])
				self.data[tuple([k1,k2][i] for i in dims)] = v2
		self.n_tiles = len(self.index)
		if self.style=='violin_special':
			res = list(set([i for j in self.index.keys() for i in j]))
			panelspec = dict(figsize=figsize if type(figsize) in [tuple,list] else (figsize,figsize),
				layout={'out':{'grid':[len(res),1],'wspace':wspace,'hspace':hspace},
				'ins':[{'grid':[1,len(res)],'wspace':wspace,'hspace':hspace} for i in res]})
			self.axes,self.fig = panelplot(**panelspec)
		else:
			self.axes,self.fig = square_tiles(
				self.n_tiles,figsize=figsize,wspace=wspace,hspace=hspace,favor_rows=True)
		self.plot(self.style)	
		picturesave('fig.lipid_lipid_bonds.%s'%self.name,work.plotdir,
			backup=False,version=True,meta={},extras=self.extras)

	def count_bonds(self,kinds,merged=False):
		"""
		Post-process the bond counts.
		"""
		#! replace with a mapping function to handle cholesterol in both leaflets
		nframes_by_sn = dict([(sn,len(data['salt_bridges'][sn]['data']['observations'])) for sn in sns])
		if merged: self.replicate_mapping_this = replicate_mapping 
		else: self.replicate_mapping_this = [(sn,[sn]) for sn in sns]
		post = {}
		for sn_group,sns_this in self.replicate_mapping_this:
			for snum,sn in enumerate(sns_this):
				for knum,kind in enumerate(kinds):
					# collect bonds
					bonds,obs = [data[kind][sn]['data'][k] for k in ['bonds','observations']]
					# collect pairs
					if len(bonds)==0: continue
					subjects = np.unique(bonds[:,rowspec.index('subject_resname')])
					targets = np.unique(bonds[:,rowspec.index('target_resname')])
					resnames = [i for i in np.unique(np.concatenate((subjects,targets)))
						if i not in self.resnames_exclude]
					combos = list(itertools.product(resnames,repeat=2))
					if residue_renamer!=None: 
						combos_renamed = [tuple([residue_renamer(i) for i in j]) for j in combos]
					else: combos_renamed = combos
					for combo,combo_renamed in zip(combos,combos_renamed):
						rows = np.where(np.all((
							bonds[:,rowspec.index('subject_resname')]==combo[0],
							bonds[:,rowspec.index('target_resname')]==combo[1],
							bonds[:,rowspec.index('subject_resid')]!=bonds[:,rowspec.index('target_resid')],
							),axis=0))[0]
						if len(rows)>0:
							if sn_group not in post: post[sn_group] = {}
							if combo_renamed not in post[sn_group]: 
								# include zeros for cases where a rare bond appears in only one replicate
								post[sn_group][combo_renamed] = [[np.zeros(nframes_by_sn[ii]) 
									for ii in sns_this] for jj in kinds]
							counts = obs[:,rows].sum(axis=1)
							post[sn_group][combo_renamed][knum][snum] = counts	
		# post-processing to sum bond types (first dimension) then concatenate replicates (second)
		for sn_group in post:
			for combo in post[sn_group]:
				counts_by_rep = []
				for rep in zip(*post[sn_group][combo]):
					counts_by_bond = []
					for i in rep:
						if len(i)>0: counts_by_bond.append(i)
					cbb = np.array(counts_by_bond)
					if len(cbb.shape)>1: cbb = cbb.sum(axis=0)
					counts_by_rep.append(cbb)
				post[sn_group][combo] = np.concatenate(counts_by_rep)
		if self.symmetric:
			post_symmetric = {}
			for sn in post:
				post_symmetric[sn] = {}
				combos_sym = list(set([tuple(sorted(i)) for i in post[sn].keys()]))
				for combo in combos_sym:
					post_symmetric[sn][combo] = np.sum([post[sn][c] 
						for c in post[sn] if set(c)==set(combo)],axis=0)
			post = post_symmetric
		# save a self-specific ordered list of simulation names
		self.sns = list(zip(*self.replicate_mapping_this)[0])
		return post

	def plot(self,style):
		"""
		A note on reading the y-axis. It counts bonds per <the left/first/corresponding lipid> with the other.
		So if you have a value of 0.2 on the left y-axis for PtdIns-DOPE and it is 0.04 on the right, this 
		means that there are 0.2 bonds for each PtdIns to a DOPE or 0.04 bonds for each DOPE to a PtdIns.
		Note also that we plot in the order donor-acceptor but you could just as easily reverse this. If you
		reverse it, then the y-axes are scaled differently depending on the composition, even if the actual
		distributions do not change when the plot is symmetric. We retain the option to reverse the order 
		as a consistency check only.
		"""
		if self.style=='violin_special':
			def get_nmols(top,strict=False):
				nmols = dict([(sn,tuple([float(nmol_counts[sn].get(c,0)) for c in top])) for sn in sns])
				if not strict: return [nmols[sn] for sn in sns]
				if len(nmols)!=1 and strict: 
					raise Exception('cannot report bonds per lipid when compositions vary')
				else: nmols = nmols[0]
				return nmols
			y_max = {}
			if self.key_residues: res = self.key_residues 
			else: res = list(set([i for j in self.index.keys() for i in j]))
			for ii,r1 in enumerate(res):
				for jj,r2 in enumerate(res):
					top = (r1,r2) if self.donor_acceptor else (r2,r1)
					ax = self.axes[ii][jj]
					if (self.symmetric and ii>jj) or (
						top not in self.index and top[::-1] not in self.index):
						ax.set_yticklabels([])
						if self.symmetric: 
							for spine in ax.spines.values(): spine.set_edgecolor('w')
							ax.set_yticks([])
						ax.set_xticks([])
						ax.set_xticklabels([])
						continue
					# custom titles come from variables-names-short in the metadata
					ax.set_title(r'%s$-$%s'%tuple([
						work.vars.get('names',{}).get('short',{}).get(p,p) for p in top]))
					top_lookup = top if top in self.index else top[::-1]
					nmols = get_nmols(top,strict=bonds_per_lipid)
					if bonds_per_lipid: valid_indices = range(len(self.sns))
					else: valid_indices = [si for si,sn in enumerate(self.sns)
						if all([nmols[si][j]>0 for j in range(2)])]
					# data are ordered by simulation
					colors = [color_by_simulation(sn) for si,sn in enumerate(self.sns) if si in valid_indices]
					palette = dict(zip(range(len(self.index)),colors))
					"""
					we norm at the last moment before plotting the data instead of norming when we count them
					we norm by the first species so you read the plot as "first lipid has <left axis> bonds
						with the second lipid" or "second lipid has <right axis> number of bonds with the 
						first lipid"
					if we reverse the direction of the donor-acceptor order then we have to flip the 
						normalization as well otherwise the axes will be wrong and the whole point
						of plotting both directions (admittedly a pedantic, completist effort) will be lost
					"""
					if not self.normed: factor_adjust = [1. for sn in self.sns]
					elif bonds_per_lipid: factor_adjust = [1./nmols[0] for sn in self.sns]
					else: factor_adjust = [1./nmols[si][0]/nmols[si][1] 
						for si,sn in enumerate(self.sns) if si in valid_indices]
					# carefully choose distributions in the canonical ordering by sns
					data_this = [self.data[(top_lookup,sn)]*factor_adjust[valid_indices.index(si)] 
						if (top_lookup,sn) in self.data else [] for si,sn in enumerate(self.sns)
						if si in valid_indices]
					outgoing = dict(data=data_this)
					vp = sb.violinplot(ax=ax,bw=1.0,cut=0,linewidth=0.0,palette=palette,
						labels=[work.meta.get(sn,{}).get('ptdins_label',sn) for sn in self.sns],
						**outgoing)
					ax.set_xticks([])
					ax.set_xticklabels([])
					try: y_max[(r1,r2)] = max([i.max() for i in data_this if len(i)>0])
					except: continue
					self.dump.append(vp)
			# standardize the y-axis maxima
			y_max_val = max(y_max.values())
			for ii,r1 in enumerate(res):
				for jj,r2 in enumerate(res):
					top = (r1,r2) if self.donor_acceptor else (r2,r1)
					if (self.symmetric and ii>jj) or (
						top not in self.index and top[::-1] not in self.index):
						continue
					ax = self.axes[ii][jj]
					ax.set_ylim((0,y_max_val))
					if not bonds_per_lipid or not self.normed: continue
					nmols = get_nmols(top,strict=bonds_per_lipid)
					twin = ax.twinx()
					yticks = ax.get_yticks()
					factors_rescale = [float(nmols[0])/nmols[1] for sn in self.sns]
					# the twin does not use the correct formatting so we just fix it here
					twin.set_yticklabels(['%.2f'%i for i in yticks*factors_rescale])
					twin.set_yticks(yticks)
					twin.set_ylim((0,y_max_val))
					ax.set_yticklabels(['%.2f'%i for i in yticks])
		else: raise Exception('invalid or deprecated style %s'%self.style)
		title_names = {'hydrogen_bonding':'hydrogen bonds','salt_bridges':'salt bridges'}
		self.title = ' and '.join([title_names[k] for k in self.dataspec['kinds']])
		if self.title_style=='suptitle':
			self.extras.append(plt.suptitle(self.title,fontsize=20,y=1.0))
		if self.legend_position=='upper right': self.make_legend(ax=self.axes[0][-1],ncol=1)
		# default legend placement is above the panels
		elif self.legend_position=='upper center': self.make_legend(ax=None)
		else: raise Exception

	def make_legend(self,ax=None,keys=None,**kwargs):
		"""Make a legend. Note that this is redundant with plot-lipid_rdfs.py."""
		ncol = kwargs.get('ncol',None)
		legendspec = []
		for sn_general,sns in self.replicate_mapping_this:
			if not keys or sn_general in keys:
				legendspec.append(dict(
					name=extra_labels.get(sn_general,work.meta.get(sn_general,{}).get('label','')),
					patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		title_this = self.title if self.title_style=='legend' else None
		if ax!=None: 
			legend = ax.legend(patches,labels,loc='upper left',
				ncol=ncol if ncol else len(legendspec),fontsize=14,bbox_to_anchor=(1.05,1.0),
				borderaxespad=0)
		else: 
			legend = self.fig.legend(patches,labels,title=title_this,
				loc='upper center',ncol=ncol if ncol else len(legendspec),fontsize=14,
				borderaxespad=self.legend_pad)
		plt.setp(legend.get_title(),fontsize='xx-large')
		frame = legend.get_frame()
		frame.set_edgecolor('k')
		frame.set_facecolor('w')
		self.extras.append(legend)

### AUTOPLOT

@autoload(plotrun)
def load():
	"""
	Load once and export settings to globals.
	"""
	plotname = 'lipid_lipid_bonds'
	data,calc = work.plotload(plotname)
	# format of the bonds data 
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	# load switch by the directory / project name
	project_name = os.path.basename(os.getcwd())
	# write to locals so the Observer can export to globals
	load_function = 'load_%s'%project_name if 'load_%s'%project_name in globals() else 'load_generic'
	# custom load function to a special variable sent to globals by the autoplot Observer
	__locals__ = globals()[load_function](data=data)


@autoplot(plotrun)
def plot():
	"""Plot hydrogen bonds with many different accounting options."""
	sb.set_style("white")
	# different types of bonds
	kind_map = {
		'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	figspec = {(key+
		{'violin':'','violin_special':''}[style]+
		{1:'.normed',0:''}[normed]+
		{1:'.merged',0:''}[merged]+
		{1:'.symmetric',0:''}[symmetric]):{
			'merged':merged,'kinds':kind,'normed':normed,'resnames_exclude':resnames_exclude,
			'key_residues':key_residues,'symmetric':symmetric,'legend_position':legend_position,
			'title_style':title_style,'style':style}
			for style in ['violin_special']
			for symmetric in [True,False]
			for normed in [True,False]
			for merged in merged_opts
			for key,kind in kind_map.items()}
	# plot with acceptors first as a consistency check
	figspec['hbonds.normed.acceptor_donor'] = dict(donor_acceptor=False,
		**copy.deepcopy(figspec['hbonds.normed']))
	test_key = ['hbonds_salt.merged.symmetric',None][-1]
	if test_key: figspec = {test_key: figspec[test_key]}
	for name,spec in figspec.items():
		pb = PlotBonds(dims=(1,0),name=name,dataspec=spec)

plotrun.routine = None
if __name__=='__main__': pass
