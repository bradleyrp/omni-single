#!/usr/bin/python env

"""
LIPID-LIPID hydrogen bond and salt-bridge analysis
supercedes lipid_lipid_bonds_analysis
"""

import seaborn as sb
import pandas as pd
import itertools,collections,copy
mpl.rcParams['hatch.linewidth'] = 1.5
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import glob

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
		self.sns = kwargs.pop('sns_custom',sns)
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
		self.custom_axes = kwargs.pop('custom_axes',{})
		self.y_axis_align = kwargs.pop('y_axis_align',True)
		self.y_axis_label = kwargs.pop('y_axis_label',False)
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
		if self.custom_axes: 
			self.fig = self.custom_axes['fig']
			self.axes_mapping = self.custom_axes['axes_mapping']
			self.axes = self.axes_reorder = self.custom_axes['axes_reorder']
		elif self.style in ['violin_special','bars']:
			res = list(set([i for j in self.index.keys() for i in j]))
			panelspec = dict(figsize=figsize if type(figsize) in [tuple,list] else (figsize,figsize),
				layout={'out':{'grid':[len(res),1],'wspace':wspace,'hspace':hspace},
				'ins':[{'grid':[1,len(res)],'wspace':wspace,'hspace':hspace} for i in res]})
			self.axes,self.fig = panelplot(**panelspec)
		else:
			self.axes,self.fig = square_tiles(
				self.n_tiles,figsize=figsize,wspace=wspace,hspace=hspace,favor_rows=True)
		self.plot(self.style)	
		meta = dict(style=self.style)
		if 'salt_bridges' in self.dataspec['kinds']:
			meta['salt_bridge_cutoff'] = calc['salt_bridges']['calcs']['specs']['distance_cutoff']
		if self.custom_axes: return
		else: picturesave('fig.lipid_lipid_bonds.%s'%self.name,work.plotdir,form='pdf',
			backup=False,version=True,meta=meta,extras=self.extras)

	def count_bonds(self,kinds,merged=False):
		"""
		Post-process the bond counts.
		"""
		#! replace with a mapping function to handle cholesterol in both leaflets
		nframes_by_sn = dict([(sn,len(data['salt_bridges'][sn]['data']['observations'])) for sn in self.sns])
		if merged: self.replicate_mapping_this = replicate_mapping 
		else: self.replicate_mapping_this = [(sn,[sn]) for sn in self.sns]
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
						#! debugging the previous code at plot-hydrogen_bonding.py 
						#! ... rows_with_intra = np.where(np.all((bonds[:,rowspec.
						#! ... index('subject_resname')]==combo[0],
						#! ... bonds[:,rowspec.index('target_resname')]==combo[1],),axis=0))[0]
						if len(rows)>0:
							if sn_group not in post: post[sn_group] = {}
							if combo_renamed not in post[sn_group]: 
								# include zeros for cases where a rare bond appears in only one replicate
								post[sn_group][combo_renamed] = [[np.zeros(nframes_by_sn[ii]) 
									for ii in sns_this] for jj in kinds]
							counts = obs[:,rows].sum(axis=1)
							post[sn_group][combo_renamed][knum][snum] = counts	
							"""
							debugging detritus
							if sn=='membrane-v538':
								print (combo)
								if combo==('PI2P','PI2P'):
									import ipdb;ipdb.set_trace()
									rows_with_intra = np.where(np.all((
										bonds[:,rowspec.index('subject_resname')]==combo[0],
										bonds[:,rowspec.index('target_resname')]==combo[1],),axis=0))[0]
									obs[:,rows_with_intra].sum(axis=1).mean()
							"""
		# post-processing to sum bond types (first dimension) then concatenate replicates (second)
		for sn_group in post:
			for combo in post[sn_group]:
				counts_by_rep = []
				for rep in zip(*post[sn_group][combo]):
					counts_by_bond = []
					for i in rep:
						if len(i)>0: counts_by_bond.append(i)
					cbb = np.array(counts_by_bond)
					# if there is more than one bond type, we sum them here
					if cbb.shape[0]>1: 
						# legacy PtdIns calculations had rare missing frames due to the infamous ckDTree
						# ... bug so we save the valid frames and only consider those
						try: cbb = cbb.sum(axis=0)
						except:
							"""
							# tabulate nframes by simulation
							#nframes_by_sn = dict([(sn,len(data['salt_bridges'][sn]['data']['observations'])) for sn in sns])
							valid_frames_by_kind = dict([(sn,dict([(k,data[k][sn]['data']['valid_frames']) for k in kinds])) for sn in sns])
							"""
							import ipdb;ipdb.set_trace()
					counts_by_rep.append(cbb)
				# concatenate replicates here
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
		#! self.sns = list(zip(*self.replicate_mapping_this)[0])
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
		if self.style=='bars':
			bar_formats = make_bar_formats(sns,work)
			palette_bars = dict(zip(range(len(self.index)),[bar_formats[sn]['c'] for sn in self.sns]))
		if self.style in ['violin_special','bars']:
			def get_nmols(top,strict=False):
				nmols = dict([(sn,tuple([float(nmol_counts[sn].get(
					c if c!='PtdIns' else work.meta[sn]['ptdins_resname'],0)) 
					for c in top])) for sn in self.sns])
				if not strict: return [nmols[sn] for sn in self.sns]
				if len(set(nmols.values()))!=1 and strict: 
					raise Exception('cannot report bonds per lipid when compositions vary')
				else: nmols = list(set(nmols.values()))[0]
				return nmols
			y_max = {}
			if self.key_residues: res = self.key_residues 
			else: res = list(set([i for j in self.index.keys() for i in j]))
			for ii,r1 in enumerate(res):
				for jj,r2 in enumerate(res):
					top = (r1,r2) if self.donor_acceptor else (r2,r1)
					# custom axis handling
					if self.custom_axes:
						if (r1,r2) not in self.axes_mapping: continue
						else: ax = self.axes[self.axes_mapping[(r1,r2)]]
					else: ax = self.axes[ii][jj]
					if (self.symmetric and ii>jj) or (
						top not in self.index and top[::-1] not in self.index):
						ax.set_yticklabels([])
						if self.symmetric: 
							for spine in ax.spines.values(): spine.set_edgecolor('w')
							ax.set_yticks([])
						ax.set_xticks([])
						ax.set_xticklabels([])
						continue
					# removing ticks
					ax.tick_params(axis=u'both',which=u'both',length=0)
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
					if self.style=='violin':
						vp = sb.violinplot(ax=ax,bw=1.0,cut=0,linewidth=0.0,
							palette=palette,labels=[work.meta.get(sn,{}).get('ptdins_label',sn) 
								for sn in self.sns],**outgoing)
					elif self.style=='bars':
						vp = sb.barplot(ax=ax,palette=palette_bars,ci='sd',capsize=0.,**outgoing)
						for bnum,bar in enumerate(vp.patches):
							# very minor error here: used sns, not self.sns and it messed up the hatching
							bar.set_hatch(bar_formats[self.sns[bnum]]['hatch'])
					else: raise Exception
					ax.set_xticks([])
					ax.set_xticklabels([])
					try: y_max[(r1,r2)] = max([i.max() for i in data_this if len(i)>0])
					except: 
						print('failing to add %s'%((r1,r2)))
						continue
					self.dump.append(vp)
			# standardize the y-axis maxima
			#! hacking for now
			y_max_val = max(y_max.values())
			for ii,r1 in enumerate(res):
				for jj,r2 in enumerate(res):
					top = (r1,r2) if self.donor_acceptor else (r2,r1)
					if (self.symmetric and ii>jj) or (
						top not in self.index and top[::-1] not in self.index):
						continue
					#! two options from a merge
					if 1:
						# custom axis handling
						if self.custom_axes:
							if (r1,r2) not in self.axes_mapping: continue
							else: ax = self.axes[self.axes_mapping[(r1,r2)]]
						else: ax = self.axes[ii][jj]
						if self.y_axis_align: ax.set_ylim((0,y_max_val))
					else:
						ax = self.axes[ii][jj]
						ax.set_ylim((0,y_max_val))
						if not self.normed: ax.set_ylabel('bonds')
					if not bonds_per_lipid or not self.normed: continue
					nmols = get_nmols(top,strict=bonds_per_lipid)
					twin = ax.twinx()
					yticks = ax.get_yticks()
					factors_rescale = [float(nmols[0])/nmols[1] for i in yticks]
					# the twin does not use the correct formatting so we just fix it here
					if self.y_axis_align:
						twin.set_yticklabels(['%.2f'%i for i in yticks*factors_rescale])
						twin.set_yticks(yticks)
						twin.set_ylim((0,y_max_val))
					else: twin.set_yticks([])
					ax.set_yticklabels(['%.2f'%i for i in yticks])
					if self.y_axis_label: ax.set_ylabel('bonds per lipid')
		else: raise Exception('invalid or deprecated style %s'%self.style)
		title_names = {'hydrogen_bonding':'hydrogen bonds','salt_bridges':'salt bridges'}
		self.title = ' and '.join([title_names[k] for k in self.dataspec['kinds']])
		if self.title_style=='suptitle':
			self.extras.append(plt.suptitle(self.title,fontsize=20,y=1.0))
		if self.legend_position=='upper right': 
			self.make_legend(ax=self.axes[0][-1],ncol=1,fancy=self.style=='bars')
		# default legend placement is above the panels
		elif self.legend_position=='upper center': self.make_legend(ax=None,fancy=self.style=='bars')
		elif self.legend_position==None: pass
		else: raise Exception

	def make_legend(self,ax=None,keys=None,fancy=False,**kwargs):
		"""Make a legend. Note that this is redundant with plot-lipid_rdfs.py."""
		if fancy: bar_formats = make_bar_formats(work.sns(),work)
		ncol = kwargs.get('ncol',None)
		legendspec = []
		for sn_general,sns in self.replicate_mapping_this:
			if not keys or sn_general in keys:
				if not fancy:
					legendspec.append(dict(
						name=extra_labels.get(sn_general,work.meta.get(sn_general,{}).get('label','')),
						patch=mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,fc=color_by_simulation(sns[0]))))
				else:
					#! fancy bars only for ptdins, which does not have replicates
					legendspec.append(dict(
						name=extra_labels.get(sn_general,work.meta.get(sn_general,{}).get('label','')),
						patch=mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,
							fc=bar_formats[sn_general]['c'],hatch=bar_formats[sn_general]['hatch'])))
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		title_this = self.title if self.title_style=='legend' else None
		extras = dict(handleheight=2.0)
		if ax!=None: 
			legend = ax.legend(patches,labels,loc='upper left',
				ncol=ncol if ncol else len(legendspec),fontsize=14,bbox_to_anchor=(1.05,1.0),
				borderaxespad=0,**extras)
		else: 
			legend = self.fig.legend(patches,labels,title=title_this,
				loc='upper center',ncol=ncol if ncol else len(legendspec),fontsize=14,
				borderaxespad=self.legend_pad,**extras)
		plt.setp(legend.get_title(),fontsize='xx-large')
		frame = legend.get_frame()
		frame.set_edgecolor('k')
		frame.set_facecolor('w')
		self.extras.append(legend)
		return legend

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
def plot_counts():
	"""Plot hydrogen bonds with many different accounting options."""
	#! note that this was an autoplot then it was main and now it is autoplot so we can do snapshots
	do_bars = os.path.basename(os.getcwd())=='ptdins'
	sb.set_style("white")
	# different types of bonds
	kind_map = {
		'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		#! off because ckdtree-bug-induced frame problems 'hbonds_salt':['hydrogen_bonding','salt_bridges']
		}
	figspec = {(key+
		{'violin_special':'','bars':''}[style]+
		{1:'.normed',0:''}[normed]+
		{1:'.merged',0:''}[merged]+
		{1:'.symmetric',0:''}[symmetric]):{
			'merged':merged,'kinds':kind,'normed':normed,'resnames_exclude':resnames_exclude,
			'key_residues':key_residues,'symmetric':symmetric,'legend_position':legend_position,
			'title_style':title_style,'style':style}
			#! do not use bars for any project other than ptdins
			for style in (['violin_special']+(['bars'] if do_bars else []))[:]
			for symmetric in [True,False]
			for normed in [True,False]
			for merged in merged_opts
			for key,kind in kind_map.items()}
	# plot with acceptors first as a consistency check
	#! figspec['hbonds.normed.acceptor_donor'] = dict(donor_acceptor=False,
	#! 	**copy.deepcopy(figspec['hbonds.normed']))
	#! keeping two lines from a merge here
	#! test_key = [None,'hbonds.normed'][0]
	#! test_key = ['hbonds_salt.merged.symmetric',None,'hbonds.symmetric'][-1]
	#! if test_key: figspec = {test_key: figspec[test_key]}
	for name,spec in figspec.items():
		status('plotting %s'%name,tag='plot')
		pb = PlotBonds(dims=(1,0),name=name,dataspec=spec)

def manuscript_plot():
	"""
	Special plot that combines some of the lipid-lipid bond counts with snapshots.
	Note that this was easier than adding the data from this script to the snapshots over in 
	plot-hydrogen_bonding.py
	"""

	import matplotlib.image as mpimg
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

	#! hacking over to hydrogen_bonding_standard. a rare connection between plots that should be removed
	#! ... then folded settings into the specs for this plot script
	if False:
		layout_name = 'layout_2'
		layout = copy.deepcopy(work.metadata.plots[
			'hydrogen_bonding_standard']['specs']['snapshot_examples'][layout_name])
		interesting = work.metadata.plots['hydrogen_bonding_standard']['specs']['interesting']
	"""
	incrementing the name
	2b: original hydrogen bonds, however
	2c: fixed the hatching
	2d: added salt bridge bar plots
	"""
	rows_extras = ['hbonds','salt'][:1]
	layout_name = 'layout_2d' # went from layout_2b to layout_2c to correct hatching problem on the former
	#! actually we wstill need "interesting" from earlier
	interesting = work.metadata.plots['hydrogen_bonding_standard']['specs']['interesting']
	arrangement = 'square_tiles_bars_below'
	mesh_inset = True

	# prepare axes based on layout
	if len(rows_extras)==1:
		figsize = [16,14]
		ntiles = 16
		axes,fig = square_tiles(ntiles=ntiles,figsize=(12,12),wspace=0.5,favor_rows=True)
		axes_snapshots = axes[4:] # bars on top
		#! hacked to get this stuff in the right order sorry
		axes_snapshots = [axes_snapshots[i] for i in [0,1,4,5,8,9]+[2,3,6,7,10,11]]
		axes_counts = [axes[:4]]
	elif set(['hbonds','salt'])==set(rows_extras):
		figsize = [18,14]
		nrows_snap = 3
		nrows_counts = 2
		nrows = nrows_snap+nrows_counts
		ncols = 4
		axes,fig = panelplot(figsize=figsize,layout={'out':{'grid':[nrows,1],'hspace':0.5},'ins':[{'grid':[1,ncols],'wspace':0.5} for i in range(nrows)],})
		axes_counts = [axes[0],axes[1]]
		axes_snapshots = [m for n in [axes[i] for i in [2,3,4]] for m in n]
	else: raise Exception

	legends = []
	matches = work.metadata.plots['lipid_lipid_bonds']['specs']['images']
	tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
	nsnaps = [len(v) for i in matches for k,v in i.items()]
	# use the hybrid list/dict structure in the YAML hence the weird structure below
	for pnum,item in enumerate(matches):
		if len(item)!=1: raise Exception('use list/dict in YAML please')
		(iname,fns) = item.items()[0]
		# first axis is the bar plot
		#! hacked
		axnum = pnum 
		ax = axes_snapshots[axnum]
		ispec = interesting[iname]
		comparison = ispec['comparison']
		bond_type = ispec['bond_type']
		#! sns_this = work.metadata.collections[comparison]
		# sometimes the legend causes a NoneType problem on get_window_extent
		#   legend = render_interesting_hydrogen_bonds(ispec=ispec,postpreps=postprep(sns),
		# sns=sns,fig=fig,ax=ax,save=False,do_legend=pnum==len(matches)-1,style='alt')
		#! discarded: if pnum==0: legends.append(legend)
		for fnum,fn in enumerate(fns):
			#---second axis holds the exemplars
			if arrangement=='square_tiles': axnum = len(matches)+sum(nsnaps[:pnum])+fnum
			elif arrangement=='square_tiles_bars_below': axnum = sum(nsnaps[:pnum])+fnum
			ax = axes_snapshots[axnum]
			# previously sent a single image to sidestack but the following routine makes it square
			image = mpimg.imread(os.path.join(work.plotdir,fn))
			border = get_blank_border(image)
			sizes = np.array([np.ptp(i) for i in border])
			padding = np.array([(max(sizes)-s)/2. for s in sizes]).astype(int)
			zooms = np.array([[border[i][0]-1*padding[i],border[i][1]+padding[i]] 
				for i in range(2)])
			# note that you have to reverse the order of dimensions here
			image_zoomed = image[zooms[1][0]:zooms[1][1],zooms[0][0]:zooms[0][1]]
			ax.imshow(image_zoomed)
			if False:
				tb_ion = ax.text(-0.1,0.2,chr(ord('A')+len(legends)),fontsize=14,
					bbox=tagbox,rotation=0,ha="center",va="top",color='k',
					transform=ax.transAxes)
			else:
				tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.35")
				label_ion = '$\\mathrm{{%s}^{2+}}$'%('Ca' if axnum<6 else 'Mg')
				tb_ion = ax.text(-0.18,0.25,label_ion,fontsize=14,
					bbox=tagbox,rotation=0,ha="center",va="top",color='k',
					transform=ax.transAxes)
			legends.append(tb_ion)
			ax.axis('off')
			# handle inset mesh
			if mesh_inset:
				axins = inset_axes(ax,width="35%",height="35%",loc=3)
				image = mpimg.imread(os.path.join(work.plotdir,
					#---! somewhat clumsy renaming scheme to get the mesh snapshots
					re.sub('.png$','.v1.png',re.sub('snapshot','snapshot_mesh_review_zoom',fn))))
				# cut the white space and add some padding
				whites = [10,10]
				border_lims_raw = get_blank_border(image)
				sizes = np.array([i.ptp() for i in np.array(border_lims_raw)])
				padding = np.array([(max(sizes)-s)/2.+whites[ss] 
					for ss,s in enumerate(sizes)]).astype(int)
				border_lims = np.array([[b[0]-padding[bb],b[1]+padding[bb]]
					for bb,b in enumerate(border_lims_raw)])
				# swap in the full image after we use the highlights to zoom
				image = mpimg.imread(os.path.join(work.plotdir,
					#! somewhat clumsy renaming scheme to get the mesh snapshots
					re.sub('.png$','.v1.png',re.sub('snapshot','snapshot_mesh_review',fn))))
				image_rect = image[slice(*(border_lims[1])),slice(*(border_lims[0]))]
				axins.imshow(image_rect)
				axins.tick_params(axis=u'both',which=u'both',length=0)
				axins.set_xticks([])
				axins.set_yticks([])

	style = 'bars'
	# loop over included salt bridge and hydrogen bond rows
	for knum,key in enumerate(rows_extras):
		kind = [{'hbonds':'hydrogen_bonding','salt':'salt_bridges'}[key]]
		axes_reorder = dict(enumerate(axes_counts[knum]))
		merged = False
		normed = True
		symmetric = True
		spec_this = {
			'merged':merged,'kinds':kind,'normed':normed,'resnames_exclude':resnames_exclude,
			'key_residues':key_residues,'symmetric':symmetric,'legend_position':legend_position,
			'title_style':title_style,'style':style}
		spec_this['legend_position'] = None
		# ensure PtdIns comes first in this instance
		spec_this['donor_acceptor'] = False
		name_this = (key+{'violin_special':'','bars':''}[style]+{1:'.normed',0:''}[normed]+
			{1:'.merged',0:''}[merged]+{1:'.symmetric',0:''}[symmetric])
		sb.set_style("white")
		axes_mapping = {('PtdIns','PtdIns'):0,('DOPS','PtdIns'):1,('DOPE','PtdIns'):2,('CHL1','PtdIns'):3,}
		if 0: sns_custom = ['membrane-v536','membrane-v538','membrane-v531',
			'membrane-v533','membrane-v599','membrane-v532','membrane-v534']
		sns_custom = ['membrane-v538','membrane-v531',
			'membrane-v533','membrane-v599','membrane-v532','membrane-v534']
		#! global pb
		pb = PlotBonds(dims=(1,0),name=name_this,dataspec=spec_this,
			y_axis_align=False,y_axis_label=True,sns_custom=sns_custom,
			custom_axes=dict(fig=fig,axes_reorder=axes_reorder,axes_mapping=axes_mapping))
		# legend first, since it runs long
		if knum==0: 
			ax_legend = axes_reorder[3]
			legends.append(pb.make_legend(ax=ax_legend,ncol=1,fancy=pb.style=='bars'))

	if False:
		key = 'hbonds'
		merged = False
		kind = ['hydrogen_bonding']
		normed = True
		style = 'bars'
		symmetric = True
		spec_this = {
			'merged':merged,'kinds':kind,'normed':normed,'resnames_exclude':resnames_exclude,
			'key_residues':key_residues,'symmetric':symmetric,'legend_position':legend_position,
			'title_style':title_style,'style':style}
		spec_this['legend_position'] = None
		# ensure PtdIns comes first in this instance
		spec_this['donor_acceptor'] = False
		name_this = (key+{'violin_special':'','bars':''}[style]+{1:'.normed',0:''}[normed]+
			{1:'.merged',0:''}[merged]+{1:'.symmetric',0:''}[symmetric])
		sb.set_style("white")
		axes_reorder = dict(enumerate(axes[:4]))
		ax_legend = axes[3]
		axes_mapping = {('PtdIns','PtdIns'):0,('DOPS','PtdIns'):1,('DOPE','PtdIns'):2,('CHL1','PtdIns'):3,}
		sns_custom = ['membrane-v536','membrane-v538','membrane-v531',
			'membrane-v533','membrane-v599','membrane-v532','membrane-v534']
		#! global pb
		pb = PlotBonds(dims=(1,0),name=name_this,dataspec=spec_this,
			y_axis_align=False,y_axis_label=True,sns_custom=sns_custom,
			custom_axes=dict(fig=fig,axes_reorder=axes_reorder,axes_mapping=axes_mapping))
		legends.append(pb.make_legend(ax=ax_legend,ncol=1,fancy=pb.style=='bars'))
	picturesave('fig.bonding_subset.%s'%layout_name,work.plotdir,form='pdf',
		backup=False,version=True,meta={},extras=legends)

plotrun.routine = None
#! plotrun.routine = []

#!!!!!!!!! manuscript_plot() is a thing you could run to get e.g. fig.bonding_subset.layout_2d.v1.svg

#! being refactored below
#! to make the big tile plots (the refactor was for the 4-panel subset), turn this block on and manually adjust the type below. also turn off the 4-tile plot
if __name__=='__main__' and False:

	"""
	Catalog of atoms participating in various bonds.
	"""

	# getting information about the various lipids from a more central source
	if 'lipid_struct_lib' not in globals():
		#! note that it might make more sense to get the ITP files from the simulations but this requires work
		#! the following connection to automacs should be formalized at some point
		import makeface
		try: mod = makeface.import_remote('amx/amx')
		except: raise Exception('please clone a copy of automacs next to omni in `amx`. '
			'you must also run `make setup all` from that directory to get force field files.')
		struct_lib_spot = 'amx/inputs/charmm/lipid-tops/*.itp'
		lipid_struct_lib = {}
		for g in glob.glob(struct_lib_spot):
			this = mod['GMXTopology'](g)
			for key,val in this.molecules.items():
				if key in lipid_struct_lib: raise Exception('refusing to overwrite %s'%key)
				lipid_struct_lib[key] = val
		# save all atom names for later filtering
		lipid_atom_names = dict([(k,[i['atom'] for i in v['atoms']]) for k,v in lipid_struct_lib.items()])

	# rerun this script with alternate salt bridge cutoffs in ptdins.yaml if desired
	data_type = ['salt_bridges','hydrogen_bonding'][1]
	post_key = 'post_%s'%data_type
	# collect atom counts for normalization
	nmol_counts = dict([(sn,dict(zip(*[data['hydrogen_bonding'][sn]['data'][k] 
		for k in ['resnames','nmols']]))) for sn in sns])
	#! note that the symmetrize function is only used for some basic review during the calculation
	def symmetrize(pairs,vals):
		"""Collapse pairs and values to set-pairs and values."""
		uni = np.unique(np.sort(pairs,axis=1),axis=0)
		toc = [(uu,[
			pp for pp,p in enumerate(pairs)
			if np.any([np.all(u[::i]==p) for i in [1,-1]])
			]) for uu,u in enumerate(uni)]
		reval = dict([(tuple(uni[tt]),
			np.array(vals)[np.array(t)].sum()) for tt,t in toc])
		return reval
	# note that after completing the plot I had to trace the units back to match lipid_lipid_bonds.kind
	#!!! THIS TASK IS STILL IN PROGRESS
	# compute post once per bond designation i.e. salt (and cutoff) vs hydrogen bonds
	if post_key not in globals():
		post,allscores = {},{}
		globals()[post_key] = post
		for snum,sn in enumerate(work.sns()):
			status('processing %s'%sn,i=snum,looplen=len(work.sns()))
			this = data[data_type][sn]['data']
			bonds,obs = [this[i] for i in ['bonds','observations']]
			# prevalence is the mean number of observed bonds of a particular type, fine-grained here
			prev = obs.mean(axis=0)
			ranks = np.argsort(prev)[::-1]
			# check for stupid bonds, which need the most attention
			inds_no_same = np.where(~np.all(bonds[:,:3]==bonds[:,3:6],axis=1))[0]
			inds_no_intra = np.where(~np.all(bonds[:,0:2]==bonds[:,3:5],axis=1))[0]
			if any([i for i in inds_no_intra if i not in inds_no_same]): 
				raise Exception('intra is not inside same!')
			# important to note that the weird self-same bonds (where one atom is bonded to itself, 
			#   (obviously a filter error) are filtered out when filter for intermolecular bonds
			ranked_inter = np.argsort(prev[inds_no_intra])[::-1]
			# subset focus for bonds. note that it is irrelevant that we ranked the bonds here
			prev_sub = prev[inds_no_intra][ranked_inter]
			bonds_sub = bonds[inds_no_intra][ranked_inter]
			# compute a prevalence denominator: the total number of bonds per frame
			denom = prev_sub.sum()
			# compute first of two filters: residue combinations
			# get combinations, ordered
			combos_ord = np.unique(bonds_sub[:,np.array([0,3])],axis=0)
			scores_pair = {}
			symmetric_residue_pairs = np.unique(np.sort(combos_ord,axis=1),axis=0)
			for subject in symmetric_residue_pairs:
				these = np.where(np.all(bonds_sub[:,np.array([0,3])]==subject,axis=1))[0]
				# note some of these scores were just for comparing entire residues
				scores_pair[tuple(subject)] = np.round(prev_sub[these].sum()/nmol_counts[sn][subject[0]],3)
				scores_pair_sym = symmetrize(*zip(*scores_pair.items()))
				#! for i,j in sorted(scores_pair_sym.items(),key=lambda x:x[1])[::-1]: print((i,j))
				# now that we have the rankings pick one type of bond and drill down
				# the following step handles symmetry in the pairs of residues
				# note that above in the for loop we have symmetrized the combos
				subsel = np.any([np.all(bonds_sub[:,np.array([0,3])]==np.array(subject[::i]),axis=1) 
					for i in [1,-1]],axis=0)
				atom_pairs = bonds_sub[subsel][:,np.array([2,5])]
				reverse = bonds_sub[subsel][:,0]!=subject[1]
				# since the salt bridges are symmetric, we flip the atom order
				atom_pairs[np.where(reverse)] = atom_pairs[np.where(reverse)][:,::-1]
				#! rank the pairs atom_pairs[] #!sub_rank = np.argsort(prev_sub[subsel])[::-1]
				if len(atom_pairs)==0: continue
				elif sn not in post: post[sn] = {}
				pairs_red,remap = np.unique(atom_pairs,axis=0,return_inverse=True)
				scores = np.zeros(len(pairs_red))
				for rr,r in enumerate(remap): scores[r] += prev_sub[subsel][rr]
				rerank = np.argsort(scores)[::-1]
				# we normalize by the *first* item in the residue which will become the X-axis in our plot
				# NO! we have to double normalize so the bars have the same meaning across all tiles
				norm = float(nmol_counts[sn][subject[0]])*nmol_counts[sn][subject[1]]
				scores_final = scores[rerank]/norm
				pairs_final = pairs_red[rerank]
				post[sn][tuple(subject)] = dict(pairs=pairs_final,scores=scores_final)
				# debugging and tracing a single example calculation here
				#   goal is to check the prevalence of DOPE-PtdIns bonds in v532
				if False and sn == 'membrane-v534' and all([i in subject for i in ['DOPE','P35P']]):
					print('DEBUG!')
					# problem is that prev_sub[subsel].sum()/nmol_counts[sn][subject[0]] = 0.217
					#   but this should be about 1, the average PE-PIP2 bonds for v532 in figure 6
					# trying to get this directly
					#! see: obs[:,np.where(np.any([np.all([bonds[:,0]==i,bonds[:,3]==j],axis=0) for i,j in [['DOPE','PI2P'],['PI2P','DOPE']]],axis=0))[0]].sum(axis=1).mean() is 43.413
					#  which is identical to prev_sub[subsel].sum()
					# so perhaps the hydrogen bond plot is wrong? time to debut that one
					import ipdb;ipdb.set_trace()
			allscores[sn] = dict(sym=scores_pair_sym,reg=scores_pair)
	"""
	Organizing the visualization.
	1. exclude residues on the outer leaflet
	2. rows are simulations and columns are pairs hence we need a canonical list of pairs
		with a mapping from the isoforms to PtdIns"
	3. develop a scheme to make the atom lists coherent across pairs so we use repetition to simplify
	"""
	if True:
		post = globals()[post_key]
		col_names = sorted(work.sns())
		ptdins_namer = lambda x: x if x not in ['SAPI','P35P','PI2P'] else 'PtdIns'
		exclude_resnames = ['POPC']
		pair_names = list(set([tuple([ptdins_namer(k) for k in i]) for j in 
			[post[sn].keys() for sn in post] for i in j if not any([m in exclude_resnames for m in i])]))
		row_names = sorted(pair_names)
		"""
		outstanding issues:
			canon_names['CHL1'] is missing OP3
			since this appears on the x axis when chl1 is on the y, this means axes are backwards
			perhaps there are no OP5X bonds for 3,5 and magnesium, but even so, this causes the 5 rows/cols to be blanked out even though they are in-play. 3 is fine.
				recommend increasing to 4.6 to see if this is the problem
		"""
		# now that we have columns and rows for the calculation, we create a canonical ordering for atoms
		#! the following canon scheme was deprecated because faulty
		if False:
			canon_names = {}
			for sn in post:
				for pair in post[sn]:
					for ii,i in enumerate(pair):
						resname = ptdins_namer(i)
						if resname in exclude_resnames: continue
						if resname not in canon_names: canon_names[resname] = []
						canon_names[resname].extend(post[sn][pair]['pairs'][:,ii])
			for k,v in canon_names.items(): canon_names[k] = sorted(list(set(v)))
			#! canon names by residue shows e.g. cholesterol versus 3,5 has no spot for OP33 which is weird
			#!   hence we just use one canon_names here
			canon_names_explicit = canon_names
			canon_names_all = sorted(list(set([i for j in canon_names.values() for i in j])))
			#! canon_names = dict([(i,canon_names_all) for i,j in canon_names_explicit.items()])
		"""
		# more debugging when we realized OP3X was missing! this identified a bug in salt_bridges.py
		sn = 'membrane-v533'
		this = data[data_type][sn]['data']
		bonds,obs = [this[i] for i in ['bonds','observations']]
		bonds[np.where(np.any([bonds[:,i]=='P35P' for i in [0,3]],axis=0))[0]][:,np.array([2,5])]
		np.unique(bonds[np.where(np.any([bonds[:,i]=='P35P' for i in [0,3]],axis=0))[0]][:,np.array([2,5])])
		# only OP52, OP53, OP54
		"""
		"""
		debugging the second problem, matching hydrogen bond peaks to 
		pb.raw['membrane-v536'][('DOPE','PtdIns')].shape
			is 801 long and averages 25.52, the average number of bonds
		for v532
		pb.raw[sn][('DOPE','PtdIns')].mean()
			43.41323345817728
		however here we find that there are far less in the scores
		post[sn][('DOPE','PI2P')]['scores'].sum()
			0.2170661672908863
		note that the difference is a factor of exactly 200 for some reason?
		>>> pb.raw[sn][('PtdIns','PtdIns')].mean()
		2.017478152309613
		>>> post[sn][('PI2P','PI2P')]['scores'].sum()
		0.050436953807740326
		>>> 2.017478152309613/0.050436953807740326
		40.0
		problem appears to be normalization. we might be normalizing twice or by the wrong number
		note that the pb.raw output is not normalized
		!!! note also that this is 0.05 which might be correct but this is the horizontal bars on fig 6
		aside from hatch issues, the numbers are correct for PIP2-PIP2 bonds
		"""
		"""
		final notes on ordering
			hydrogen bond is directional which means order matters for donor/acceptor
			if we are comparing hbonds between two residues with the same name, there will be no
				reversing during the symmetrize step obviously, which means one of the axes
				is the donor and one is the acceptor and you can probably infer it
			for example you see a vertical streak at OP52 on hbonds for 4,5 
			and this is because this has the hydrogen
		final notes on comparison to the hydrogen bond subset plot
			post['membrane-v538'][('DOPE','PI2P')]['scores'].sum()*200./40. is 1.535718750000002
			post['membrane-v538'][('DOPE','PI2P')]['scores'].sum() is 0.3071437500000004
			while pb.raw['membrane-v538'][('DOPE','PtdIns')].mean() is 61.42875
			61.42875/40. is 1.53571875
			hence the 0.3 is normed to PE while 1.5 is normed to PIP2
			basically we cannot do any norming on the bit plot if we want the red bars to all have the same
				meaning. hence I will do the double-normalization
		note that I should have noticed that the maximum bond was 0.307 which is 200/40. away from the max 
			on the hbonds plot!
		"""
		# aesthetics
		fs_small = 8
		fs_huge = 36
		fs_med = 14
		fs_large = 26
		extras = []
		figsize_factor = 4 #! should be 4
		side_bar_mag = True
		cmap = 'Blues'
		axes,fig = panelplot(figsize=(figsize_factor*len(col_names),figsize_factor*len(row_names)),
			layout={'out':{'grid':[len(row_names),1],'hspace':0.5},
			'ins':[{'grid':[1,len(col_names)],'wspace':0.5} for i in row_names],})
		resname_label = lambda r: work.vars['names']['short'][r]
		# important to take the sum of scores and not the max
		max_score = max([post[sn][p]['scores'].sum() for sn in post for p in post[sn]])

		# naming schemes
		atoms_all = sorted(np.unique(np.concatenate([np.unique(
			data[data_type][sn]['data']['bonds'][:,np.array([2,5])]) for sn in work.sns()])))
		canon_names = {}
		resnames_all = list(set([m for n in [i for sn in work.sns() for i in post[sn]] for m in n]))
		for r in resnames_all:
			canon_names[r] = list(set([i for i in atoms_all if i in lipid_atom_names[r]]))

		for snum,sn in enumerate(work.sns()):
			for pnum,rpair in enumerate(pair_names):
				ax = axes[pnum][snum]
				rpair_explicit = tuple([(i if i!='PtdIns' else 
					work.meta[sn]['ptdins_resname']) for i in rpair])
				# standard atom names ensures that all tiles are square, but leaves extra whitespace
				raw = np.ones([len(atoms_all) for i in rpair]+[4])*\
					mpl.cm.__dict__[cmap](0.0)
				# not that if you use zeros above, you get white where there were never any 
				#   observations. we set the background to the zero for our colormap instead, so that
				#   zero observations is nearly seamless with nearly zero
				names_x,names_y = [canon_names[i] for i in rpair_explicit]
				if rpair_explicit in post[sn]:
					names_x_this,names_y_this = [sorted(list(set(i))) 
						for i in np.transpose(post[sn][rpair_explicit]['pairs'])]
				else: names_x_this,names_y_this = [],[]
				invalid_x,invalid_y = [np.array([ii for ii,i in enumerate(atoms_all) 
					if i not in n]) for n in [names_x,names_y]]
				if rpair_explicit in post[sn]:
					for ii,(a0,a1) in enumerate(post[sn][rpair_explicit]['pairs']):
						this_max = post[sn][rpair_explicit]['scores'].max()
						raw[atoms_all.index(a0),atoms_all.index(a1)] = \
							mpl.cm.__dict__[cmap](post[sn][rpair_explicit]['scores'][ii]/this_max)
				if not side_bar_mag: imshow_kwargs = dict(vmax=max_score,vmin=0)
				else: imshow_kwargs = dict(vmax=1.0,vmin=0)
				# set the invalid color
				if len(invalid_y)>0: raw[invalid_y,:] = (0.,0.,0.,0.15)
				if len(invalid_x)>0: raw[:,invalid_x] = (0.,0.,0.,0.15)
				raw = np.transpose(raw,(1,0,2))
				ax.imshow(raw,interpolation='nearest',origin='lower',**imshow_kwargs)
				ax.set_xticks(np.arange(len(atoms_all)))
				ax.set_xticklabels(atoms_all,fontsize=fs_small-2,rotation=-90)
				ax.set_yticks(np.arange(len(atoms_all)))
				ax.set_yticklabels(atoms_all,fontsize=fs_small-2)
				#! mpl finally changes from 'off' to a proper boolean
				ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
				ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
				#! note that labels were backwards until thankfully I noticed a discrepancy
				#!   because cholesterol literally never made a bond with OP3X on 3,5 hence it was missing
				#!   and ultimately this led me to reformulate canon names
				ax.set_ylabel(resname_label(rpair_explicit[0]),fontsize=fs_med)
				ax.set_xlabel(resname_label(rpair_explicit[1]),labelpad=10,fontsize=fs_med)
				if side_bar_mag:
					axins = inset_axes(ax,width="10%",height="100%",loc=3,
						bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
					axins.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
					axins.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
					axins.set_xlim((0,1))
					axins.set_ylim((0,1))
					axins.set_xticks([])
					axins.set_yticks([])
					if rpair_explicit in post[sn]:
						# important to take the sum and not the max when comparing tiles
						axins.bar([0.5],[post[sn][rpair_explicit]['scores'].sum()/max_score],
							color='red',width=1.0)
					extras.append(axins)
				if pnum==0: 
					ax.set_title('%s\n%s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label']),
						fontsize=fs_huge,y=1.1)
				if snum==0:
					tagbox = dict(facecolor='w',lw=0,alpha=1.0,boxstyle="round,pad=0.35")
					label_pair = '%s\n%s'%tuple([resname_label(r) if r!='PtdIns' else r for r in rpair])
					extras.append(ax.text(-1.0,0.5,label_pair,fontsize=fs_large,
						bbox=tagbox,rotation=0,ha="center",va="center",color='k',
						transform=ax.transAxes))
		if data_type=='salt_bridges':
			tag = '.cutoff%.1f'%calc['salt_bridges']['calcs']['specs']['distance_cutoff']
		else: tag = ''
		picturesave('fig.atomwise.%s'%data_type+tag,work.plotdir,form='pdf',
			backup=False,version=True,meta={},extras=extras)
		#! document the careful matching where you confirmed that these results match the 
		#!   results on the counts. possibly described in a comment above
		print('max score is %s'%max_score)

def prep_atomwise():
	global lipid_struct_lib
	# getting information about the various lipids from a more central source
	if 'lipid_struct_lib' not in globals():
		#! note that it might make more sense to get the ITP files from the simulations but this requires work
		#! the following connection to automacs should be formalized at some point
		import makeface
		try: mod = makeface.import_remote('amx/amx')
		except: raise Exception('please clone a copy of automacs next to omni in `amx`. '
			'you must also run `make setup all` from that directory to get force field files.')
		struct_lib_spot = 'amx/inputs/charmm/lipid-tops/*.itp'
		lipid_struct_lib = {}
		for g in glob.glob(struct_lib_spot):
			this = mod['GMXTopology'](g)
			for key,val in this.molecules.items():
				if key in lipid_struct_lib: raise Exception('refusing to overwrite %s'%key)
				lipid_struct_lib[key] = val
		# save all atom names for later filtering
		lipid_atom_names = dict([(k,[i['atom'] for i in v['atoms']]) for k,v in lipid_struct_lib.items()])
	# rerun this script with alternate salt bridge cutoffs in ptdins.yaml if desired
	data_type = ['salt_bridges','hydrogen_bonding'][1]
	post_key = 'post_%s'%data_type
	# collect atom counts for normalization
	nmol_counts = dict([(sn,dict(zip(*[data['hydrogen_bonding'][sn]['data'][k] 
		for k in ['resnames','nmols']]))) for sn in sns])
	#! note that the symmetrize function is only used for some basic review during the calculation
	def symmetrize(pairs,vals):
		"""Collapse pairs and values to set-pairs and values."""
		uni = np.unique(np.sort(pairs,axis=1),axis=0)
		toc = [(uu,[
			pp for pp,p in enumerate(pairs)
			if np.any([np.all(u[::i]==p) for i in [1,-1]])
			]) for uu,u in enumerate(uni)]
		reval = dict([(tuple(uni[tt]),
			np.array(vals)[np.array(t)].sum()) for tt,t in toc])
		return reval
	# note that after completing the plot I had to trace the units back to match lipid_lipid_bonds.kind
	#!!! THIS TASK IS STILL IN PROGRESS
	# compute post once per bond designation i.e. salt (and cutoff) vs hydrogen bonds
	if post_key not in globals():
		global allscores
		post,allscores = {},{}
		globals()[post_key] = post
		for snum,sn in enumerate(work.sns()):
			status('processing %s'%sn,i=snum,looplen=len(work.sns()))
			this = data[data_type][sn]['data']
			bonds,obs = [this[i] for i in ['bonds','observations']]
			# prevalence is the mean number of observed bonds of a particular type, fine-grained here
			prev = obs.mean(axis=0)
			ranks = np.argsort(prev)[::-1]
			# check for stupid bonds, which need the most attention
			inds_no_same = np.where(~np.all(bonds[:,:3]==bonds[:,3:6],axis=1))[0]
			inds_no_intra = np.where(~np.all(bonds[:,0:2]==bonds[:,3:5],axis=1))[0]
			if any([i for i in inds_no_intra if i not in inds_no_same]): 
				raise Exception('intra is not inside same!')
			# important to note that the weird self-same bonds (where one atom is bonded to itself, 
			#   (obviously a filter error) are filtered out when filter for intermolecular bonds
			ranked_inter = np.argsort(prev[inds_no_intra])[::-1]
			# subset focus for bonds. note that it is irrelevant that we ranked the bonds here
			prev_sub = prev[inds_no_intra][ranked_inter]
			bonds_sub = bonds[inds_no_intra][ranked_inter]
			# compute a prevalence denominator: the total number of bonds per frame
			denom = prev_sub.sum()
			# compute first of two filters: residue combinations
			# get combinations, ordered
			combos_ord = np.unique(bonds_sub[:,np.array([0,3])],axis=0)
			scores_pair = {}
			symmetric_residue_pairs = np.unique(np.sort(combos_ord,axis=1),axis=0)
			for subject in symmetric_residue_pairs:
				these = np.where(np.all(bonds_sub[:,np.array([0,3])]==subject,axis=1))[0]
				# note some of these scores were just for comparing entire residues
				scores_pair[tuple(subject)] = np.round(prev_sub[these].sum()/nmol_counts[sn][subject[0]],3)
				scores_pair_sym = symmetrize(*zip(*scores_pair.items()))
				#! for i,j in sorted(scores_pair_sym.items(),key=lambda x:x[1])[::-1]: print((i,j))
				# now that we have the rankings pick one type of bond and drill down
				# the following step handles symmetry in the pairs of residues
				# note that above in the for loop we have symmetrized the combos
				subsel = np.any([np.all(bonds_sub[:,np.array([0,3])]==np.array(subject[::i]),axis=1) 
					for i in [1,-1]],axis=0)
				atom_pairs = bonds_sub[subsel][:,np.array([2,5])]
				reverse = bonds_sub[subsel][:,0]!=subject[1]
				# since the salt bridges are symmetric, we flip the atom order
				atom_pairs[np.where(reverse)] = atom_pairs[np.where(reverse)][:,::-1]
				#! rank the pairs atom_pairs[] #!sub_rank = np.argsort(prev_sub[subsel])[::-1]
				if len(atom_pairs)==0: continue
				elif sn not in post: post[sn] = {}
				pairs_red,remap = np.unique(atom_pairs,axis=0,return_inverse=True)
				scores = np.zeros(len(pairs_red))
				for rr,r in enumerate(remap): scores[r] += prev_sub[subsel][rr]
				rerank = np.argsort(scores)[::-1]
				# we normalize by the *first* item in the residue which will become the X-axis in our plot
				# NO! we have to double normalize so the bars have the same meaning across all tiles
				norm = float(nmol_counts[sn][subject[0]])*nmol_counts[sn][subject[1]]
				scores_final = scores[rerank]/norm
				pairs_final = pairs_red[rerank]
				post[sn][tuple(subject)] = dict(pairs=pairs_final,scores=scores_final)
				# debugging and tracing a single example calculation here
				#   goal is to check the prevalence of DOPE-PtdIns bonds in v532
				if False and sn == 'membrane-v534' and all([i in subject for i in ['DOPE','P35P']]):
					print('DEBUG!')
					# problem is that prev_sub[subsel].sum()/nmol_counts[sn][subject[0]] = 0.217
					#   but this should be about 1, the average PE-PIP2 bonds for v532 in figure 6
					# trying to get this directly
					#! see: obs[:,np.where(np.any([np.all([bonds[:,0]==i,bonds[:,3]==j],axis=0) for i,j in [['DOPE','PI2P'],['PI2P','DOPE']]],axis=0))[0]].sum(axis=1).mean() is 43.413
					#  which is identical to prev_sub[subsel].sum()
					# so perhaps the hydrogen bond plot is wrong? time to debut that one
					import ipdb;ipdb.set_trace()
			allscores[sn] = dict(sym=scores_pair_sym,reg=scores_pair)
	global post
	post = globals()[post_key]

	# naming schemes
	global atoms_all,canon_names
	atoms_all = sorted(np.unique(np.concatenate([np.unique(
		data[data_type][sn]['data']['bonds'][:,np.array([2,5])]) for sn in work.sns()])))
	canon_names = {}
	resnames_all = list(set([m for n in [i for sn in work.sns() for i in post[sn]] for m in n]))

	lipid_atom_names = dict([(k,[i['atom'] for i in v['atoms']]) 
		for k,v in lipid_struct_lib.items()])
	for r in resnames_all:
		canon_names[r] = list(set([i for i in atoms_all if i in lipid_atom_names[r]]))

class SqueezedNorm(mpl.colors.Normalize):
	# via https://stackoverflow.com/questions/44432693
	def __init__(self, vmin=None, vmax=None, mid=0, s1=2, s2=2, clip=False):
		self.vmin = vmin # minimum value
		self.mid  = mid  # middle value
		self.vmax = vmax # maximum value
		self.s1=s1; self.s2=s2
		f = lambda x, zero,vmax,s: np.abs((x-zero)/(vmax-zero))**(1./s)*0.5
		self.g = lambda x, zero,vmin,vmax, s1,s2: f(x,zero,vmax,s1)*(x>=zero) - \
			f(x,zero,vmin,s2)*(x<zero)+0.5
		mpl.colors.Normalize.__init__(self, vmin, vmax, clip)
	def __call__(self, value, clip=None):
		r = self.g(value, self.mid,self.vmin,self.vmax, self.s1,self.s2)
		return np.ma.masked_array(r)

def plot_atomwise(ax,sn,rpair,fig,do_title=False,
	do_left_label=False,max_score=1.0,alt_bar=False,cbar_width="10%",
	residue_labels=True,alt_bar_bpl=True,annotated_spines=False,
	filter_names=None,annotated_atoms=False):

	from matplotlib.offsetbox import AnchoredText

	# aesthetics
	fs_small = 8
	fs_huge = 36
	fs_med = 14
	fs_large = 26
	extras = []
	figsize_factor = 4 #! should be 4
	side_bar_mag = True
	cmap = 'Blues'
	atoms_all_this = [i for i in atoms_all if i not in ([] or filter_names)]

	resname_label = lambda r: work.vars['names']['short'][r]
	#! ax = axes[pnum][snum]
	rpair_explicit = tuple([(i if i!='PtdIns' else 
		work.meta[sn]['ptdins_resname']) for i in rpair])
	# standard atom names ensures that all tiles are square, but leaves extra whitespace
	raw = np.ones([len(atoms_all_this) for i in rpair]+[4])*\
		mpl.cm.__dict__[cmap](0.0)
	# not that if you use zeros above, you get white where there were never any 
	#   observations. we set the background to the zero for our colormap instead, so that
	#   zero observations is nearly seamless with nearly zero
	names_x,names_y = [canon_names[i] for i in rpair_explicit]
	names_x,names_y = [[i for i in j if i not in 
		(filter_names if filter_names else [])] for j in [names_x,names_y]]
	if rpair_explicit in post[sn]:
		names_x_this,names_y_this = [sorted(list(set(i))) 
			for i in np.transpose(post[sn][rpair_explicit]['pairs'])]
	else: names_x_this,names_y_this = [],[]
	invalid_x,invalid_y = [np.array([ii for ii,i in enumerate(atoms_all_this) 
		if i not in n]) for n in [names_x,names_y]]
	if rpair_explicit in post[sn]:
		this_max = post[sn][rpair_explicit]['scores'].max()
		for ii,(a0,a1) in enumerate(post[sn][rpair_explicit]['pairs']):
			raw[atoms_all_this.index(a0),atoms_all_this.index(a1)] = \
				mpl.cm.__dict__[cmap](post[sn][rpair_explicit]['scores'][ii]/this_max)
	if not side_bar_mag: imshow_kwargs = dict(vmax=max_score,vmin=0)
	else: imshow_kwargs = dict(vmax=1.0,vmin=0)
	this_sum = float(post[sn][rpair_explicit]['scores'].sum())
	#if this_sum<max_score:
	#	imshow_kwargs['norm'] = SqueezedNorm(vmin=0.,vmax=1.,mid=this_sum/max_score,s1=1.7,s2=4)
	# set the invalid color
	if len(invalid_y)>0: raw[invalid_y,:] = (0.,0.,0.,0.15)
	if len(invalid_x)>0: raw[:,invalid_x] = (0.,0.,0.,0.15)
	raw = np.transpose(raw,(1,0,2))
	#imshow_kwargs['extent'] = [0.,len(atoms_all_this)-1.,0.,len(atoms_all_this)-1.]
	im = ax.imshow(raw,interpolation='nearest',origin='lower',
		**imshow_kwargs)
	#! highlight parts of the name
	atoms_all_this_labels = [dict(
		[('OP%d%d'%(n,m),r'OP$\mathbf{%d}%d$'%(n,m)) 
			for m in [2,3,4]
			for n in [3,4,5]]+
		[('O1%d'%(k),r'OP$\mathbf{1}%d$'%(k)) 
			for k in [1,2,3,4]]
		).get(i,i) for i in atoms_all_this]
	ax.set_xticks(np.arange(len(atoms_all_this)))
	ax.set_xticklabels(atoms_all_this_labels,fontsize=fs_small,rotation=-90)
	ax.set_yticks(np.arange(len(atoms_all_this)))
	ax.set_yticklabels(atoms_all_this_labels,fontsize=fs_small)
	#! mpl finally changes from 'off' to a proper boolean
	ax.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
	ax.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
	#! note that labels were backwards until thankfully I noticed a discrepancy
	#!   because cholesterol literally never made a bond with OP3X on 3,5 hence it was missing
	#!   and ultimately this led me to reformulate canon names
	if residue_labels:
		ax.set_ylabel(resname_label(rpair_explicit[0]),fontsize=fs_med)
		ax.set_xlabel(resname_label(rpair_explicit[1]),labelpad=10,fontsize=fs_med)

	if annotated_spines:
		ax2 = ax.twiny()
		ax2.spines["bottom"].set_position(("axes", -0.2))
		ax2.xaxis.set_ticks_position("bottom")
		ax2.xaxis.set_label_position("bottom")
		# Offset the twin axis below the host
		ax2.spines["bottom"].set_position(("axes", -0.2))
		# Turn on the frame for the twin axis, but then hide all 
		# but the bottom spine
		ax2.set_frame_on(True)
		ax2.patch.set_visible(False)
		for sp in ax2.spines.itervalues(): sp.set_visible(False)
		ax2.spines["bottom"].set_visible(True)
		#ax2.set_xlim((0,len(atoms_all_this)))

		ax3 = ax.twinx()
		ax3.spines["bottom"].set_position(("axes", -0.2))
		ax3.yaxis.set_ticks_position("left")
		ax3.yaxis.set_label_position("left")
		# Offset the twin axis below the host
		ax3.spines["bottom"].set_position(("axes", -0.2))
		# Turn on the frame for the twin axis, but then hide all 
		# but the bottom spine
		ax3.set_frame_on(True)
		ax3.patch.set_visible(False)
		for sp in ax3.spines.itervalues(): sp.set_visible(False)
		ax3.spines["left"].set_visible(True)
		#ax3.set_ylim((0,len(atoms_all_this)))

		#! getting adjustable='datalim' is not allowed when both axes are shared.
		#!   so probably cannot do the extra spine thing on two axes

	#! abandoned, second attempt
	#! see https://matplotlib.org/users/annotations.html
	elif annotated_atoms:

		"""
		this very nearly works. you get a bracket for two points on the x-axis
		but you cannot move the data off the axis or the bracket disappears
		see https://stackoverflow.com/questions/18537879/matplotlib-how-to-write-annotation-outside-the-drawing-in-data-coords
		possibly see https://matplotlib.org/users/transforms_tutorial.html#blended-transformations
		"""

		#ax.set_xlim((0,len(atoms_all_this)))
		#ax.set_ylim((0,len(atoms_all_this)))

		#fac = float(len(atoms_all_this))
		#fac = 10.0
		#x1, y1 = 0./fac,0./fac
		#x2, y2 = 0./fac,3.0/fac
		# if you go a little negative then it disappears
		# if you make the 
		x1, y1 = -0.0,0.0
		x2, y2 = -0.0,3.0
		connectionstyle = "bar,fraction=1.0"
		ax.plot([x1, x2],[y1, y2], ".", clip_on=False,)
		#! do not use text here or it screws it up a bit
		extras.append(ax.annotate("",
			xy=(x1, y1), xycoords='data',
			xytext=(x2, y2), textcoords='data',
			arrowprops=dict(arrowstyle="-", #linestyle="dashed",
				color="0.6",
				shrinkA=10, shrinkB=10,
				patchA=None,
				patchB=None,
				connectionstyle=connectionstyle,
				),
			))
		#ax.add_artist(AnchoredText('hi',loc=2,prop=dict(size=8)))
		#ax.set_yticks([])
		#ax.set_yticklabels([])

	if side_bar_mag:
		axins = inset_axes(ax,width=cbar_width,height="100%",loc=3,
			bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
		axins.tick_params(axis='y',which='both',left=False,right=False,labelleft=True)
		axins.tick_params(axis='x',which='both',top=False,bottom=False,labelbottom=True)
		axins.set_xlim((0,1))
		axins.set_ylim((0,1))
		axins.set_xticks([])
		axins.set_yticks([])
		if rpair_explicit in post[sn]:
			if not alt_bar:
				# important to take the sum and not the max when comparing tiles
				axins.bar([0.5],[post[sn][rpair_explicit]['scores'].sum()/max_score],
					color='red',width=1.0)
			else:
				if this_sum<max_score:
					clist = [(0.,mpl.cm.__dict__['Blues'](0.0)),
						(this_sum/max_score,mpl.cm.__dict__['Blues'](1.0)),
						(this_sum/max_score+10**-5,"white"),(1,"white")]
				else: clist = [(0.,mpl.cm.__dict__['Blues'](0.0)),(1.,mpl.cm.__dict__['Blues'](1.0))]
				dummy_im = ax.imshow(raw,
					cmap=mpl.colors.LinearSegmentedColormap.from_list("name",clist),
					interpolation='nearest',origin='lower')
				dummy_im.set_visible(False)
				cbar = plt.colorbar(dummy_im,cax=axins)
				# we can score bonds per lipid in some cases if the counts are the same for all tiles
				#   hence we rely on the calling function to decide this
				if alt_bar_bpl:
					these_yticks = cbar.ax.get_yticks()
					# note that the max score corresponds to the normalized score: the number of 
					#   bonds divided by the nmol of both lipids hence we multiply by one here
					if rpair_explicit[0]!=rpair_explicit[1]: raise Exception
					cbar.ax.set_yticklabels(['%.1f'%(max_score*t*nmol_counts[sn][rpair_explicit[0]]) 
						for t in these_yticks])
					cbar.ax.set_title('bonds\nper\nlipid',fontsize=10)
					cbar.ax.tick_params(axis='y',which='both',left=False,
						right=False,labelleft=False,labelright=True)
				else:
					cbar.ax.tick_params(axis='y',which='both',left=False,
						right=False,labelleft=False,labelright=False)
		extras.append(axins)
	if do_title: # if pnum==0 
		ax.set_title('%s\n%s'%(work.meta[sn]['ion_label'],work.meta[sn]['ptdins_label']),
			fontsize=fs_huge,y=1.1)
	if do_left_label: # if snum==0
		tagbox = dict(facecolor='w',lw=0,alpha=1.0,boxstyle="round,pad=0.35")
		label_pair = '%s\n%s'%tuple([resname_label(r) if r!='PtdIns' else r for r in rpair])
		extras.append(ax.text(-1.0,0.5,label_pair,fontsize=fs_large,
			bbox=tagbox,rotation=0,ha="center",va="center",color='k',
			transform=ax.transAxes))

# pulling together atomwise
if __name__=='__main__' and False:

	plt.rcParams['text.latex.preamble']=[r"\usepackage{xcolor}"]
	
	targets = {
		0:{'sn':'membrane-v599','rpair':('PtdIns','PtdIns')},
		#! note v538 is the same pattern as v530 but the latter just has more
		1:{'sn':'membrane-v538','rpair':('PtdIns','PtdIns')},
		2:{'sn':'membrane-v532','rpair':('PtdIns','PtdIns')},
		3:{'sn':'membrane-v534','rpair':('PtdIns','PtdIns')},}
	filter_names = 'N O13A O13B'.split()
	data_type = ['salt_bridges','hydrogen_bonding'][1]
	prep_atomwise() #! uses globals a lot
	post_key = 'post_%s'%data_type
	post = globals()[post_key]

	extras = []
	if data_type=='salt_bridges':
		tag = '.cutoff%.1f'%calc['salt_bridges']['calcs']['specs']['distance_cutoff']
	else: tag = ''
	#! better way to plot this?
	max_score = max([post[t['sn']][
		tuple([i if i!='PtdIns' else work.meta[t['sn']]['ptdins_resname'] 
			for i in  t['rpair']])]['scores'].sum() for t in targets.values()])

	axes,fig = square_tiles(ntiles=4,figsize=(8,8),wspace=0.4,hspace=0.4,favor_rows=True)
	for anum,target in targets.items():
		ax = axes[anum]
		rpair = target['rpair']
		sn = target['sn']
		plot_atomwise(ax,sn,fig=fig,rpair=rpair,max_score=max_score,alt_bar=True,cbar_width="5%",
			residue_labels=False,alt_bar_bpl=True,filter_names=filter_names)
		rpair_explicit = tuple([(i if i!='PtdIns' else 
			work.meta[sn]['ptdins_label']) for i in rpair])
		if rpair_explicit[0]!=rpair_explicit[1]: raise Exception
		ax.set_title(rpair_explicit[0]+
			' with '+work.meta[target['sn']]['ion_label'])
	picturesave('fig.atomwise.subset.%s'%data_type+tag,work.plotdir,form='svg',
		backup=False,version=True,meta={},extras=extras)
