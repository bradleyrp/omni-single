#!/usr/bin/env python

"""
Successor to plot-actinlink_bonds.py.
"""

import time,copy,collections
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper,dictsum
from base.compute_loop import basic_compute_loop
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import brewer2mpl

# SETTINGS
#! note cannot put these in main because functions called by the loader might need them
residue_types = {'ARG':'basic','HIS':'basic','LYS':'basic',
	'ASP':'acidic','GLU':'acidic','SER':'polar','THR':'polar','ASN':'polar','GLN':'polar',
	'ALA':'hydrophobic','VAL':'hydrophobic','ILE':'hydrophobic','LEU':'hydrophobic',
		'MET':'hydrophobic','PHE':'hydrophobic','TYR':'hydrophobic','TRP':'hydrophobic',
	'CYS':'special','SEC':'special','GLY':'special','PRO':'special'}
residue_codes = {'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E',
	'SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','SEL':'U','GLY':'G','PRO':'P',
	'ALA':'A','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V'}
residue_codes_reverse = dict([(j,i) for i,j in residue_codes.items()])
residue_type_colors = {'basic':'Blues','acidic':'Reds','hydrophobic':'Greens',
	'polar':'Purples','special':'Oranges'}
residue_colors = dict([(name,residue_type_colors[residue_types[name]]) for name in residue_types])
ticks_font = mpl.font_manager.FontProperties(family='Latin Modern Mono',style='normal',
	size=14,weight='normal',stretch='normal')
#! hardcoded replicate mapping for reducing simulations into one
replicate_mapping = [('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
	('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
	('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
	('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])]
# augment the metadata with better labels
extra_labels = {
	'pip2_20_no_chol':r'mDia2, 20% $PIP_2$, no CHOL ($\times2$)',
	'pip2_20':r'mDia2, 20% $PIP_2$ ($\times2$)',
	'pip2_30':r'mDia2, 30% $PIP_2$ ($\times2$)',
	'pip2_10':r'mDia2, 10% $PIP_2$ ($\times2$)',}
extra_labels_short = {
	'pip2_20_no_chol':r'mDia2, 20% $PIP_2$'+'\n'+r'no CHOL ($\times2$)',
	'pip2_20':r'mDia2'+'\n'+r'20% $PIP_2$ ($\times2$)',
	'pip2_30':r'mDia2'+'\n'+r'30% $PIP_2$ ($\times2$)',
	'pip2_10':r'mDia2'+'\n'+r'10% $PIP_2$ ($\times2$)',}
for key,val in extra_labels.items(): 
	work.metadata.meta[key] = dict(label=val)
for key,val in extra_labels_short.items(): 
	work.metadata.meta[key]['label_compact'] = val
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
# define salt bridges
valid_salt_bridges = [
	{'resname':'ARG','atoms':['NH1','NH2']},
	{'resname':'LYS','atoms':['NZ']},
	{'resname':'HIS','atoms':['ND1','NE2']},]
# special subselections
special_protein_parts = {'nwaspbilayernochl':np.arange(0,22+1).astype(int)}
rowspec = ['subject_resname','subject_resid','subject_atom',
	'target_resname','target_resid','target_atom']
# common colors
colors = {'DOPC':'blue','DOPS':'red','POP2':'magenta','PI2P':'magenta',
	'all lipids':'black','DOPE':'blue'}
colors = {'DOPE':'#808080','DOPS':'#000080','PI2P':'#FF0000'}
lipid_label = lambda x: dict([(i,r'$\mathrm{{PIP}_{2}}$') 
	for i in work.vars.get('selectors',{}).get('resnames_PIP2',{})]).get(x,x)
sn_title = lambda sn,which='label': work.meta.get(sn,{}).get(which,re.sub('_','-',sn))
#! hardcoding lipid types here
lipid_types_ref = ['PI2P','DOPS','DOPE']

@autoload(plotrun)
def load():
	sns_contacts,(data_contacts,calc_contacts) = work.sns(),plotload('contacts',work)
	sns_hbonds,(data_hbonds,calc_hbonds) = work.sns(),plotload('hydrogen_bonding',work)
	sns_mdia2_ordering = ['mdia2bilayer_nochl2','mdia2bilayer_nochl3','mdia2bilayer10','mdia2bilayer10_2',
		'mdia2bilayerphys','mdia2bilayerphys2','mdia2bilayer30','mdia2bilayer30_2']
	if sns_hbonds!=sns_contacts: 
		raise Exception('collections for hydrogen_bonding and contacts are not equal')
	else: sns = sns_mdia2_ordering
	colors.update(**dict([(sn,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[sns.index(sn)]) 
		for sn in sns]))

class ActinlinkPlotter:
	"""Supervise the mapping from data to plot."""
	def __init__(self,name,order):
		"""Receive a key value pair from the request list."""
		self.name = name
		order_this = copy.deepcopy(order)
		self.sns = order_this.pop('sns')
		self.plot = order_this.pop('plot',{})
		self.namer = order_this.pop('namer',None)
		self.ps = self.plotspec = order_this.pop('plotspec',{})
		# get data
		self.ds = self.dataspec = order_this.pop('data')
		self.get_data(**self.ds)
		self.tags = order_this.pop('tags',{})
		if order_this: raise Exception('unprocessed arguments %s'%order_this)
	def get_data_bonds(self,kind,combine_replicates=False):
		"""Get the bond data from globals."""
		# filtering steps
		if kind in ['data_salt','data_salt_hbonds'] and 'data_salt' not in globals(): 
			compute_salt_bridges_from_contacts()
		if kind=='data_salt_hbonds' and 'data_salt_hbonds' not in globals(): 
			merge_salt_hbonds()
		self.data_bonds = globals()[kind]
		if kind=='data_contacts':
			# detect and select multiple cutoffs
			if ('contacts',0) in self.data_bonds.keys():
				# get the calc
				keys = [key for key in self.data_bonds if calc_contacts[key].get('calcs',{}).get(
					'specs',{}).get('cutoff',None)==self.contacts_cutoff]
				if len(keys)!=1: 
					raise Exception(
						'cannot find contacts with cutoff %s in data_contacts'%self.contacts_cutoff)
				else: self.data_bonds = data_contacts[keys[0]]
		# combine replicates
		if combine_replicates:
			#! soft code this in the yaml
			self.replicate_mapping = replicate_mapping
			#! assume full coverage in the mapping
			for supkey,sns_sub in self.replicate_mapping:
				offs = [filter_protein_lipid_bonds(self.data_bonds[sn]['data']) for sn in sns_sub]
				reduced = bond_combinator(*offs,combine_replicates=combine_replicates)
				self.data_bonds[supkey] = {'data':
					dict([(key,reduced[key]) for key in ['bonds','observations','observations_std']])}
				if 'valid_frames' in self.data_bonds[sn]['data']:
					self.data_bonds[supkey]['data']['valid_frames'] = \
						self.data_bonds[sn]['data']['valid_frames']
				else: valid_frames = range(len(reduced['observations']))
			self.sns = zip(*self.replicate_mapping)[0]
	def get_function(self,name):
		"""Get compute funcitons from globals."""
		if name not in globals(): raise Exception('cannot find function %s in globals'%name)
		else: return globals()[name]
	def get_data(self,**kwargs):
		"""Send self to a function to get data."""
		self.contacts_cutoff = kwargs.pop('contacts_cutoff',None)
		# get upstream bond data
		self.get_data_bonds(kwargs.pop('data_bonds'),
			combine_replicates=kwargs.pop('combine_replicates',False))
		# run the compute function
		function_name = kwargs.pop('function')
		self.get_function(function_name)(self,**kwargs.get('kwargs',{}))
	def make_output_name(self):
		"""Name the figure."""
		#! intervene here to actually generate a report of what the figure is
		mods = self.namer.get('mods',[])
		return 'fig.%s%s'%(self.namer['base'],'.'+'.'.join(mods) if mods else '')
	def render(self): globals()[self.plot['function']](self,**self.plot.get('kwargs',{}))

def bond_counter(resid,resname_set):
	"""
	Adapted from the explicit section of the contacts.py function called "count_reduced_contact".
	The hydrogen bonding currently uses the same rowspec as the contacts.
	"""
	global bonds,obs,rowspec
	#---filter the observations by the protein residue (subject_resid) and target resname
	#---...providing a result
	which = np.where(np.all((bonds[:,rowspec.index('subject_resid')].astype(int)==resid,
		np.in1d(bonds[:,rowspec.index('target_resname')],resname_set)),axis=0))
	result = obs.T[which].sum(axis=0)
	return result

def compute_reduced_bonds(bonds,obs):
	"""Convert atom-wise bonds into residue-wise bonds."""
	rowspec_this = ['subject_resname','subject_resid','target_resname','target_resid']
	cols_red = np.array([rowspec.index(r) for r in rowspec_this])
	bonds_this = bonds[:,cols_red]
	idx,counts = uniquify(bonds_this)
	# note that you missed the innermost all in the following resulting in confoundingly high multivalency
	obs_this = np.array([obs.T[np.where(np.all(bonds[:,cols_red]==b,axis=1))[0]].any(axis=0) for 
		b in bonds[idx][:,cols_red]]).T
	return bonds_this[idx],obs_this

def compute_contact_maps(self,explicit=False):
	"""Prepare data for contact maps by breaking things down by protein residue."""
	# postdat is indexed by simulation, then bond type
	postdat = dict([(sn,{}) for sn in sns])
	# loop over simulations
	for sn in sns:
		# common features
		lipid_resnames = lipid_types_ref
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		lipid_resnames = list(zip(*resname_combos))[0]
		postdat[sn]['lipid_resnames'] = lipid_resnames
		#! get data from salt. note that we sometimes combine salt and hbonds and they must be identical
		#! ... anyway. admittedly the following is ugly
		try: resids = postdat[sn]['resids'] = self.data_bonds[sn]['data']['subject_residues_resids']
		except: resids = postdat[sn]['resids'] = data_salt[sn]['data']['subject_residues_resids']
		try: protein_resnames = postdat[sn]['subject_residues_resnames'] = \
			self.data_bonds[sn]['data']['subject_residues_resnames']
		except: protein_resnames = postdat[sn]['subject_residues_resnames'] = \
			data_salt[sn]['data']['subject_residues_resnames']
		try: postdat[sn]['times'] = self.data_bonds[sn]['data']['times']
		except: postdat[sn]['times'] = data_contacts.values()[0][sn]['data']['times']
		bonds,obs = [self.data_bonds[sn]['data'][k] for k in ['bonds','observations']]
		nframes = len(obs)
		# ensure bond list is peptide then lipid
		#! note that this might flip lipid-lipid bonds but they will be filtered out later
		if any(np.in1d(bonds[:,rowspec.index('subject_resname')],lipid_types_ref)): 
			bonds = bond_strict_order(bonds)
		# reduce the explicit bonds to a residue-lipid bond thereby ignoring the specific atoms
		if not explicit: 
			bonds_this,obs_this = compute_reduced_bonds(bonds,obs)
			rowspec_this = ['subject_resname','subject_resid','target_resname','target_resid']
		else: bonds_this,obs_this,rowspec_this = bonds,obs,rowspec
		# get bonds for different combinations
		compactor = dict([(combo_name,np.zeros((len(resids),nframes))) 
			for combo_name,combos in resname_combos])
		# loop over lipid combinations
		for combo_name,lipid_resnames in resname_combos:
			# loop over protein resids
			for rnum,resid in enumerate(resids):
				# get lipid-peptide bonds (subject/target or vis versa) with lipids in the resnames
				inds = np.where(np.any([np.all([np.in1d(bonds_this[:,rowspec_this.index(i)],j) 
					for i,j in zip(k,[lipid_resnames,[str(resid)]])],axis=0) for k in 
					[['target_resname','subject_resid'],['subject_resname','target_resid']]],axis=0))[0]
				if False and sn=='mdia2bilayerphys' and rnum==0 and combo_name=='all lipids':
					import ipdb;ipdb.set_trace()
				# explicit counting versus one bond-per-lipid i.e. reduced is handled above
				compactor[combo_name][rnum] = np.sum(obs_this.T[inds],axis=0)
				#if len(inds)>1:
				#	import ipdb;ipdb.set_trace()
		postdat[sn]['compacted'] = compactor
	# package
	self.data = postdat

def compute_salt_bridges_from_contacts():
	"""
	Generate salt bridges from the contacts list.
	"""
	#! we hard-code the definition of a salt bridge here and request the right contacts list
	#! ... move this to the plot specs?
	global cutoff_salt_bridge
	cutoff_salt_bridge = work.metadata.plots['actinlink_bonds_analysis'].get(
		'specs',{}).get('cutoff_salt_bridge',2.2)
	# get the key in data_contacts corresponding to the right cutoff
	data_contacts_cursor = [key for key in calc_contacts 
		if calc_contacts[key].get('calcs',{}).get('specs',{}).get('cutoff',None)==cutoff_salt_bridge]
	if len(data_contacts_cursor)!=1: raise Exception(
		'cannot find the right contacts data for salt with cutoff %s. '%cutoff_salt_bridge+
		'you have to add it to the compute loop and to the loop inside of the plot metadata for this plot.')
	else: data_contacts_cursor = data_contacts_cursor[0]
	data_contacts_this = dict([(sn,data_contacts[data_contacts_cursor][sn]) for sn in sns])
	# share variables with a parallel function
	global bonds,obs,data_salt
	# package the data in the same manner as data_contacts and data_hbonds
	data_salt = dict([(sn,{'data':{}}) for sn in sns])
	for sn in sns:
		#---filter the bonds and observations from contact maps
		bonds_all = data_contacts_this[sn]['data']['bonds']
		obs_all = data_contacts_this[sn]['data']['observations']
		nframes = len(obs_all)
		salt_bridge_inds =[]
		#---loop over frames in the simulation
		for fr in range(nframes):
			status('filtering salt bridges from contact data',i=fr,looplen=nframes,tag='compute')
			#---find observed bonds for that frame
			bonds_inds = np.where(obs_all[fr]==1.0)[0]
			frame = bonds_all[bonds_inds]
			hits_over_salt_bridges = []
			for definition in valid_salt_bridges:
				matches_resname = frame[:,0]==definition['resname']
				matches_atom = np.in1d(frame[:,2],definition['atoms'])
				matches_lipid_oxygen = np.array([i[0] for i in frame[:,5]])=='O'
				matches = np.all((matches_resname,matches_atom,matches_lipid_oxygen),axis=0)
				hits_over_salt_bridges.append(matches)
			frame_matches = np.where(np.any(hits_over_salt_bridges,axis=0))
			#---save the observed salt bridges by index number for the master bond list
			salt_bridge_inds.append(bonds_inds[frame_matches])
		#---get unique indices for the observed salt bridges
		salt_inds = np.unique(np.concatenate(salt_bridge_inds))
		#---set global bonds and obs so they only contain salt bridges and then run the bond_counter
		bonds = bonds_all[salt_inds]
		obs = obs_all[:,salt_inds]
		status('salt nbonds for %s is %d'%(sn,len(salt_inds)),tag='note')
		#---! get resids for the protein and lipid_resnames from contact maps
		lipid_resnames = np.unique(
			data_contacts_this[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resids = data_contacts_this[sn]['data']['subject_residues_resids']
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		#---compute loop
		looper = [{'resid':resid,'resname_set':resname_set} 
			for resid in resids for resname_name,resname_set in resname_combos]
		compute_function = bond_counter
		incoming = basic_compute_loop(compute_function,looper,run_parallel=True)
		data_salt[sn]['data']['bonds'] = bonds
		data_salt[sn]['data']['observations'] = obs
		data_salt[sn]['data']['valid_frames'] = data_contacts_this[sn]['data']['valid_frames']
		data_salt[sn]['data']['times'] = data_contacts_this[sn]['data']['times']
		for key in ['subject_residues_resnames','subject_residues_resids']:
			data_salt[sn]['data'][key] = data_contacts_this[sn]['data'][key]

def color_light(color,percent):
	"""Make RGB colors lighter."""
	color = np.array(mpl.colors.to_rgb(color))
	white = np.array([1.,1.,1.])
	return color+(white-color)*percent

def compute_summary(self,explicit=False):
	"""
	Compute peptide-lipid bonds by lipid type.
	"""
	# unpack
	sns = self.sns
	data_bonds = self.data_bonds
	# compute
	counts,counts_std = {},{}
	for sn in sns:
		get_std = 'observations_std' in data_bonds[sn]['data']
		#! hacking away
		get_std = False
		counts[sn],counts_std[sn] = {},{}
		bonds,obs = [data_bonds[sn]['data'][k] for k in ['bonds','observations']]
		# reduce the explicit bonds to a residue-lipid bond thereby ignoring the specific atoms
		if not explicit: 
			bonds_this,obs_this = compute_reduced_bonds(bonds,obs)
			rowspec_this = ['subject_resname','subject_resid','target_resname','target_resid']
		else: bonds_this,obs_this,rowspec_this = bonds,obs,rowspec
		if get_std: obs_std = data_bonds[sn]['data']['observations_std']
		# loop over lipid types
		for lipid in lipid_types_ref:
			lipid_this = [lipid]
			# filter peptide-lipid bonds
			inds = np.where(np.any((
				np.all((np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],lipid_this),
					np.in1d(bonds_this[:,rowspec_this.index('target_resname')],
						residue_codes.keys())),axis=0),
				np.all((np.in1d(bonds_this[:,rowspec_this.index('target_resname')],lipid_this),
					np.in1d(bonds_this[:,rowspec_this.index('subject_resname')],
						residue_codes.keys())),axis=0),
				),axis=0))
			counts[sn][lipid] = obs_this.T[inds].sum(axis=0)
			#! unclear if this is the right way to do error bars
			if get_std: counts_std[sn][lipid] = np.sqrt((obs_std.T[inds]**2).sum(axis=0))
	# package
	self.data = dict(counts=counts,counts_std=counts_std if get_std else None)

def compute_chargings(self,explicit=False):
	"""
	Compute the lipids bound per residue.
	Acts on instances of ActinlinkPlotter.
	"""
	# unpack
	sns = self.sns
	data_bonds = self.data_bonds
	# compute
	chargings,residue_listings = {},{}
	for sn in sns:
		bonds,obs = [data_bonds[sn]['data'][k] for k in ['bonds','observations']]
		# reduce the explicit bonds to a residue-lipid bond thereby ignoring the specific atoms
		if not explicit: 
			bonds_this,obs_this = compute_reduced_bonds(bonds,obs)
			rowspec_this = ['subject_resname','subject_resid','target_resname','target_resid']
		else: bonds_this,obs_this,rowspec_this = bonds,obs,rowspec
		residues = np.unique(bonds_this[:,rowspec.index('subject_resid')])
		residue_listings[sn] = residues
		lipid_resnames = np.unique(bonds_this[:,rowspec_this.index('target_resname')])
		chargings_this = {}
		for resnum in residues:
			for target in lipid_resnames:
				partners = np.unique(bonds_this[np.where(np.all((
					bonds_this[:,rowspec_this.index('subject_resid')]==resnum,
					bonds_this[:,rowspec_this.index('target_resname')]==target),
					axis=0))[0]][:,rowspec_this.index('target_resid')])
				bond_this = np.zeros(len(obs_this))
				for partner in partners:
					rows = np.where(np.all((
						bonds_this[:,rowspec_this.index('subject_resid')]==resnum,
						bonds_this[:,rowspec_this.index('target_resname')]==target,
						bonds_this[:,rowspec_this.index('target_resid')]==partner,
						),axis=0))[0]
					bond_this += obs_this[:,rows].sum(axis=1)
				chargings_this[(resnum,target)] = bond_this
		chargings[sn] = chargings_this
	# package
	self.data = dict(chargings=chargings,residues=residues)
	return

def subdivide_trajectory(segnum,n_segments,nframes=None):
	"""Evenly subdivide a trajectory."""
	return np.where(segnum==np.floor(np.arange(nframes)/(nframes/float(n_segments))).astype(int))[0]

def bond_combinator(*offs,**kwargs):
	"""Combine bonds and observations between replicates."""
	# note that we discard the final column of hydrogens 
	bonds_all = np.concatenate([o['bonds'][:,:6] for o in offs])
	details = kwargs.pop('combine_replicates',True)
	# all incoming data must be filtered before
	idx,counts = uniquify(bonds_all)
	bonds_master = bonds_all[idx]
	# for histograms and non-dynamic quantities we just concatenate the frames
	if details==True:
		#! only works for combining two
		valid_frames = np.intersect1d(*[o.get('valid_frames',
			range(o['observations'].shape[0])) for o in offs])
		mappings_vf = [np.concatenate([np.where(valid_frames==i)[0] for i in 
			o.get('valid_frames',range(o['observations'].shape[0]))]) for o in offs]
		nframes = len(valid_frames)
		obs = np.zeros((len(offs),nframes,len(bonds_master)))
		for onum,off in enumerate(offs):
			# sum over the incoming lists (instead of concatenating them)
			for bnum,bond in enumerate(off['bonds']):
				status('combining bonds from replicate %d'%onum,i=bnum,looplen=len(off['bonds']))
				inds = np.where(np.all(bonds_master==bond[:6],axis=1))[0]
				if len(inds)!=1: raise Exception
				obs[onum,:,inds[0]] = off['observations'][mappings_vf[onum],bnum]
		# we could compute standard deviation over the two incoming summands but that is qualitatively 
		# ... different than the kind of framewise variation we wish to capture
		#!!! document the various std methods somewhere in this code?
		obs = obs.sum(axis=0)
		obs_std = False
	# alternate merging schemes for dynamics
	elif type(details)==dict:
		# segmentation mode where we divide a trajectory into parts and then take the average between them
		if details.keys()==['n_segments']:
			n_segments = details['n_segments']
			obs = np.zeros((n_segments,len(bonds_master)))
			obs_std = np.zeros((n_segments,len(bonds_master)))
			for onum,off in enumerate(offs):
				nframes = off['observations'].shape[0]
				for bnum,bond in enumerate(off['bonds']):
					status('combining bonds from replicate %d'%onum,i=bnum,looplen=len(off['bonds']))
					inds = np.where(np.all(bonds_master==bond[:6],axis=1))[0]
					#! taking mean here!!!
					obs[:,inds[0]] += np.array([
						off['observations'][subdivide_trajectory(nn,n_segments,nframes),bnum].mean()
						for nn in range(n_segments)])
					obs_std[:,inds[0]] += np.array([
						off['observations'][subdivide_trajectory(nn,n_segments,nframes),bnum].std()
						for nn in range(n_segments)])
		else: raise Exception
	if kwargs.get('debug',False):
		import ipdb;ipdb.set_trace()
	return dict(observations=obs,bonds=bonds_master,observations_std=obs_std)

def bond_strict_order(bonds):
	"""Reverse subject/target to ensure we are considering lipid-peptide and peptide-lipid."""
	#! wherever you see ":6" we are ignoring the hydrogen atom name
	bonds_out = bonds[:,:6]
	reorder_cols = np.array([3,4,5,0,1,2])
	switchers = np.where(np.in1d(bonds[:,rowspec.index('subject_resname')],lipid_types_ref))[0]
	if len(switchers)>0: bonds_out[switchers] = bonds_out[switchers][:,reorder_cols]
	return bonds_out

def filter_protein_lipid_bonds(incoming):
	"""Ensure we are only considering bonds between the proteins and lipids."""
	inds = np.where(np.any((
		np.all((np.in1d(incoming['bonds'][:,rowspec.index('subject_resname')],residue_codes.keys()),
			np.in1d(incoming['bonds'][:,rowspec.index('target_resname')],lipid_types_ref)),axis=0),
		np.all((np.in1d(incoming['bonds'][:,rowspec.index('target_resname')],residue_codes.keys()),
			np.in1d(incoming['bonds'][:,rowspec.index('subject_resname')],lipid_types_ref)),axis=0),
		),axis=0))[0]
	if 'valid_frames' not in incoming:
		valid_frames = np.arange(incoming['observations'].shape[0]).astype(int)
	else: valid_frames = incoming['valid_frames']
	bonds_out = bond_strict_order(incoming['bonds'][inds,:6])
	result = dict(bonds=bonds_out,observations=incoming['observations'].T[inds].T,valid_frames=valid_frames)
	return result

def merge_salt_hbonds():
	"""Combine salt bridges and hydrogen bonds."""
	global data_salt_hbonds
	data_salt_hbonds = {}
	for sn in sns:
		#! apply filters here
		reduced = bond_combinator(
			filter_protein_lipid_bonds(data_salt[sn]['data']),
			filter_protein_lipid_bonds(data_hbonds[sn]['data']),)
		data_salt_hbonds[sn] = {'data':
			dict([(key,reduced[key]) for key in ['bonds','observations']])}
		(data_salt_hbonds[sn]['data']['subject_residues_resnames'],
			data_salt_hbonds[sn]['data']['subject_residues_resids']) = zip(
			*[(i,j) for i,j in zip(data_salt[sn]['data']['subject_residues_resnames'],
			data_salt[sn]['data']['subject_residues_resids']) if i not in lipid_types_ref])

def plot_charging_histograms_by_lipid(sup):
	"""
	Simple histogram view. Superceded by plot_charging_histograms_stacked
	"""
	raise Exception('deprecated. need to add tags to picture save. not recently tested.')
	# unpack
	chargings = sup.data['chargings']
	residues = sup.data['residues']
	target_resname = sup.plot.get('target_resname')
	sup_name = sup.plot.get('sup_name')
	if sup.dataspec['data_bonds']=='data_salt': cutoff_this = cutoff_salt_bridge
	elif sup.dataspec['data_bonds']=='data_contacts': cutoff_this = sup.dataspec['contacts_cutoff']
	else: raise Exception 
	resnames,resids = [[sup.data_bonds[sn]['data'][key] for sn in sns]
		for key in ['subject_residues_resnames','subject_residues_resids']]
	if (any([np.all(resnames[0]!=i) for i in resnames[1:]]) or 
		any([np.all(resids[0]!=i) for i in resids[1:]])):
		raise Exception('inconsistent residues between simulations')
	else: resnames,resids = resnames[0],resids[0]
	resid_to_title = dict([(str(j),'%s%s'%(residue_codes[i],j)) for i,j in zip(resnames,resids)])
	# plot square tiles with overlaid histograms
	# note that the charging data only contains entries for observed bonds hence the resulting
	# ... square tiles plot only has the R, K, Q residues without any intervention
	peak = max([np.concatenate(chargings[sn].values()).max() for sn in sup.sns])
	bins = np.arange(0,peak+1)-0.5
	axes,fig = square_tiles(len(residues),figsize=sup.ps.get('figsize',(8,8)),
		hspace=sup.ps.get('hspace',0.5),wspace=sup.ps.get('wspace',0.5))
	for rnum,res in enumerate(residues):
		ax = axes[rnum]
		for snum,sn in enumerate(sup.sns):
			for key in chargings[sn]:
				if key[0]==res and key[1]==target_resname:
					counts,bins_this = np.histogram(chargings[sn][key],bins=bins,density=True)
					mids = (bins_this[1:]+bins_this[:-1])/2.
					ax.plot(mids,counts,color=color_by_simulation(sn),lw=sup.plotspec.get('lw',2),
						label=sn_title(sn))
		# decoration
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		ax.set_ylabel('probability')
		ax.set_xlabel(sup.ps['sup_name'])
		ax.set_xticks(np.arange(peak))
		ax.set_title(resid_to_title[res])
	legend = axes[-1].legend(loc='upper left',bbox_to_anchor=(1.05,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor(sup.ps.get('legend_edge_color','k'))
	frame.set_facecolor('white')
	suptitle = fig.suptitle('mDia2 %s ($%.1f\AA$) with %s'%(
		sup_name,cutoff_this,work.vars['names']['short'].get(target_resname,target_resname)),fontsize=16)
	meta = dict(tags=sup.tags) if sup.tags else {}
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta=meta,extras=[legend,suptitle])

def get_cutoff(sup):
	if sup.dataspec['data_bonds'] in ['data_salt','data_salt_hbonds','data_hbonds']: 
		cutoff_this = cutoff_salt_bridge
	elif sup.dataspec['data_bonds']=='data_contacts': cutoff_this = sup.dataspec['contacts_cutoff']
	#! get this from the specs for the hydrogen bonds calculation instead
	elif sup.dataspec['data_bonds']=='data_hbonds': cutoff_this = 3.4
	else: raise Exception(sup.dataspec['data_bonds'])
	return cutoff_this

def plot_charging_histograms_stacked(sup,show_zero=False,total=False,
	emphasis=True,emphasis_base=0.35,square_axes=False,invert=False,mode='landscape'):
	"""
	Panels by residue and simulation with stacked histograms for lipid types.
	"""
	total_label = 'total'
	sup_name = sup.plot.get('sup_name')
	cutoff_this = get_cutoff(sup)
	try: resnames = [sup.data_bonds[sn]['data']['subject_residues_resnames'] for sn in sns]
	except: resnames = [data_salt[sn]['data']['subject_residues_resnames'] for sn in sns]
	try: resids = [sup.data_bonds[sn]['data']['subject_residues_resids'] for sn in sns]
	except: resids = [data_salt[sn]['data']['subject_residues_resids'] for sn in sns]
	if (any([np.all(resnames[0]!=i) for i in resnames[1:]]) or 
		any([np.all(resids[0]!=i) for i in resids[1:]])):
		raise Exception('inconsistent residues between simulations')
	else: resnames,resids = resnames[0],resids[0]
	chargings = sup.data['chargings']
	residues = sup.data['residues']
	if total: residues = np.concatenate((residues,['total']))
	axes,fig = panelplot(figsize=sup.ps['figsize'],
		layout={'out':{'hspace':sup.ps.get('hspace',0.0),'wspace':sup.ps.get('wspace',0.0),
		'grid':[1,len(residues)]},
		'ins':[{'hspace':sup.ps.get('hspace',0.0),'wspace':sup.ps.get('wspace',0.0),
		'grid':[len(sup.sns),1]} for i in residues]})
	extras,axes_alt = [],[]
	resid_to_title = dict([(str(j),'%s%s'%(residue_codes[i],j)) for i,j in zip(resnames,resids)])
	peak = max([np.concatenate(chargings[sn].values()).max() for sn in sup.sns])
	bins = np.arange(peak+2)-0.5
	if show_zero: xbin_filt = slice(None,None)
	else: xbin_filt = slice(1,None)
	bond_values = np.arange(peak+1)
	if emphasis:
		alpha_func = lambda i,m,base=emphasis_base: round(base+(1.0-base)*(float(i)/m),3)
	else: alpha_func = lambda i,m,base=emphasis_base: 1.0
	ymax,ymax_total = 0,0
	for rnum,residue in enumerate(residues):
		for snum,sn in enumerate(sup.sns):
			ax = axes[rnum][snum]
			lipid_types = list(set([l for r,l in chargings[sn].keys() if r==residue or residue==total_label]))
			# very rarely there are no lipids that match
			if set(lipid_types)==set(lipid_types_ref): 
				counts_prev = np.zeros(len(bond_values))[xbin_filt]
				for lipid in lipid_types_ref:
					traj = np.concatenate([v for (r,l),v in chargings[sn].items() 
						if l==lipid and (residue==total_label or r==residue)])
					counts,bins_this = np.histogram(traj,bins=bins)
					# normalize
					counts_norm = counts.astype(float)/len(traj)/len(lipid_types_ref)/(
						1.0 if r!=total_label else len(residues)-1.)
					for ii,(x,y) in enumerate(zip(bond_values+0.5,counts_norm)[xbin_filt]):
						color = color_light(colors[lipid],1.-alpha_func(ii,peak))
						if invert:
							ax.barh([x],[y],left=counts_prev[ii],height=1.0,color=color,zorder=2,alpha=1.0,
								label=work.vars['names']['short'][lipid] if ii==len(bond_values)-1 else None)
						else:
							ax.bar([x],[y],bottom=counts_prev[ii],width=1.0,color=color,zorder=2,alpha=1.0,
								label=work.vars['names']['short'][lipid] if ii==len(bond_values)-1 else None)
					counts_prev += counts_norm.astype(float)[xbin_filt]
					ymax = max([ymax,max(counts_prev)])
					if residue==total_label: ymax_total = max([ymax_total,max(counts_prev)])
			# decoration
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			if invert and mode=='landscape' and rnum==len(residues)-1:
				ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
			# note brute force method on the layout because it is more time-efficient. sorry to disappoint
			if invert and mode=='landscape':
				if snum==len(sup.sns)-1: ax.set_xlabel('rate')
				if rnum>0: ax.set_yticklabels([])
				else: 
					ax.set_yticks(bond_values+0.5)
					ax.set_yticklabels(bond_values.astype(int))
					ax.set_ylabel('bonds')
				ax.set_ylim((0 if show_zero else 1,peak+1))
				ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=1,prune='lower'))
				for i in bond_values[1:]: ax.axhline(i,c=sup.ps.get('frame_color','k'),lw=1,zorder=1)
			elif not invert and mode=='landscape':
				if rnum==0: ax.set_ylabel('rate')
				if snum<len(sup.sns)-1: ax.set_xticklabels([])
				else: 
					ax.set_xticks(bond_values+0.5)
					ax.set_xticklabels(bond_values.astype(int))
					ax.set_xlabel('bonds')
				ax.set_xlim((0 if show_zero else 1,peak+1))
				for i in bond_values[1:]: ax.axvline(i,c=sup.ps.get('frame_color','k'),lw=1,zorder=1)
			else: raise Exception
			if mode=='landscape':
				if snum==0: ax.set_title(resid_to_title.get(residue,residue),fontsize=14)
			# deprecated the use of twinx because bugs
			if False and mode=='landscape' and rnum==len(residues)-1:
				ax_alt = ax.twinx()
				ax_alt.set_ylabel(sn_title(sn,which='label_compact'),
					color='k',rotation=0,horizontalalignment='left',fontsize=12)
				plt.setp(ax_alt.spines.values(),color=sup.ps.get('frame_color','k'))
				ax_alt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
				ax_alt.set_yticks([])
				axes_alt.append((ax_alt,rnum,snum))
			# manual labels on the right
			if mode=='landscape' and rnum==len(residues)-1:
				#! would prefer to use padding instead of 1.1 from the origin
				tagbox = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.5")
				tb = ax.text(1.1,0.5,sn_title(sn,which='label_compact'),
					bbox=tagbox,rotation=0,ha="left",va="center",color='k',
					fontsize=14,transform=ax.transAxes)
				extras.append(tb)
			# decoration general
			plt.setp(ax.spines.values(),color=sup.ps.get('frame_color','k'))
		if snum==len(sup.sns)-1 and rnum==0:
			# sometimes the legend is empty if nothing was plotted for a particular bar and panel
			# ... so we make a custom legend here according to the lipid names
			legendspec = []
			for lipid in lipid_types_ref:
				legendspec.append(dict(name=work.vars['names']['short'][lipid],
					patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=colors[lipid])))
			patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
			legend = ax.legend(patches,labels,loc='lower left',
				bbox_to_anchor=sup.ps.get('legend_bbox',(0.0,-0.75,1.,1.)),
				ncol=len(lipid_types_ref),fontsize=14)
			frame = legend.get_frame()
			frame.set_edgecolor(sup.ps.get('legend_edge_color','k'))
			frame.set_facecolor('white')
			extras.append(legend)
	# alternate axes are included in the following decorations (axes_alt deprecated with twinx)
	for ax,rnum,snum,residue in [(axes[rnum][snum],rnum,snum,residue) 
		for rnum,residue in enumerate(residues) for snum,sn in enumerate(sup.sns)]+axes_alt:
		if invert and mode=='landscape': 
			if residue==total_label: ax.set_xlim((0,ymax_total))
			else: ax.set_xlim((0,ymax))
			if snum!=len(sup.sns)-1: ax.set_xticks([]) #ax.set_xticks(np.arange(0,ymax,0.2))
		elif not invert and mode=='landscape':
			if residue==total_label: ax.set_ylim((0,ymax_total))
			else: ax.set_ylim((0,ymax))
			if rnum!=0: ax.set_yticks([])
			#! note that we cannot use square with the twinx used to put labels on the right
			#! ... see https://github.com/matplotlib/matplotlib/issues/7654
			#! note that I gutted the twinx method for putting the label on the right because 
			#! ... it was absolutely impossible to remove tick marks after twinx in the invert layout
			#! ... mode which is an extremely annoying feature and caused me to discard twinx entirely
			if square_axes: ax.set_aspect(peak/ymax)
		else: raise Exception
	meta = dict(tags=sup.tags) if sup.tags else {}
	extras.append(fig.suptitle('%s'%sup_name,fontsize=16))
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta=meta,extras=extras)

def plot_summary(sup,**kwargs):
	"""
	Summarize the peptide-lipid bonds counts.
	"""
	offset = 0.1
	show_error_bars = kwargs.get('show_error_bars',True)
	show_vlines = kwargs.get('show_vlines',False)
	mode = kwargs.get('mode','bars')
	bar_alpha = kwargs.get('bar_alpha',0.6)
	n_segments = kwargs.get('n_segments',1)
	if n_segments==1 and mode!='bars': raise Exception('must use bars for one segment')
	counts,counts_std = [sup.data[key] for key in ['counts','counts_std']]
	extras = []
	axes,fig = square_tiles(1,figsize=sup.ps.get('figsize',(8,8)),
		hspace=sup.ps.get('hspace',0.5),wspace=sup.ps.get('wspace',0.5))
	ax = axes[0]
	mark = 0.0
	for snum,sn in enumerate(sup.sns):
		nframes = list(set([counts[sn][lipid].shape[0] for lipid in lipid_types_ref]))
		if len(nframes)!=1: raise Exception('inconsistent number of frames')
		else: nframes = nframes[0]
		# if we bin in trajectories here, then we use the intra-bin variation as the std
		if nframes!=n_segments:
			traj_raw = np.concatenate(([np.zeros(nframes)],np.array([counts[sn][l] 
				for l in lipid_types_ref]))).cumsum(axis=0)
			traj = np.array([[t[subdivide_trajectory(n,n_segments,nframes=nframes)].mean() 
				for n in range(n_segments)] for t in traj_raw])
			traj_std = np.array([[t[subdivide_trajectory(n,n_segments,nframes=nframes)].std() 
				for n in range(n_segments)] for t in traj_raw])
		else: 
			traj_raw = np.concatenate(([np.zeros(nframes)],np.array([counts[sn][l] 
				for l in lipid_types_ref]))).cumsum(axis=0)
			traj = np.array([[t[subdivide_trajectory(n,n_segments,nframes=nframes)].mean() 
				for n in range(n_segments)] for t in traj_raw])
			traj_std_raw = np.concatenate(([np.zeros(nframes)],np.array([counts_std[sn][l] 
				for l in lipid_types_ref]))).cumsum(axis=0)
			traj_std = np.array([[t[subdivide_trajectory(n,n_segments,nframes=nframes)].mean() 
				for n in range(n_segments)] for t in traj_std_raw])
			#! cannot develop error bars here because we have already combined replicates into n_segments???
		for lnum,lipid in enumerate(lipid_types_ref):
			color = color_light(colors[lipid],1-bar_alpha)
			if mode=='bars':
				# use bars with unity width or they will not be flush!
				xs = np.arange(n_segments)+mark+1.0
				#width_bar = (1.-2*offset)/(n_segments)
				ax.bar(xs,bottom=traj[lnum],height=traj[lnum+1]-traj[lnum],
					width=1.0,color=color,edgecolor=color,lw=0,align='edge',)
				ecolor = color_light('k',0.5)
				if show_error_bars:
					ax.errorbar(xs+0.5,traj[lnum+1],yerr=traj_std[lnum+1],
						zorder=3,fmt='none',ecolor=ecolor,lw=0.5)
			else: raise Exception(mode)
		mark += n_segments+2.
		if snum<len(sup.sns)-1 and show_vlines: ax.axvline(mark,c='k',lw=1)
	ax.set_xlim((0,(2+n_segments)*len(sup.sns)))
	ax.set_ylim((0,ax.get_ylim()[1]))
	ax.set_xticks([(n_segments+2)*(snum+0.5) for snum in range(len(sup.sns))])
	ax.set_xticklabels([sn_title(sn,which='label_compact') for sn in sup.sns],rotation=0)
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	if 'title' in sup.plot: ax.set_title(sup.plot['title'],fontsize=16)
	ax.set_xlabel(r'time $\rightarrow$',fontsize=16)
	ax.set_ylabel(sup.plot['label'],fontsize=16)
	# custom legend
	legendspec = []
	for lipid in lipid_types_ref:
		legendspec.append(dict(name=work.vars['names']['short'][lipid],
			patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=colors[lipid])))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='lower left',
		bbox_to_anchor=sup.ps.get('legend_bbox',(0.0,1.0,1.,1.)),ncol=len(lipid_types_ref),fontsize=14)
	frame = legend.get_frame()
	frame.set_edgecolor(sup.ps.get('legend_edge_color','k'))
	frame.set_facecolor('white')
	extras.append(legend)
	meta = dict(tags=sup.tags) if sup.tags else {}
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta=meta,extras=extras)

def plot_colorstreak_contact_map(sup,bar_style=False,**kwargs):
	"""
	Final resting place for the contact map code. Ported in from plot-actinlink_bonds.py.
	"""
	if len(sup.sns)<3: raise Exception('we need three simulations to plot the title, colorbar, and legend!')
	# custom plotspec from globals
	plotspec = {'fs_xlabel':14,'fs_ylabel':12,'fs_title':20,
		'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,
		'label_color_strength':1.0,'fs_legend_title':20,
		'binary_color_intensity':0.5,'time_tick_interval':20,
		'fs_xticks':11,'fs_yticks':6,'fs_xticks_bars':6}
	plotspec_small = {'fs_legend_title':14,'fs_legend':12,'fs_xticks':12,'fs_yticks':8}
	# other settings
	do_crazy_colors = kwargs.get('do_crazy_colors',True)
	cmap = kwargs.get('cmap','Greys')
	# same residue comparisons
	lipid_resnames = lipid_types_ref[::-1]+['all lipids']
	ceiling = float(max([max([i.max() for i in sup.data[sn]['compacted'].values()]) for sn in sns]))
	ceiling_bar = 0.0
	extras = []
	axes,fig = panelplot(
		layout={'out':{'hspace':0.05,'wspace':0.05,'grid':[1,len(lipid_resnames)]},
		'ins':[{'hspace':0.2,'wspace':0.05,'grid':[len(sns),1]} 
		for i in lipid_resnames]},figsize=sup.ps['figsize'])
	for ss,sn in enumerate(sns):
		for rr,resname in enumerate(lipid_resnames):
			ax = axes[rr][ss]
			# override the protein selections
			resnums = special_protein_parts.get(sn,np.arange(sup.data[sn]['compacted'][resname].shape[0]))
			try: resnames = np.array(sup.data[sn]['subject_residues_resnames'])
			except: resnames = np.array(data_salt[sn]['subject_residues_resnames'])
			if rr==0:
				ax.set_ylabel(sn_title(sn,which='label_compact'),fontsize=plotspec['fs_ylabel'])
				if True: #! hack work.plots[plotname].get('settings',{}).get('show_residue_names',True):
					# never set the ticklabels before the ticks
					ax.set_yticks(np.arange(len(resnames[resnums]))+0.5)
					ax.set_yticklabels([residue_codes[r] for r in resnames[resnums]],
						fontsize=plotspec['fs_yticks'])
					for label in ax.get_yticklabels():
						label.set_color(mpl.cm.__dict__[
							residue_colors.get(residue_codes_reverse[label._text],
							'Greys')](plotspec['label_color_strength'])[:-1])
			else: ax.set_yticks([])
			if ss==0:
				ax.set_title(lipid_label(resname),
					fontsize=plotspec['fs_title'])
			resids = sup.data[sn]['resids']
			if do_crazy_colors:
				scores = np.array([sup.data[sn]['compacted'][resname][rnum].astype(float)/ceiling 
					for rnum in resnums])
				image_data = np.array([mpl.cm.__dict__[residue_colors.get(resnames[rnum],'Greys')](
					sup.data[sn]['compacted'][resname][rnum].astype(float)/ceiling).T[:-1].T 
					for rnum in resnums])
				# note that the crazy color color maps do not go identically to zero
				# ...so we white them out
				image_data[np.where(scores==0.0)] = [1.,1.,1.]
			else: image_data = sup.data[sn]['compacted'][resname][resnums]
			duration = (sup.data[sn]['times'].max()-sup.data[sn]['times'].min())*10**-3
			xtick_interval = plotspec['time_tick_interval'] 
			ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=4,prune='lower'))
			
			if bar_style:
				ax.set_xlabel('bonds',fontsize=plotspec['fs_xlabel'])
				image_data = sup.data[sn]['compacted'][resname][resnums]
				colors = np.array([mpl.cm.__dict__[residue_colors.get(resnames[rnum],'Greys')](0.7)
					for rnum in resnums])
				xs = np.arange(len(resids))
				ys = image_data.mean(axis=1)
				ys_err = image_data.std(axis=1)
				ax.barh(xs,ys,align='edge',height=1.0,color=colors)
				ecolor = 'k' # color_light('k',0.5)
				ax.errorbar(ys,xs+0.5,xerr=ys_err,zorder=3,fmt='none',ecolor=ecolor,lw=0.5)
				ceiling_bar = max([ceiling_bar,max(ys+ys_err)])
			else:
				ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
				im = ax.imshow(image_data,
					origin='lower',interpolation='nearest',extent=[0,duration,0,len(resnums)],
					aspect=float(duration)/len(resnums),
					vmax=ceiling,vmin=0,cmap=mpl.cm.__dict__[cmap])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		# color bar on the second row
		if not bar_style and ss==1 and rr==len(lipid_resnames)-1:
			axins = inset_axes(ax,width="10%",height="100%",loc=3,
				bbox_to_anchor=(1.3,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
			if ceiling>15: cbar_ticks = np.arange(0,ceiling+1,ceiling/10).astype(int)
			else: cbar_ticks = np.arange(ceiling+1).astype(int)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical",ticks=cbar_ticks)
			axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
			if False: cbar.ax.set_ylabel('score',rotation=-90,labelpad=20,fontsize=14)
			#! it would be nice to make the colorbar discrete
		if ss==0 and rr==len(lipid_resnames)-1:
			# title on the first row instead
			tagbox = dict(facecolor='w',lw=2,alpha=1.0,boxstyle="round,pad=0.5")
			tb = ax.text(1.2,0.5,sup.ps['sup_name'],
				bbox=tagbox,rotation=0,ha="left",va="center",color='k',
				fontsize=14,transform=ax.transAxes)
			extras.append(tb)
	if bar_style:
		for ax in [i for j in axes for i in j]:
			ax.set_xlim((0,ceiling_bar))
	# title on the legend
	cutoff = get_cutoff(sup)
	patches,labels = zip(*[(mpl.patches.Rectangle((0,0),1.0,1.0,
		fc=mpl.cm.__dict__[residue_type_colors[kind]](plotspec['legend_color_strength'])),kind) 
		for kind in residue_type_colors])
	legend = axes[-1][1 if bar_style else 2].legend(
		patches,labels,loc='upper left',fontsize=plotspec['fs_legend'],
		bbox_to_anchor=(1.1,0.0,1.,1.),title=None,shadow=False,fancybox=False)
	extras.append(legend)
	frame = legend.get_frame()
	frame.set_edgecolor('white')
	frame.set_facecolor('white')
	if False: plt.setp(legend.get_title(),fontsize=plotspec['fs_legend_title'])
	meta = dict(tags=sup.tags) if sup.tags else {}
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta=meta,extras=extras)

if __name__=='__main__':

	if False:
		# deprecated plots
		orders_deprecated = {
			'charging_simple PIP2 cut_5.0':{
				'plot':{'function':'plot_charging_histograms_by_lipid',
				'target_resname':'PI2P','sup_name':'contacts'},
				'plotspec':{'lw':2,'figsize':(8,8),'legend_edge_color':'w','hspace':0.8},
				'namer':{'base':'charging','mods':['cut_5.0']},
				'sns':sns_mdia2_ordering,'data':{'function':'compute_chargings',
				'data_bonds':'data_contacts','contacts_cutoff':5.0},},
			'charging_simple PIP2 salt':{
				'plot':{'function':'plot_charging_histograms_by_lipid',
				'target_resname':'PI2P','sup_name':'salt bridges'},
				'plotspec':{'lw':3,'figsize':(8,8),'legend_edge_color':'w'},
				'namer':{'base':'charging','mods':['salt','overlay']},
				'sns':sns_mdia2_ordering,'data':{'function':'compute_chargings',
				'data_bonds':'data_salt'},},
			'charging_simple PIP2 salt merged':{
				'plot':{'function':'plot_charging_histograms_by_lipid',
				'target_resname':'PI2P','sup_name':'salt bridges'},
				'plotspec':{'lw':3,'figsize':(8,8),'legend_edge_color':'w'},
				'namer':{'base':'charging','mods':['salt','reps_combo','overlay']},
				'sns':sns_mdia2_ordering,'data':{'function':'compute_chargings',
				'data_bonds':'data_salt','combine_replicates':True},},}

		orders_explicit = {
			# panels of sns by residue with stacked bar histograms
			'charging salt merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','salt','merged']},
				'tags':{'analysis':'charging histograms','merged':True,'bonds':'salt bridges'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'salt bridges',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_salt','combine_replicates':True},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'charging cut_2.2 merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','cut_2.2','merged']},
				'tags':{'analysis':'charging histograms','merged':True,'bonds':'contacts 2.2A'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_contacts',
					'contacts_cutoff':2.2,'combine_replicates':True},
				'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
					'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9'}},
			'charging cut_5.0 merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','cut_5.0','merged']},
				'tags':{'analysis':'charging histograms','merged':True,'bonds':'contacts 5.0A'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				# somewhat slow because of high cutoff hence larger bond list
				'data':{'function':'compute_chargings','data_bonds':'data_contacts',
					'contacts_cutoff':5.0,'combine_replicates':True},
				'plotspec':{'figsize':(16,8),'legend_edge_color':'w',
					'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9'}},
			# panels of sns by residue with stacked bar histograms, not merged
			'charging salt':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','salt']},
				'tags':{'analysis':'charging histograms','merged':False,'bonds':'salt bridges'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'salt bridges',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_salt'},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.25,1.,1.)}},
			'charging cut_2.2':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','cut_2.2']},
				'tags':{'analysis':'charging histograms','merged':False,'bonds':'contacts 2.2A'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_contacts',
					'contacts_cutoff':2.2},
				'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
					'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.5,1.,1.)}},
			'charging cut_5.0':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','cut_5.0']},
				'tags':{'analysis':'charging histograms','merged':False,'bonds':'contacts 5.0A'},
				'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				# somewhat slow because of high cutoff hence larger bond list
				'data':{'function':'compute_chargings','data_bonds':'data_contacts',
					'contacts_cutoff':5.0},
				'plotspec':{'figsize':(16,8),'legend_edge_color':'w',
					'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.5,1.,1.)}},
			'charging salt_hbonds':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','salt_hbonds']},
				'tags':{'analysis':'charging histograms','merged':False,'bonds':'salt bridges + hydrogen bonds'},
				'plot':{'function':'plot_charging_histograms_stacked',
					'sup_name':'hydrogen bonds + salt bridges','sup_name_y':'bonds',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_salt_hbonds'},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.25,1.,1.)}},
			'charging salt_hbonds merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'charging','mods':['stacked','salt_hbonds','merged']},
				'tags':{'analysis':'charging histograms','merged':True,'bonds':'salt bridges + hydrogen bonds'},
				'plot':{'function':'plot_charging_histograms_stacked',
					'sup_name':'hydrogen bonds + salt bridges','sup_name_y':'bonds',
					'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
				'data':{'function':'compute_chargings','data_bonds':'data_salt_hbonds','combine_replicates':True},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9','legend_bbox':(0.0,-0.75,1.,1.)}},
			# peptide-lipid bond totals
			'summary salt':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'salt']},
				'tags':{'analysis':'bonds summary','merged':False,'bonds':'salt bridges'},
				'plot':{'function':'plot_summary',
					'kwargs':{'n_segments':20},'label':'salt bridges'},
				'data':{'function':'compute_summary','data_bonds':'data_salt'},
				'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'summary hbonds':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'hbonds']},
				'tags':{'analysis':'bonds summary','merged':False,'bonds':'hydrogen bonds'},
				'plot':{'function':'plot_summary',
					'kwargs':{'n_segments':20},'label':'hydrogen bonds'},
				'data':{'function':'compute_summary','data_bonds':'data_hbonds'},
				'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'summary salt_hbonds':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'salt_hbonds']},
				'tags':{'analysis':'bonds summary','merged':False,'bonds':'salt bridges + hydrogen bonds'},
				'plot':{'function':'plot_summary',
					'kwargs':{'n_segments':20},'label':'salt bridges + hydrogen bonds'},
				'data':{'function':'compute_summary','data_bonds':'data_salt_hbonds'},
				'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'summary salt merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'salt','merged']},
				'tags':{'analysis':'bonds summary','merged':True,'bonds':'salt bridges'},
				'plot':{'function':'plot_summary','sup_name':'salt bridges',
					'kwargs':{'n_segments':20,'show_error_bars':True},'label':'salt bridges'},
				'data':{'function':'compute_summary','data_bonds':'data_salt',
					'combine_replicates':{'n_segments':20}},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'summary hbonds merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'hbonds','merged']},
				'tags':{'analysis':'bonds summary','merged':True,'bonds':'hydrogen bonds'},
				'plot':{'function':'plot_summary',
					'kwargs':{'n_segments':20,'show_error_bars':True},'label':'hydrogen bonds'},
				'data':{'function':'compute_summary','data_bonds':'data_hbonds',
					'combine_replicates':{'n_segments':20}},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			'summary salt_hbonds merged':{
				'sns':sns_mdia2_ordering,
				'namer':{'base':'summary','mods':['nbins_%d'%20,'salt_hbonds','merged']},
				'tags':{'analysis':'bonds summary','merged':True,'bonds':'salt bridges + hydrogen bonds'},
				'plot':{'function':'plot_summary',
					'kwargs':{'n_segments':20,'show_error_bars':True},'label':'salt bridges + hydrogen bonds'},
				'data':{'function':'compute_summary','data_bonds':'data_salt_hbonds',
					'combine_replicates':{'n_segments':20}},
				'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
					'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
			} # end of the plot requests list. moving to a loop scheme below

	routine = ['contact_maps','charging'][-1:]
	orders = collections.OrderedDict()

	# master types
	bond_types = ['salt','salt_hbonds','hbonds','cut_2.2','cut_5.0']
	bond_types_labels = ['salt bridges','salt bridges + hydrogen bonds',
		'hydrogen bonds',r'contacts $({\vec{r}\leq2.2\AA})$',r'contacts $({\vec{r}\leq5.0\AA})$']
	bond_types_labels_text = ['salt bridges','salt bridges + hydrogen bonds',
		'hydrogen bonds','contacts 2.2A','contacts 5.0A']
	bond_types_labels_newlines = ['salt bridges','salt bridges +\nhydrogen bonds','hydrogen bonds']+[
		'contacts\n($r\leq2.2\AA$)','contacts\n($r\leq5.0\AA$)']
	labelmaker_newlines = dict(zip(bond_types,bond_types_labels_newlines))
	labelmaker_text = dict(zip(bond_types,bond_types_labels_text))
	labelmaker = dict(zip(bond_types,bond_types_labels))

	if 'charging' in routine:
		charging_base = {
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':[]},
			'tags':{'analysis':'charging histograms','merged':True,'bonds':'salt bridges','explicit':True},
			'plot':{'function':'plot_charging_histograms_stacked',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			'data':{'function':'compute_chargings','kwargs':{'explicit':True},
				'data_bonds':'data_salt','combine_replicates':True},
			'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
				'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}}
		for plot_style in ['charging timeseries','charging histograms'][1:]:
			#! 5.0 is too slow with reduced I think?
			for bond_type in ['salt','salt_hbonds','hbonds','cut_2.2','cut_5.0'][:-1]:
			 	for merged in [False,True]:
			 		for explicit in [False,True]:
						key = 'charging %s'%bond_type+\
							(' explicit' if explicit else '')+(' merged' if merged else '')
						orders[key] = copy.deepcopy(charging_base)
						orders[key]['namer']['mods'] = [bond_type]+\
							(['merged'] if merged else [])+(['explicit'] if explicit else [])
						orders[key]['plot']['sup_name'] = labelmaker[bond_type]+\
							(' (explicit)' if explicit else '')
						if re.match('^cut_',bond_type):
							orders[key]['data']['data_bonds'] = 'data_contacts'
							orders[key]['data']['contacts_cutoff'] = float(
								re.match('^cut_(.+)$',bond_type).group(1))
						else: orders[key]['data']['data_bonds'] = 'data_%s'%bond_type
						orders[key]['data']['kwargs']['explicit'] = explicit
						orders[key]['tags']['merged'] = merged
						orders[key]['tags']['explicit'] = explicit
						orders[key]['tags']['bonds'] = labelmaker_text[bond_type]
						orders[key]['tags']['analysis'] = plot_style
						if plot_style=='charging timeseries': 
							orders[key]['namer']['base'] = 'charging_timeseries'
							orders[key]['plot']['function'] = 'plot_summary'
							orders[key]['plot']['kwargs'] = {'n_segments':20}
							orders[key]['plot']['label'] = labelmaker[bond_type]
							orders[key]['data']['function'] = 'compute_summary'
						elif plot_style=='charging histograms':
							orders[key]['namer']['base'] = 'charging_histograms'
						else: raise Exception

	if 'contact_maps' in routine:
		# contact_maps and interacting_residues
		contact_base = {
			'sns':sns_mdia2_ordering,
			'namer':{'base':'contact_maps','mods':[]},
			'tags':{'analysis':'contact maps','merged':False,'bonds':'hydrogen bonds + salt bridges'},
			'plot':{'function':'plot_colorstreak_contact_map',
				'kwargs':{'n_segments':20,'show_error_bars':True,'bar_style':False},
				'label':'hydrogen bonds + salt bridges'},
			'data':{'function':'compute_contact_maps','data_bonds':'data_salt_hbonds',
				'combine_replicates':False},
			'plotspec':{'figsize':(8,16),'legend_edge_color':'w',
				'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}}
		for plot_style in ['contact_map','interacting_residues']:
			for bond_type in ['salt','salt_hbonds','cut_5.0','cut_2.2','hbonds']:
				for merged in [False]: # impossible to combine replicates for this precise dynamics style
					for explicit in [True,False]:
						key = '%s %s'%(plot_style,bond_type)+\
							(' explicit' if explicit else '')+(' merged' if merged else '')
						orders[key] = copy.deepcopy(contact_base)
						# contact_base is set for contact_maps
						if plot_style=='interacting_residues':
							orders[key]['plot']['kwargs']['bar_style'] = True
							orders[key]['namer']['base'] = 'interacting_residues'
							orders[key]['tags']['analysis'] = 'interacting residues'
						orders[key]['namer']['mods'] = [bond_type]+\
							(['explicit'] if explicit else [])+(['merged'] if merged else [])
						if re.match('^cut_',bond_type):
							orders[key]['data']['data_bonds'] = 'data_contacts'
							orders[key]['data']['contacts_cutoff'] = float(
								re.match('^cut_(.+)$',bond_type).group(1))
						else: orders[key]['data']['data_bonds'] = 'data_%s'%bond_type
						orders[key]['data']['kwargs'] = dict(explicit=explicit)
						orders[key]['plotspec']['sup_name'] = labelmaker_newlines[bond_type]+(
							'\n(explicit)' if explicit else '')
						orders[key]['tags']['merged'] = merged
						orders[key]['tags']['bonds'] = labelmaker_text[bond_type]
						orders[key]['tags']['explicit'] = explicit

	#! to develop a single plot, skip the render function and trim the orders list
	#! i = 'contact_map salt' ; orders = {i:orders[i]}
	# loop over requested plots
	for name,spec in orders.items(): 
		orders[name]['plotter'] = ActinlinkPlotter(name,spec)
		orders[name]['plotter'].render()

"""
tasks:
	check flipping on the histograms?
"""
