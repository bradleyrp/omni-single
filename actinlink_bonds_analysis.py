#!/usr/bin/env python

"""
Successor to plot-actinlink_bonds.py.
"""

import time,copy
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
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,
	'label_color_strength':1.0,'fs_legend_title':20,
	'binary_color_intensity':0.5,'time_tick_interval':20,
	'fs_xticks':11,'fs_yticks':6,'fs_xticks_bars':6}
plotspec_small = {'fs_legend_title':14,'fs_legend':12,'fs_xticks':12,'fs_yticks':8}
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
		if order_this: raise Exception('unprocessed arguments %s'%order_this)
	def get_data_bonds(self,kind,combine_replicates=False):
		"""Get the bond data from globals."""
		# filtering steps
		if kind=='data_salt' and 'data_salt' not in globals():
			compute_salt_bridges_from_contacts()
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
				offs = [self.data_bonds[sn]['data'] for sn in sns_sub]
				reduced = bond_combinator(*offs)
				self.data_bonds[supkey] = {'data':
					dict([(key,reduced[key]) for key in ['bonds','observations']])}
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
		self.get_function(function_name)(self)
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
		for key in ['subject_residues_resnames','subject_residues_resids']:
			data_salt[sn]['data'][key] = data_contacts_this[sn]['data'][key]

def compute_chargings(self):
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
		residues = np.unique(bonds[:,rowspec.index('subject_resid')])
		residue_listings[sn] = residues
		lipid_resnames = np.unique(
			data_bonds[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		chargings_this = {}
		for resnum in residues:
			for target in lipid_resnames:
				partners = np.unique(bonds[np.where(np.all((
					bonds[:,rowspec.index('subject_resid')]==resnum,
					bonds[:,rowspec.index('target_resname')]==target),
					axis=0))[0]][:,rowspec.index('target_resid')])
				bond_this = np.zeros(len(obs))
				for partner in partners:
					rows = np.where(np.all((
						bonds[:,rowspec.index('subject_resid')]==resnum,
						bonds[:,rowspec.index('target_resname')]==target,
						bonds[:,rowspec.index('target_resid')]==partner,
						),axis=0))[0]
					bond_this += obs[:,rows].any(axis=1)*1.0
				chargings_this[(resnum,target)] = bond_this
		chargings[sn] = chargings_this
	# package
	self.data = dict(chargings=chargings,residues=residues)
	return

def bond_combinator(*offs):
	"""Combine bonds and observations between replicates."""
	bonds_all = np.concatenate([o['bonds'] for o in offs])
	idx,counts = uniquify(bonds_all)
	bonds_master = bonds_all[idx]
	nframes = np.cumsum([o['observations'].shape[0] for o in offs])
	obs = np.zeros((sum(nframes),len(bonds_master)))
	for onum,off in enumerate(offs):
		# frames are just contatenated
		frameslice = slice(0 if onum==0 else nframes[onum-1],nframes[onum])
		for bnum,bond in enumerate(off['bonds']):
			status('combining bonds from replicate %d'%onum,i=bnum,looplen=len(off['bonds']))
			inds = np.where(np.all(bonds_master==bond,axis=1))[0]
			if len(inds)!=1: raise Exception
			try: obs[frameslice,inds[0]] += off['observations'][:,bnum]
			except:
				import ipdb;ipdb.set_trace()
	return dict(observations=obs,bonds=bonds_master)

def plot_charging_histograms_by_lipid(sup):
	"""
	"""
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
		ax.set_xlabel('salt bridges')
		ax.set_xticks(np.arange(peak))
		ax.set_title(resid_to_title[res])
	legend = axes[-1].legend(loc='upper left',bbox_to_anchor=(1.05,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor(sup.ps.get('legend_edge_color','k'))
	frame.set_facecolor('white')
	suptitle = fig.suptitle('mDia2 %s ($%.1f\AA$) with %s'%(
		sup_name,cutoff_this,work.vars['names']['short'].get(target_resname,target_resname)),fontsize=16)
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta={},extras=[legend,suptitle])

def plot_charging_histograms_stacked(sup,show_zero=False,total=False,
	emphasis=True,emphasis_base=0.35,square_axes=False,invert=False,mode='landscape'):
	"""
	Panels by residue and simulation with stacked histograms for lipid types.
	"""
	total_label = 'total'
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
	#! hardcoding lipid types here
	lipid_types_ref = ['PI2P','DOPS','DOPE']
	ymax = 0
	for rnum,residue in enumerate(residues):
		for snum,sn in enumerate(sup.sns):
			ax = axes[rnum][snum]
			lipid_types = list(set([l for r,l in chargings[sn].keys() 
				if r==residue or residue==total_label]))
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
						if invert:
							ax.barh([x],[y],left=counts_prev[ii],height=1.0,color=colors[lipid],zorder=2,
								alpha=alpha_func(ii,peak),
								label=work.vars['names']['short'][lipid] if ii==len(bond_values)-1 else None)
						else:
							ax.bar([x],[y],bottom=counts_prev[ii],width=1.0,color=colors[lipid],zorder=2,
								alpha=alpha_func(ii,peak),
								label=work.vars['names']['short'][lipid] if ii==len(bond_values)-1 else None)
					counts_prev += counts_norm.astype(float)[xbin_filt]
					ymax = max([ymax,max(counts_prev)])
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
					ax.set_ylabel(sup_name)
				ax.set_ylim((0 if show_zero else 1,peak+1))
				for i in bond_values[1:]: ax.axhline(i,c=sup.ps.get('frame_color','k'),lw=1,zorder=1)
			elif not invert and mode=='landscape':
				if rnum==0: ax.set_ylabel('rate')
				if snum<len(sup.sns)-1: ax.set_xticklabels([])
				else: 
					ax.set_xticks(bond_values+0.5)
					ax.set_xticklabels(bond_values.astype(int))
					ax.set_xlabel(sup_name)
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
	for ax,rnum,snum in [(axes[rnum][snum],rnum,snum) 
		for rnum,residue in enumerate(residues) for snum,sn in enumerate(sup.sns)]+axes_alt:
		if invert and mode=='landscape': 
			ax.set_xlim((0,ymax))
			if snum!=len(sup.sns)-1: ax.set_xticks([]) #ax.set_xticks(np.arange(0,ymax,0.2))
		elif not invert and mode=='landscape':
			ax.set_ylim((0,ymax))
			if rnum!=0: ax.set_yticks([])
			#! note that we cannot use square with the twinx used to put labels on the right
			#! ... see https://github.com/matplotlib/matplotlib/issues/7654
			#! note that I gutted the twinx method for putting the label on the right because 
			#! ... it was absolutely impossible to remove tick marks after twinx in the invert layout
			#! ... mode which is an extremely annoying feature and caused me to discard twinx entirely
			if square_axes: ax.set_aspect(peak/ymax)
		else: raise Exception
	extras.append(fig.suptitle('mDia2 %s ($%.1f\AA$)'%(sup_name,cutoff_this),fontsize=16))
	picturesave(sup.make_output_name(),work.plotdir,
		backup=False,version=True,meta={},extras=extras)

if __name__=='__main__':

	# requested plots
	orders = {
		# DEPRECATED BECAUSE UGLY
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
			'data_bonds':'data_salt','combine_replicates':True},},
		# panels of sns by residue with stacked bar histograms
		'charging salt merged':{
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':['stacked','salt','merged']},
			'plot':{'function':'plot_charging_histograms_stacked','sup_name':'salt bridges',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			'data':{'function':'compute_chargings','data_bonds':'data_salt','combine_replicates':True},
			'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
				'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9'}},
		'charging cut_2.2 merged':{
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':['stacked','cut_2.2','merged']},
			'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			'data':{'function':'compute_chargings','data_bonds':'data_contacts',
				'contacts_cutoff':2.2,'combine_replicates':True},
			'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
				'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9'}},
		'charging cut_5.0 merged':{
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':['stacked','cut_5.0','merged']},
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
			'plot':{'function':'plot_charging_histograms_stacked','sup_name':'salt bridges',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			'data':{'function':'compute_chargings','data_bonds':'data_salt'},
			'plotspec':{'figsize':(8,8),'legend_edge_color':'w',
				'wspace':0.1,'hspace':0.1,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.25,1.,1.)}},
		'charging cut_2.2':{
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':['stacked','cut_2.2']},
			'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			'data':{'function':'compute_chargings','data_bonds':'data_contacts',
				'contacts_cutoff':2.2},
			'plotspec':{'figsize':(12,8),'legend_edge_color':'w',
				'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.5,1.,1.)}},
		'charging cut_5.0':{
			'sns':sns_mdia2_ordering,
			'namer':{'base':'charging','mods':['stacked','cut_5.0']},
			'plot':{'function':'plot_charging_histograms_stacked','sup_name':'contacts',
				'kwargs':{'show_zero':False,'emphasis':True,'emphasis_base':0.5,'invert':True,'total':True}},
			# somewhat slow because of high cutoff hence larger bond list
			'data':{'function':'compute_chargings','data_bonds':'data_contacts',
				'contacts_cutoff':5.0},
			'plotspec':{'figsize':(16,8),'legend_edge_color':'w',
				'wspace':0.2,'hspace':0.2,'frame_color':'#A9A9A9','legend_bbox':(0.0,-1.5,1.,1.)}},
		}

	#! dev
	i = 'charging salt' ; orders = {i:orders[i]}
	# loop over requested plots
	for name,spec in orders.items(): 
		orders[name]['plotter'] = ActinlinkPlotter(name,spec)
		orders[name]['plotter'].render()
