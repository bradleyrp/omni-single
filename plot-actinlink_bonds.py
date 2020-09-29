#!/usr/bin/env python

"""
ACTINLINK BONDS
see lab notebook for details
"""

import time
from joblib import Parallel,delayed
#from joblib.pool import has_shareable_memory
from base.tools import status,framelooper,dictsum
from base.compute_loop import basic_compute_loop
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---block: SETTINGS
routine = ['contact_map','histograms','bars','histogram_summary','charging_curves'][:3]
avail = ['contact_maps','histograms','bars']
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
def color_by_simulation(sn):
	four_colors = ['#B2D732','#FE2712','#347B98','#70d09f']
	# tetradic colors via https://www.sessions.edu/color-calculator/
	four_colors = ['#ff6c28','#28c6ff','#a928ff','#ffd828']
	refs = ['^mdia2bilayer_nochl','^mdia2bilayer(?!(_nochl|[0-9]))','^mdia2bilayer10','^mdia2bilayer30']
	ref_to_color = dict(zip(refs,four_colors))
	matches = [color for ref,color in ref_to_color.items() if re.match(ref,sn)]
	if len(matches)>1: raise Exception
	else: return matches[0]
sns_mdia2_ordering = ['mdia2bilayer_nochl2','mdia2bilayer_nochl3','mdia2bilayerphys',
	'mdia2bilayerphys2','mdia2bilayer10','mdia2bilayer10_2','mdia2bilayer30','mdia2bilayer30_2']
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,
	'label_color_strength':1.0,'fs_legend_title':20,
	'binary_color_intensity':0.5,'time_tick_interval':20,
	'fs_xticks':11,'fs_yticks':6,'fs_xticks_bars':6}
plotspec_small = {'fs_legend_title':14,'fs_legend':12,'fs_xticks':12,'fs_yticks':8}
#---define salt bridges
valid_salt_bridges = [
	{'resname':'ARG','atoms':['NH1','NH2']},
	{'resname':'LYS','atoms':['NZ']},
	{'resname':'HIS','atoms':['ND1','NE2']},]
#---special subselections
special_protein_parts = {'nwaspbilayernochl':np.arange(0,22+1).astype(int)}

#---which combinations of simulations to plot, and how to name them 
specs = {
	'all':{'plotspec':dictsum(plotspec,dict(figsize=(8,16),fs_ylabel=10))},}
#! focusing on mDia2, the only component of the focus collection, for now
specs_others = {'mDia2':{'lipid_resnames':['DOPS','PI2P'],'tag':'.PS_PIP2',
		'sns':['mdia2bilayer_nochl2','mdia2bilayerphys'],
		'plotspec':dictsum(plotspec_small,dict(figsize=(8,8),fs_ylabel=10))},
	'gelsolin':{'lipid_resnames':['DOPS','PI2P'],'tag':'.PS_PIP2',
		'sns':['gelbilayer_nochl','gelbilayerphys'],
		'plotspec':dictsum(plotspec_small,dict(figsize=(8,8),fs_ylabel=10))},
	'nwasp':{'lipid_resnames':['DOPS','PI2P'],'tag':'.PS_PIP2',
		'sns':['nwasppeptide2','nwaspbilayermut1','nwaspbilayernochl','nwasppeptideb'],
		'plotspec':dictsum(plotspec_small,dict(figsize=(8,14),fs_ylabel=10))},}

#---which types of data to plot
bond_mappings = [
	{'name':'reduced','post_key':'counts_resid_resname_singleton','upstream':'contacts'},
	{'name':'explicit','post_key':'counts_resid_resname','upstream':'contacts'},
	{'name':'hbonds','post_key':'hbonds_compacted','upstream':'hbonds'},
	{'name':'salt','post_key':'salt_compacted','upstream':'contacts'},
	#---! do not change the order here. hbonds comes first to make the merge for raws below
	{'name':'combo_salt_hbonds','post_key':['hbonds_compacted','salt_compacted'],
		'upstream':['hbonds','contacts']},][:]

#---block: mimic the coda in contacts.py
def hydrogen_bond_compactor():
	global data_hbonds,bonds,obs
	for sn in sns:
		#---custom mapping for collecting hydrogen bonds
		bonds,obs = [data_hbonds[sn]['data'][i] for i in ['bonds','observations']]
		#---! get resids for the protein and lipid_resnames from contact maps
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resids = data_contacts[sn]['data']['subject_residues_resids']
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		#---compute loop
		looper = [{'resid':resid,'resname_set':resname_set} 
			for resid in resids for resname_name,resname_set in resname_combos]
		compute_function = bond_counter
		incoming = basic_compute_loop(compute_function,looper,run_parallel=True)
		#---tacking on compacted data to mimic the form of the contact maps
		data_hbonds[sn]['data']['hbonds_compacted'] = np.array(incoming)
		data_hbonds[sn]['data']['pairs_resid_resname'] = np.array([(resid,resname_name) 
			for resid in resids for resname_name,resname_set in resname_combos]).astype(str)

#---block: filter out the salt bridges
def salt_bridge_filter():
	global data_contacts,bonds,obs,valid_salt_bridges
	for sn in sns:
		#---filter the bonds and observations from contact maps
		bonds_all = data_contacts[sn]['data']['bonds']
		obs_all = data_contacts[sn]['data']['observations']
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
		status('salt nbonds for %s is %d'%(sn,len(salt_inds)))
		#---! get resids for the protein and lipid_resnames from contact maps
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resids = data_contacts[sn]['data']['subject_residues_resids']
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		#---compute loop
		looper = [{'resid':resid,'resname_set':resname_set} 
			for resid in resids for resname_name,resname_set in resname_combos]
		compute_function = bond_counter
		incoming = basic_compute_loop(compute_function,looper,run_parallel=True)
		#---tacking on compacted data to mimic the form of the contact maps
		data_contacts[sn]['data']['salt_compacted'] = np.array(incoming)
		if False: data_contacts[sn]['data']['pairs_resid_resname'] = np.array([(resid,resname_name) 
			for resid in resids for resname_name,resname_set in resname_combos]).astype(str)

#---block: post-post processing
def make_postdat():
	global sns,data_hbonds,data_contacts,bond_mappings
	#---easy lookup for multiple upstream data types
	pointers = {'hbonds':data_hbonds,'contacts':data_contacts}
	#---postdat is indexed by simulation, then bond type
	postdat = dict([(sn,{}) for sn in sns])
	#---loop over simulations
	for sn in sns:
		#---common features
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		lipid_resnames = list(zip(*resname_combos))[0]
		postdat[sn]['lipid_resnames'] = lipid_resnames
		postdat[sn]['resids'] = data_contacts[sn]['data']['subject_residues_resids']
		postdat[sn]['subject_residues_resnames'] = data_contacts[sn]['data']['subject_residues_resnames']
		postdat[sn]['times'] = data_contacts[sn]['data']['times']
		#---loop over bond mappings
		for bond in bond_mappings:
			if type(bond['post_key']) in str_types: bond_post_keys = [bond['post_key']]
			else: bond_post_keys = bond['post_key']
			if type(bond['upstream']) in str_types: bond_upstreams = [bond['upstream']]
			else: bond_upstreams = bond['upstream']
			raws = []
			if len(bond_post_keys)!=len(bond_upstreams): raise Exception('invalid bond mapping')
			for bond_post_key,bond_upstream in zip(bond_post_keys,bond_upstreams):
				this_data = pointers[bond_upstream]
				#---raw data are timewise bond counts for resid,lipid_name pairs
				raw = dict([((int(i[0]),i[1]),j)
					for i,j in zip(
						this_data[sn]['data']['pairs_resid_resname'],
						#---this line specifies the two kinds of data to extract from contacts: either
						#---...either the explicit or reduced bonds
						this_data[sn]['data'][bond_post_key])])
				raws.append(raw)
			#---! valid frames may be slightly different for some reason
			if len(bond_post_keys)==2:
				#---! this fix only applies for merging salt bridges and hbonds, obv
				valid_frames = np.intersect1d(data_hbonds[sn]['data']['valid_frames'],
					data_contacts[sn]['data']['valid_frames'])
				#---produce a mapping to the common frame indices
				mappings_vf = [np.concatenate([np.where(valid_frames==i)[0] for i in 
					#---! order must match bond_mappings! code an exception to ensure this
					d[sn]['data']['valid_frames']]) for d in [data_hbonds,data_contacts]]
			else: 
				valid_frames = this_data[sn]['data']['valid_frames']
				mappings_vf = [slice(None,None)]
			#---post is the raw data organized by lipid_name into matrices suitable for plotting
			#---note that a merge is performed below if we have two data types
			postdat[sn][bond['name']] = dict([(resname,
				np.array([np.sum([raw[(r,resname)][mappings_vf[rr]] for rr,raw in enumerate(raws)],axis=0)
				for r in postdat[sn]['resids']])) for resname in postdat[sn]['lipid_resnames']])
	return postdat

#---block: contact map code
def colorstreak_contact_map(sns,postdat,bond_name,plotspec,fn,**kwargs):
	"""
	Final resting place for the contact map code.
	"""
	do_crazy_colors = kwargs.get('do_crazy_colors',True)
	cmap = kwargs.get('cmap','Greys')

	#---same residue comparisons
	lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
	if not len(lipid_resnames)==1: raise Exception('cannot use different resnames across simulations')
	else: lipid_resnames = lipid_resnames[0]
	if kwargs.get('lipid_resnames',None):
		lipid_resnames = [i for i in lipid_resnames if i in kwargs['lipid_resnames']]
	ceiling = float(max([max([i.max() for i in postdat[sn][bond_name].values()]) for sn in sns]))

	#---plot
	axes,fig = panelplot(
		layout={'out':{'hspace':0.05,'wspace':0.05,'grid':[1,len(lipid_resnames)]},
		'ins':[{'hspace':0.2,'wspace':0.05,'grid':[len(sns),1]} 
		for i in lipid_resnames]},figsize=plotspec['figsize'])
	for ss,sn in enumerate(sns):
		for rr,resname in enumerate(lipid_resnames):
			ax = axes[rr][ss]
			#---override the protein selections
			resnums = special_protein_parts.get(sn,np.arange(postdat[sn][bond_name][resname].shape[0]))
			resnames = postdat[sn]['subject_residues_resnames']
			if rr==0: 
				ax.set_ylabel(sn_title(sn,which='label_compact'),fontsize=plotspec['fs_ylabel'])
				if True: #! hack work.plots[plotname].get('settings',{}).get('show_residue_names',True):
					#---never set the ticklabels before the ticks
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
			resids = postdat[sn]['resids']
			if do_crazy_colors:
				scores = np.array([postdat[sn][bond_name][resname][rnum].astype(float)/ceiling 
					for rnum in resnums])
				image_data = np.array([mpl.cm.__dict__[residue_colors.get(resnames[rnum],'Greys')](
					postdat[sn][bond_name][resname][rnum].astype(float)/ceiling).T[:-1].T 
					for rnum in resnums])
				#---note that the crazy color color maps do not go identically to zero
				#---...so we white them out
				image_data[np.where(scores==0.0)] = [1.,1.,1.]
			else: image_data = postdat[sn][bond_name][resname][resnums]
			duration = (postdat[sn]['times'].max()-postdat[sn]['times'].min())*10**-3
			xtick_interval = plotspec['time_tick_interval'] 
			"""
			removing xtick_interval settings
			... work.plots[plotname].get('settings',{}).get('xtick_interval',plotspec['time_tick_interval'])
			ax.set_xticks(np.arange(xtick_interval,duration,xtick_interval))
			"""
			ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=4,prune='lower'))
			ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
			#---! removed division by the following, designed to make things lighter when reduced
			#---! ...{'explicit':1.0,'reduced':plotspec['binary_color_intensity']}[mode]
			im = ax.imshow(image_data,
				origin='lower',interpolation='nearest',extent=[0,duration,0,len(resnums)],
				aspect=float(duration)/len(resnums),
				vmax=ceiling,vmin=0,cmap=mpl.cm.__dict__[cmap])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		#---color bar on the second row
		if ss==len(sns)-1:
			axins = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.3,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
			cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
			cbar_title = kwargs.get('cbar_title','score')
			cbar.ax.set_title(cbar_title)

	#---title on the legend
	title = kwargs.get('legend_title','contacts\n'+r'$\mathrm{r \leq %.1f \AA}$'%(cutoff))
	patches,labels = zip(*[(mpl.patches.Rectangle((0,0),1.0,1.0,
		fc=mpl.cm.__dict__[residue_type_colors[kind]](plotspec['legend_color_strength'])),kind) 
		for kind in residue_type_colors])
	legend = axes[-1][0].legend(patches,labels,loc='upper left',fontsize=plotspec['fs_legend'],
		bbox_to_anchor=(1.05,0.0,1.,1.),title=title,shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	plt.setp(legend.get_title(),fontsize=plotspec['fs_legend_title'])
	picturesave(fn,work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: histogram code
def plot_histograms_or_occupancy(sns,postdat,bond_name,plotspec,fn,**kwargs):
	"""
	Final resting place for the contact map code.
	"""
	max_nbins = kwargs.get('max_nbins',40)
	landscape = kwargs.get('landscape',False)
	style = kwargs.get('style','histogram')
	#---same residue comparisons
	lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
	if not len(lipid_resnames)==1: raise Exception('cannot use different resnames across simulations')
	else: lipid_resnames = lipid_resnames[0]
	if kwargs.get('lipid_resnames',None):
		lipid_resnames = [i for i in lipid_resnames if i in kwargs['lipid_resnames']]
	ceiling = float(max([max([i.max() for i in postdat[sn][bond_name].values()]) for sn in sns]))

	#---plot
	maxcount,maxx = 0,0
	axes,fig = square_tiles(len(sns),figsize=plotspec['figsize'],hspace=0.4,wspace=0.4,favor_rows=landscape)
	for snum,sn in enumerate(sns):
		ax = axes[snum]
		ax.set_title(sn_title(sn))
		for rr,resname in enumerate(lipid_resnames):
			raw = postdat[sn][bond_name][resname]
			if style=='histograms':
				#---histograms show all bonds (summed over protein residues)
				raw2 = raw.sum(axis=0)
				#---! deprecated: step_size = 10**(np.ceil(np.log10(raw2.max()))-2)
				step_size = 1.0
				bins = np.arange(0,raw2.max()+step_size,step_size)
				#---limit the number of bins
				if len(bins)>max_nbins: bins = np.linspace(0,raw2.max()+step_size,max_nbins)
				maxx = max([maxx,raw2.max()+step_size])
				counts,bins = np.histogram(raw2,bins=bins)
				subs = counts>0
				mids = (bins[1:]+bins[:-1])/2.0
				ys = counts[subs].astype(float)/counts[subs].sum()
				ax.plot(bins[1:][subs],ys,'.-',
					label=lipid_label(resname),color=colors[resname])
				if bond_name=='hbonds': ax.set_xticks(bins[1:])
				maxcount = max([max(ys),maxcount])
			elif style=='bars':
				lipid_subselection = kwargs.get('lipid_bar_selection',None)
				if not lipid_subselection: 
					raise Exception('bar style requires a lipid_bar_selection')
				if resname==lipid_subselection:
					#---bars show the mean number of contacts per frame
					#---override the protein selections
					resnums = special_protein_parts.get(sn,np.arange(
						postdat[sn][bond_name][resname].shape[0]))
					raw2 = raw.mean(axis=1)
					colors_this = [mpl.cm.__dict__[residue_type_colors[residue_types[name]]](
						plotspec['legend_color_strength']) 
						for name in postdat[sn]['subject_residues_resnames']]
					ax.bar(resnums,raw2[resnums],color=colors_this)
					ys = raw2
					subs = raw2>0
					ax.errorbar(np.arange(raw.shape[0])[subs],ys[subs],yerr=raw.std(axis=1)[subs],
						fmt='none',color='k',alpha=0.5)
					maxcount = max([max(raw2),maxcount])
			else: raise Exception('unclear style')
	if style=='histograms':
		for snum,sn in enumerate(sns):
			ax = axes[snum]
			ax.set_ylim(0,maxcount*1.05)
			ax.set_xlim(0,maxx)
			ax.set_xlabel('bonds per frame')
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_ylabel('probability')
		legend = axes[-1].legend(loc='upper left',fontsize=plotspec['fs_legend'],
			bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True)
	elif style=='bars':
		for snum,sn in enumerate(sns):
			ax = axes[snum]
			ax.set_ylim(0,maxcount*1.05)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			#---override the protein selections
			resnums = special_protein_parts.get(sn,np.arange(postdat[sn][bond_name][resname].shape[0]))
			#---never set the ticklabels before the ticks
			resnames = postdat[sn]['subject_residues_resnames']
			ax.set_xticks(np.arange(len(resnums)))
			ax.set_xticklabels([residue_codes[r] for r in resnames[resnums]],
				fontsize=plotspec['fs_xticks_bars'])
			for label in ax.get_xticklabels():
				label.set_color(mpl.cm.__dict__[
					residue_colors.get(residue_codes_reverse[label._text],
						'Greys')](plotspec['label_color_strength'])[:-1])
			ax.set_ylabel('bonds per frame')
		patches,labels = zip(*[(mpl.patches.Rectangle((0,0),1.0,1.0,
			fc=mpl.cm.__dict__[residue_type_colors[kind]](plotspec['legend_color_strength'])),kind) 
			for kind in residue_type_colors])
		legend = axes[-1].legend(patches,labels,loc='upper left',fontsize=plotspec['fs_legend'],
			bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True)
		frame = legend.get_frame()
		frame.set_edgecolor('black')
		frame.set_facecolor('white')
	else: raise Exception('unclear style')
	if not legend: raise Exceptin('dev')
	picturesave(fn,work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: summary plots
def plot_histogram_summary(fn_base,bar_width=0.9):
	"""
	Bar plots over simulations of hydrogen bonds plus salt bridges with one panel for each lipid (or all).
	"""
	sns = sns_mdia2_ordering
	resnames = list(set([postdat[sn]['lipid_resnames'] for sn in sns]))[0]
	plotspec = dict(figsize=(12,8))

	# prepare metadata for different plot styles
	plotspecs = {
		'':{'figsize':(12,8),'time':'all'},
		'.nbins%d'%2:{'figsize':(12,8),'time':'all','n_segments':2},
		'.nbins%d'%20:{'figsize':(12,8),'time':'all','n_segments':20,'error_bars':False},
		'.compare':{'figsize':(12,8),'time':'all','breakout':True},}

	def seg_to_slice(sn,segnum,n_segments,nframes=None):
		if not nframes: nframes = data_contacts[sn]['data']['nframes']
		# fishing with dynamite. divide all frames into groups and make sure we use each frame
		return np.where(segnum==np.floor(np.arange(nframes)/(nframes/float(n_segments))).astype(int))[0]
	# loop over plot styles
	for name_suffix,plotspec in plotspecs.items():
		axes,fig = square_tiles(len(resnames),figsize=plotspec['figsize'],hspace=1.0,wspace=0.4)
		n_segments = plotspec.get('n_segments',1)
		if n_segments>1 and plotspec.get('breakout',False): raise Exception('invalid plotspec')
		for snum in range(len(sns)-1):
			for rr,resname in enumerate(resnames):
				axes[rr].axvline(snum+1.0-0.5,lw=1,c='k',zorder=0)
		for segnum in range(n_segments):
			means,stds = np.array([[[
				getattr(postdat[sn]['combo_salt_hbonds'][resname].sum(axis=0)[
					seg_to_slice(sn,segnum,n_segments,
						# we need to send along nframes because of valid frames jitter from kDTree
						nframes=len(postdat[sn]['combo_salt_hbonds'][resname].T))],attr)() 
				for sn in sns] for resname in resnames] for attr in ['mean','std']])
			# alternate computation for the salt bridges only in case we want the breakout
			means_salt,stds_salt = np.array([[[
				getattr(postdat[sn]['salt'][resname].sum(axis=0)[
					seg_to_slice(sn,segnum,n_segments,
						nframes=len(postdat[sn]['salt'][resname].T))],attr)() 
				for sn in sns] for resname in resnames] for attr in ['mean','std']])
			for rr,resname in enumerate(resnames):
				for snum,sn in enumerate(sns):
					width = bar_width/n_segments
					axes[rr].bar(snum+segnum*width-bar_width/2.+width/2.,means[rr][snum],
						color=color_by_simulation(sn),width=width,zorder=2)
					if plotspec.get('breakout',False): 
						axes[rr].plot([snum-0.4,snum+0.2],[means_salt[rr][snum] 
							for i in range(2)],'--',c='k',zorder=2)
					if plotspec.get('error_bars',True):
						axes[rr].errorbar(snum+segnum*width-bar_width/2.+width/2.,
							means[rr][snum],yerr=stds[rr][snum],fmt='none',color='k',alpha=1.0,zorder=3)
					if plotspec.get('breakout',False):
						eb = axes[rr].errorbar(snum-0.1,means_salt[rr][snum],yerr=stds[rr][snum],
							fmt='none',color='k',alpha=1.0,zorder=3)
						eb[-1][0].set_linestyle('-.')
				axes[rr].set_ylabel('bonds per frame\n(salt bridges + hydrogen bonds)')
		# decoration
		for anum,ax in enumerate(axes): ax.set_title(
			work.vars['names']['short'].get(resnames[anum],resnames[anum]))
		for ax in axes:
			ax.set_xlim((-0.5,len(sns)-0.5))
			ax.set_ylim((0,axes[rr].get_ylim()[1]))
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		for axnum,ax in enumerate(axes): 
			ax.set_xticks(np.arange(len(sns)))
			ax.set_xticklabels([work.meta[sn]['label'] for sn in sns],rotation=90)
		picturesave('fig.%s%s'%(fn_base,name_suffix),
			work.plotdir,backup=False,version=True,meta={},extras=[])

def plot_charging_curves():
	"""
	Charging curve equivalents for the mDia2 simulations
	"""
	#! attempting to reproduce an equivalent to the charging curves
	chargings = {}
	sns = sns_mdia2_ordering
	for sn in sns:
		bonds,obs = [data_contacts[sn]['data'][k] for k in ['bonds','observations']]
		residues = np.unique(bonds[:,rowspec.index('subject_resid')])
		lipid_resnames = np.unique(
			data_contacts[sn]['data']['bonds'][:,rowspec.index('target_resname')])
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
	#! assume same protein for clarity in rows
	resnames = postdat[sn]['subject_residues_resnames']
	axes,fig = panelplot(
		layout={'out':{'hspace':0.0,'wspace':0.0,'grid':[len(residues),1]},
		'ins':[{'hspace':0.,'wspace':0.0,'grid':[1,len(sns)]} 
		for i in residues]},figsize=(20,10))
	peak = max([np.concatenate(chargings[sn].values()).max() for sn in sns])
	for snum,sn in enumerate(sns):
		for rnum,res in enumerate(residues):
			ax = axes[rnum][snum]
			for key in chargings[sn]:
				if key[0]==res:
					ax.plot(chargings[sn][key],color=colorize(work.meta[sn],resname=key[1]),
						zorder=['DOPE','DOPS','PI2P'].index(key[1]),
						label=work.vars['names']['short'].get(key[1],key[1]))
	for snum,sn in enumerate(sns):
		for rnum,res in enumerate(residues):
			ax = axes[rnum][snum]
			if snum==0 and rnum==0:
				legend = ax.legend(loc='upper left',bbox_to_anchor=(0.0,2.5),ncol=len(residues))
				legend.get_frame().set_linewidth(2.0)
				legend.get_frame().set_edgecolor('k')
			ax.set_ylim((0,peak))
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_xticks([])
			ax.set_yticks([])
			#! assumes same ordering of resnames and residue indices
			if snum in [0,len(sns)-1]: 
				if snum==len(sns)-1: ax.yaxis.set_label_coords(1.1,0.0)
				else: ax.yaxis.set_label_coords(-0.1,0.0)
				ax.set_ylabel(residue_codes[resnames[rnum]],fontsize=26,rotation=0,
					color=mpl.cm.__dict__[residue_colors[resnames[rnum]]](0.8),labelpad=20,
					family='monospace')
			if rnum==0: ax.set_title(work.meta[sn]['label'])
			plt.setp(ax.spines.values(),color='k',alpha=0.5)
	picturesave('fig.charging_curves_mdia2',work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: counting hydrogen bonds
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

#---block: import the post-processed data	
if 'data_contacts' not in globals(): 
	sns_contacts,(data_contacts,calc_contacts) = work.sns(),plotload('contacts',work)
	sns_hbonds,(data_hbonds,calc_hbonds) = work.sns(),plotload('hydrogen_bonding',work)
	if sns_hbonds!=sns_contacts: 
		raise Exception('collections for hydrogen_bonding and contacts are not equal')
	else: sns = sns_mdia2_ordering #! sns = sns_hbonds
	#---set the cutoff in the yaml file to make this plot because we pull from multiple upstream sources
	cutoff = calc_contacts['contacts']['calcs']['specs']['cutoff']
	#---specify the upstream data and naming scheme
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	#---common colors
	import brewer2mpl
	colors = {'DOPC':'blue','DOPS':'red','POP2':'magenta','PI2P':'magenta',
		'all lipids':'black','DOPE':'blue'}
	colors.update(**dict([(sn,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[sns.index(sn)]) 
		for sn in sns]))
	lipid_label = lambda x: dict([(i,r'$\mathrm{{PIP}_{2}}$') 
		for i in work.vars.get('selectors',{}).get('resnames_PIP2',{})]).get(x,x)
	#! not sure where the weird cholesterol naming is used
	sn_title = lambda sn,which='label': '%s%s'%(work.meta[sn].get(which,re.sub('_','-',sn)),
		'' if not work.meta[sn].get('cholesterol',False) else '\n(cholesterol)')
	sn_title = lambda sn,which='label': work.meta[sn].get(which,re.sub('_','-',sn))
	#---post-post processing
	if 'postdat' not in globals(): 
		salt_bridge_filter()
		hydrogen_bond_compactor()
		postdat = make_postdat()

# block: MAIN PLOTTING LOOP
for bond in bond_mappings:
	for specname,spec in specs.items():
		this_plotspec = dict(plotspec)
		this_plotspec.update(**spec.get('plotspec',{}))
		kwargs = dict(postdat=postdat,bond_name=bond['name'],plotspec=this_plotspec)
		kwargs.update(sns=spec.get('sns',sns))
		kwargs.update(lipid_resnames=spec.get('lipid_resnames',None))
		figname = {'explicit':'contacts','reduced':'contacts',
			'hbonds':'hbonds','salt':'salt','combo_salt_hbonds':'combo'}[bond['name']]
		if 'contact_map' in routine: 
			fn = 'fig.%s.%s%s.%s%s'%(figname,specname,spec.get('tag',''),bond['name'],
				'.cutoff_%.1f'%cutoff if bond['name'] in ['explicit','reduced'] else '')
			if bond['name']=='hbonds': kwargs.update(legend_title='hydrogen bonds',cbar_title='counts')
			elif bond['name']=='salt': kwargs.update(legend_title='salt bridges',cbar_title='counts')
			colorstreak_contact_map(fn=fn,**kwargs)
		if 'histograms' in routine and specname=='all':
			plotspec_this = dict(kwargs['plotspec'])
			kwargs.update(plotspec=plotspec_this)
			fn = 'fig.%s_%s.%s%s.%s%s'%(figname,'histogram',specname,spec.get('tag',''),bond['name'],
				'.cutoff_%.1f'%cutoff if bond['name'] in ['explicit','reduced'] else '')
			plot_histograms_or_occupancy(style='histograms',fn=fn,**kwargs)
		if 'bars' in routine and specname=='all':
			plotspec_this = dict(kwargs['plotspec'])
			kwargs.update(plotspec=plotspec_this)
			#---! be careful with resnames
			#---! ALSO WHERE IS THE CHOLESTEROL?
			resnames = list(set([postdat[sn]['lipid_resnames'] for sn in sns]))	
			if len(resnames)!=1: raise Exception('too many sets of resnames')
			else: resnames = resnames[0]
			for resname in resnames:
				fn = 'fig.%s_%s.%s%s.%s.%s%s'%(figname,'bars',specname,spec.get('tag',''),
					re.sub(' ','_',resname),
					bond['name'],'.cutoff_%.1f'%cutoff if bond['name'] in ['explicit','reduced'] else '')
				plot_histograms_or_occupancy(lipid_bar_selection=resname,style='bars',fn=fn,**kwargs)
# independent plots
if 'histogram_summary' in routine: plot_histogram_summary('combo_summary')
if 'charging_curves' in routine: plot_charging_curves()
