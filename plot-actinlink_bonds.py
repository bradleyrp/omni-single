#!/usr/bin/env python

"""
ACTINLINK BONDS
see lab notebook for details
"""

import time
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from base.tools import status,framelooper

#---block: what to plot
routine = [
	'contact_map',
	'histograms','bars'][-1:]
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
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,'fs_legend_title':20,
	'binary_color_intensity':0.5}
residue_colors = dict([(name,residue_type_colors[residue_types[name]]) for name in residue_types])
ticks_font = mpl.font_manager.FontProperties(family='Latin Modern Mono',style='normal',
	size=14,weight='normal',stretch='normal')

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

#---block: post-post processing
def make_postdat():
	global sns,data_hbonds,data_salt,data_contacts,bond_mappings
	#---easy lookup for multiple upstream data types
	pointers = {'hbonds':data_hbonds,'contacts':data_contacts,'salt':data_salt}
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
			this_data = pointers[bond['upstream']]
			#---raw data are timewise bond counts for resid,lipid_name pairs
			raw = dict([((int(i[0]),i[1]),j)
				for i,j in zip(
					this_data[sn]['data']['pairs_resid_resname'],
					#---this line specifies the two kinds of data to extract from contacts: either
					#---...either the explicit or reduced bonds
					this_data[sn]['data'][bond['post_key']])])
			#---post is the raw data organized by lipid_name into matrices suitable for plotting
			postdat[sn][bond['name']] = dict([(resname,
				np.array([raw[(r,resname)] 
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
			resnames = postdat[sn]['subject_residues_resnames']
			if rr==0: 
				ax.set_ylabel(sn_title(sn),fontsize=plotspec['fs_ylabel'])
				if work.plots[plotname].get('settings',{}).get('show_residue_names',True):
					#---never set the ticklabels before the ticks
					ax.set_yticks(np.arange(len(resnames))+0.5)
					ax.set_yticklabels([residue_codes[r] for r in resnames],
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
				image_data = np.array([mpl.cm.__dict__[
					residue_colors.get(resnames[rnum],'Greys')](
						postdat[sn][bond_name][resname][rnum].astype(float)/ceiling).T[:-1].T 
					for rnum in np.arange(postdat[sn][bond_name][resname].shape[0])])
			else: image_data = postdat[sn][bond_name][resname]
			duration = (postdat[sn]['times'].max()-postdat[sn]['times'].min())*10**-3
			xtick_interval = work.plots[plotname].get('settings',{}).get('xtick_interval',
				plotspec['time_tick_interval'])
			ax.set_xticks(np.arange(xtick_interval,duration,xtick_interval))
			ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
			#---! removed division by the following, designed to make things lighter when reduced
			#---! ...{'explicit':1.0,'reduced':plotspec['binary_color_intensity']}[mode]
			im = ax.imshow(image_data,
				origin='lower',interpolation='nearest',extent=[0,duration,0,len(resids)],
				aspect=float(duration)/len(resids),
				vmax=ceiling,vmin=0,cmap=mpl.cm.__dict__[cmap])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		#---color bar on the second row
		if ss==1:
			axins = inset_axes(ax,width="5%",height="100%",loc=3,
				bbox_to_anchor=(1.1,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
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
	axes,fig = square_tiles(len(sns),figsize=plotspec['figsize'],hspace=0.4,wspace=0.4)
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
					raw2 = raw.mean(axis=1)
					colors_this = [mpl.cm.__dict__[residue_type_colors[residue_types[name]]](
						plotspec['legend_color_strength']) 
						for name in postdat[sn]['subject_residues_resnames']]
					ax.bar(np.arange(raw.shape[0]),raw2,color=colors_this)
					ys = raw2
					subs = raw2>0
					ax.errorbar(np.arange(raw.shape[0])[subs],ys[subs],yerr=raw.std(axis=1)[subs],
						fmt='none',color='k',alpha=0.5)
					maxcount = max([max(raw2),maxcount])
			else: raise Exception('unclear style')
	if style=='bars' and bond_name=='hbonds' and False:
		import ipdb;ipdb.set_trace()
	if style=='histograms':
		for snum,sn in enumerate(sns):
			ax = axes[snum]
			ax.set_ylim(0,maxcount*1.05)
			if False: ax.set_xlim(0,maxx)
			legend = axes[-1].legend(loc='upper left',fontsize=plotspec['fs_legend'],
				bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True)
			ax.set_xlabel('bonds per frame')
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_ylabel('probability')
	elif style=='bars':
		for snum,sn in enumerate(sns):
			ax = axes[snum]
			ax.set_ylim(0,maxcount*1.05)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			#---never set the ticklabels before the ticks
			resnames = postdat[sn]['subject_residues_resnames']
			ax.set_xticks(np.arange(len(resnames)))
			ax.set_xticklabels([residue_codes[r] for r in resnames],fontsize=plotspec['fs_xticks'])
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
	picturesave(fn,work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: counting hydrogen bonds
def bond_counter(resid,resname_set):
	"""
	Adapted from the explicit section of the contacts.py function called "count_reduced_contact".
	The hydrogen bonding currently uses the same rowspec as the contacts.
	"""
	global bonds,obs,rowspec
	#---filter the observations by the protein residue (subject_resid) and target resname
	#---...providing a result
	which = np.where(np.all((bonds[:,rowspec.index('subject_resid')].astype(int)==resid,np.in1d(bonds[:,rowspec.index('target_resname')],resname_set)),axis=0))
	result = obs.T[which].sum(axis=0)
	return result

#---block: stolen from contacts.py
def basic_compute_loop(compute_function,looper,run_parallel=True,debug=False):
	"""Canonical form of the basic compute loop."""
	start = time.time()
	if run_parallel:
		incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
			delayed(compute_function,has_shareable_memory)(**looper[ll]) 
			for ll in framelooper(len(looper),start=start))
	else: 
		incoming = []
		for ll in framelooper(len(looper)):
			incoming.append(compute_function(**looper[ll]))
	return incoming

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns,(data_contacts,calc) = work.sns(),plotload('contacts',work)
	_,(data_hbonds,_) = work.sns(),plotload('hydrogen_bonding',work)
	_,(data_salt,_) = work.sns(),plotload('salt_bridges',work)
	
	#---set the cutoff in the yaml file to make this plot because we pull from multiple upstream sources
	cutoff = calc['calcs']['specs']['cutoff']
	#---map data type onto keys
	bond_mappings = [
		{'name':'explicit','post_key':'counts_resid_resname','upstream':'contacts'},
		{'name':'reduced','post_key':'counts_resid_resname_singleton','upstream':'contacts'},
		{'name':'hbonds','post_key':'hbonds_compacted','upstream':'hbonds'},]
	#---specify the upstream data and naming scheme
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	#---common colors
	import brewer2mpl
	colors = {'DOPC':'blue','DOPS':'red','POP2':'magenta','PI2P':'magenta',
		'all lipids':'black','DOPE':'blue'}
	colors.update(**dict([(sn,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[sns.index(sn)]) 
		for sn in sns]))
	lipid_label = lambda x: dict([(i,'$$\mathrm{{PIP}_{2}}$$') 
		for i in work.vars['selectors']['resnames_PIP2']]).get(x,x)
	sn_title = lambda sn: '%s%s'%(work.meta[sn].get('label',re.sub('_','-',sn)),
		'' if not work.meta[sn].get('cholesterol',False) else '\n(cholesterol)')

	#---post-post processing
	if 'postdat' not in globals(): 
		hydrogen_bond_compactor()
		postdat = make_postdat()

#---block: all chemical bonds on one contact map
plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
	'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,
	'label_color_strength':1.0,'fs_legend_title':20,
	'binary_color_intensity':0.5,'figsize':(14,14),'time_tick_interval':20,
	'fs_xticks':11,'fs_yticks':11}
plotspec_small = {'figsize':(8,8),'fs_legend_title':14,'fs_legend':12,'fs_xticks':12,'fs_yticks':12}
specs = {'all':{},
	'mDia2':{'lipid_resnames':['DOPS','PI2P'],'tag':'.PS_PIP2',
		'sns':['mdia2bilayer_nochl2','mdia2bilayerphys'],
		'plotspec':dict(plotspec_small)},
	'gelsolin':{'lipid_resnames':['DOPS','PI2P'],'tag':'.PS_PIP2',
		'sns':['gelbilayer_nochl','gelbilayerphys'],'plotspec':dict(plotspec_small)},}
for bond in bond_mappings:
	for specname,spec in specs.items():
		this_plotspec = dict(plotspec)
		this_plotspec.update(**spec.get('plotspec',{}))
		kwargs = dict(postdat=postdat,bond_name=bond['name'],plotspec=this_plotspec)
		kwargs.update(sns=spec.get('sns',sns))
		kwargs.update(lipid_resnames=spec.get('lipid_resnames',None))
		figname = {'explicit':'contacts','reduced':'contacts','hbonds':'hbonds'}[bond['name']]
		if 'contact_map' in routine: 
			fn = 'fig.%s.%s%s.%s%s'%(figname,specname,spec.get('tag',''),bond['name'],
				'.cutoff_%.1f'%cutoff if bond['name'] in ['explicit','reduced'] else '')
			if bond['name']=='hbonds':
				kwargs.update(legend_title='hydrogen bonds',cbar_title='counts')
			colorstreak_contact_map(fn=fn,**kwargs)
		if 'histograms' in routine and spec=={}:
			plotspec_this = dict(kwargs['plotspec'])
			plotspec_this.update(**{'figsize':(8,8)})
			kwargs.update(plotspec=plotspec_this)
			fn = 'fig.%s_%s.%s%s.%s%s'%(figname,'histogram',specname,spec.get('tag',''),bond['name'],
				'.cutoff_%.1f'%cutoff if bond['name'] in ['explicit','reduced'] else '')
			plot_histograms_or_occupancy(style='histograms',fn=fn,**kwargs)
		if 'bars' in routine and spec=={}:
			plotspec_this = dict(kwargs['plotspec'])
			plotspec_this.update(**{'figsize':(8,8)})
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
