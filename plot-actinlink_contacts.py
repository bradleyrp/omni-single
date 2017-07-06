#!/usr/bin/env python

"""
"""

"""
workshopping the different reductions we can do
	1. occupancy
		take the mean of the sum of the bonds over the trajectory to get the occupancy as a score
		relevant dimensions
			lipid id
			residue id
			lipid resname
			residue resname
		plot styles
			1. bar style
				bar style is the score for the entire simulation broken down along a particular dimension
				but the combinations are too numerous. see the hydrogen bonds plots where we just added proteins and they are definitely too busy right now
			2. residue style might be easier
				x-axis is residue number and then a score for different lipid types
	2. residence time
"""

routine = ['contact_map','alt'][-1:]
if is_live: 
	from ipywidgets import *
	from IPython.display import display

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns,(data,calc) = work.sns(),plotload(plotname,work)
	#---specify the upstream data
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	#---set the cutoff in the yaml file to make this plot
	cutoff = calc['calcs']['specs']['cutoff']

#---block: prepare contact map data
if 'postdat' not in globals() and 'contact_map' in routine:
	postdat = {}
	for sn in sns:
		dat = data[sn]['data']
		nframes = dat['nframes']
		bonds,obs = [dat[i] for i in ['bonds','observations']]
		#---get the list of residues
		resnames,resids = [dat[i] for i in ['subject_residues_resnames','subject_residues_resids']]
		#---collect the contacts score for each residue and each lipid type
		lipid_resnames = np.unique(bonds[:,rowspec.index('target_resname')])
		resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
			('all lipids',np.array(lipid_resnames))]
		postdat[sn] = dict(residue_lipid_time={},lipid_resnames=list(zip(*resname_combos))[0],
			resids=resids,nframes=nframes)
		#---double loop over residues and lipid types
		#---! note that this is somewhat slow and perhaps might be worthy of another calculation
		for rr,resid in enumerate(resids):
			status('counting contacts for %s'%sn,i=rr,looplen=len(resids),tag='compute')
			for resname_name,resname_set in resname_combos:
				#---! rpb notices high timewise variation in chunks so this might be really dynamic
				postdat[sn]['residue_lipid_time'][(resid,resname_name)] = obs.T[
					np.where(np.all((
						bonds[:,rowspec.index('subject_resid')].astype(int)==resid,
						np.in1d(bonds[:,rowspec.index('target_resname')],resname_set)),axis=0))
					].sum(axis=0)
	#---chompdown
	postpostdat = dict([(sn,dict([(resname,np.array([postdat[sn]['residue_lipid_time'][(r,resname)] 
		for r in postdat[sn]['resids']])) for resname in postdat[sn]['lipid_resnames']])) for sn in sns])

#---block: draw the contact maps
if 'contact_map' in routine:

	#---plot settings
	figsize = (14,14)
	cmap = 'Greys'
	do_crazy_colors = True
	residue_types = {'ARG':'basic','HIS':'basic','LYS':'basic',
		'ASP':'acidic','GLU':'acidic','SER':'polar','THR':'polar','ASN':'polar','GLN':'polar',
		'ALA':'hydrophobic','VAL':'hydrophobic','ILE':'hydrophobic','LEU':'hydrophobic',
			'MET':'hydrophobic','PHE':'hydrophobic','TYR':'hydrophobic','TRP':'hydrophobic',
		'CYS':'special','SEC':'special','GLY':'special','PRO':'special'}
	residue_type_colors = {'basic':'Blues','acidic':'Reds','hydrophobic':'Greens',
		'polar':'Purples','special':'Oranges'}
	plotspec = {'fs_xlabel':14,'fs_ylabel':20,'fs_title':20,
		'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,'fs_legend_title':20}
	residue_colors = dict([(name,residue_type_colors[residue_types[name]]) for name in residue_types])
	ticks_font = mpl.font_manager.FontProperties(family='Latin Modern Mono',style='normal',
		size=14,weight='normal',stretch='normal')
	#---same residue comparisons
	lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
	if not len(lipid_resnames)==1: raise Exception('cannot use different resnames across simulations')
	else: lipid_resnames = lipid_resnames[0]
	ceiling = max([max([i.max() for i in postpostdat[sn].values()]) for sn in sns])
	#---plot
	axes,fig = panelplot(
		layout={'out':{'hspace':0.05,'wspace':0.05,'grid':[1,len(lipid_resnames)]},
		'ins':[{'hspace':0.2,'wspace':0.05,'grid':[len(sns),1]} for sn in sns]},figsize=figsize)
	for ss,sn in enumerate(sns):
		for rr,resname in enumerate(lipid_resnames):
			ax = axes[rr][ss]
			resnames = data[sn]['data']['subject_residues_resnames']
			if rr==0: 
				sn_title = '%s%s'%(work.meta[sn].get('label',re.sub('_','-',sn)),
					'' if not work.meta[sn].get('cholesterol',False) else '\n(cholesterol)')
				ax.set_ylabel(sn_title,fontsize=plotspec['fs_ylabel'])
				#---never set the ticklabels before the ticks
				ax.set_yticks(np.arange(len(resnames))+0.5)
				ax.set_yticklabels(resnames)
				for label in ax.get_yticklabels():
					label.set_color(mpl.cm.__dict__[residue_colors.get(str(label._text),'Greys')](1.0)[:-1])
			else: ax.set_yticks([])
			if ss==0:
				ax.set_title({'PI2P':'$$\mathrm{{PIP}_{2}}$$'}.get(resname,resname),
					fontsize=plotspec['fs_title'])
			ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
			resids = postdat[sn]['resids']
			if do_crazy_colors:
				image_data = np.array([mpl.cm.__dict__[
					residue_colors.get(resnames[rnum],'Greys')](
						postpostdat[sn][resname][rnum]/ceiling).T[:-1].T 
					for rnum in np.arange(postpostdat[sn][resname].shape[0])])
			else: image_data = postpostdat[sn][resname]
			duration = (data[sn]['data']['times'].max()-data[sn]['data']['times'].min())*10**-3
			ax.imshow(image_data,origin='lower',interpolation='nearest',
				extent=[0,duration,0,len(resids)],
				aspect=float(duration)/len(resids),
				vmax=ceiling,vmin=0,cmap=mpl.cm.__dict__[cmap])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')

	#---title on the legend
	title = 'contacts\n'+r'$\mathrm{r \leq %.1f \AA}$'%(cutoff)
	patches,labels = zip(*[(mpl.patches.Rectangle((0,0),1.0,1.0,
		fc=mpl.cm.__dict__[residue_type_colors[kind]](plotspec['legend_color_strength'])),kind) 
		for kind in residue_type_colors])
	legend = axes[-1][0].legend(patches,labels,loc='upper left',fontsize=plotspec['fs_legend'],
		bbox_to_anchor=(1.05,0.0,1.,1.),title=title,shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	plt.setp(legend.get_title(),fontsize=plotspec['fs_legend_title'])
	picturesave('fig.contacts_occupancy.residue_lipid_time.cutoff_%.1f'%cutoff,work.plotdir,
		backup=False,version=True,meta={},extras=[legend])

#---block: prepare contact map data
if 'alt' in routine:

	sn = sns[0]
	dat = data[sn]['data']
	bonds,obs = dat['bonds'],dat['observations']
	pass