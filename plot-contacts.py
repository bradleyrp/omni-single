#!/usr/bin/env python

"""
"""

routine = ['contact_map','contact_by_lipid','contact_by_lipid_by_restype',
	'lipid_capture'][:1]
mode_list = ['explicit','reduced']
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

#---block: settings and definitions
residue_types = {'ARG':'basic','HIS':'basic','LYS':'basic',
	'ASP':'acidic','GLU':'acidic','SER':'polar','THR':'polar','ASN':'polar','GLN':'polar',
	'ALA':'hydrophobic','VAL':'hydrophobic','ILE':'hydrophobic','LEU':'hydrophobic',
		'MET':'hydrophobic','PHE':'hydrophobic','TYR':'hydrophobic','TRP':'hydrophobic',
	'CYS':'special','SEC':'special','GLY':'special','PRO':'special'}
residue_type_colors = {'basic':'Blues','acidic':'Reds','hydrophobic':'Greens',
	'polar':'Purples','special':'Oranges'}

#---common colors
import brewer2mpl
colors = {'DOPC':'blue','DOPS':'red','POP2':'magenta','PI2P':'magenta',
	'all lipids':'black','DOPE':(0.0,1.0,1.0)}
colors.update(**dict([(sn,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[sns.index(sn)]) 
	for sn in sns]))
lipid_label = lambda x: dict([(i,'$$\mathrm{{PIP}_{2}}$$') 
	for i in work.vars['selectors']['resnames_PIP2']]).get(x,x)

#---block: transform incoming data
if any([i in routine for i in ['contact_map','contact_by_lipid','contact_by_lipid_by_restype']]):

	postdat_super,postpostdat_super = {},{}
	mode_spec = [('explicit','counts_resid_resname'),('reduced','counts_resid_resname_singleton')]
	for mode,post_key in mode_spec:	
		#---assemble postdat from results
		postdat = dict([(sn,dict(residue_lipid_time=dict([((int(i[0]),i[1]),j) for i,j in 
			zip(data[sn]['data']['pairs_resid_resname'],data[sn]['data'][post_key])]))) 
			for sn in sns])
		postdat_super[mode] = postdat
		for sn in sns:
			lipid_resnames = np.unique(
				data[sn]['data']['bonds'][:,rowspec.index('target_resname')])
			resname_combos = [(r,np.array([r])) for r in lipid_resnames]+[
				('all lipids',np.array(lipid_resnames))]
			lipid_resnames = list(zip(*resname_combos))[0]
			postdat[sn]['lipid_resnames'] = lipid_resnames
			postdat[sn]['resids'] = data[sn]['data']['subject_residues_resids']
		#---further chompdown
		postpostdat_super[mode] = dict([(sn,
			dict([(resname,np.array([postdat[sn]['residue_lipid_time'][(r,resname)] 
			for r in postdat[sn]['resids']])) for resname in postdat[sn]['lipid_resnames']])) 
			for sn in sns])

#---block: transform incoming data
if 'contact_map' in routine:

	#---consider singletons or multiplicity
	for mode in mode_list:
		postdat,postpostdat = postdat_super[mode],postpostdat_super[mode]

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
			'legend_loc':'upper right','fs_legend':14,'legend_color_strength':0.5,'fs_legend_title':20,
			'binary_color_intensity':0.5}
		residue_colors = dict([(name,residue_type_colors[residue_types[name]]) for name in residue_types])
		ticks_font = mpl.font_manager.FontProperties(family='Latin Modern Mono',style='normal',
			size=14,weight='normal',stretch='normal')
		#---same residue comparisons
		lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
		if not len(lipid_resnames)==1: raise Exception('cannot use different resnames across simulations')
		else: lipid_resnames = lipid_resnames[0]
		ceiling = float(max([max([i.max() for i in postpostdat[sn].values()]) for sn in sns]))
		#---plot
		axes,fig = panelplot(
			layout={'out':{'hspace':0.05,'wspace':0.05,'grid':[1,len(lipid_resnames)]},
			'ins':[{'hspace':0.2,'wspace':0.05,'grid':[len(sns),1]} 
			for i in lipid_resnames]},figsize=figsize)
		for ss,sn in enumerate(sns):
			for rr,resname in enumerate(lipid_resnames):
				ax = axes[rr][ss]
				resnames = data[sn]['data']['subject_residues_resnames']
				if rr==0: 
					sn_title = '%s%s'%(work.meta[sn].get('label',re.sub('_','-',sn)),
						'' if not work.meta[sn].get('cholesterol',False) else '\n(cholesterol)')
					ax.set_ylabel(sn_title,fontsize=plotspec['fs_ylabel'])
					if work.plots[plotname].get('settings',{}).get('show_residue_names',True):
						#---never set the ticklabels before the ticks
						ax.set_yticks(np.arange(len(resnames))+0.5)
						ax.set_yticklabels(resnames)
						for label in ax.get_yticklabels():
							label.set_color(mpl.cm.__dict__[
								residue_colors.get(str(label._text),'Greys')](1.0)[:-1])
				else: ax.set_yticks([])
				if ss==0:
					ax.set_title(lipid_label(resname),
						fontsize=plotspec['fs_title'])
				resids = postdat[sn]['resids']
				if do_crazy_colors:
					image_data = np.array([mpl.cm.__dict__[
						residue_colors.get(resnames[rnum],'Greys')](
							postpostdat[sn][resname][rnum].astype(float)/ceiling).T[:-1].T 
						for rnum in np.arange(postpostdat[sn][resname].shape[0])])
				else: image_data = postpostdat[sn][resname]
				duration = (data[sn]['data']['times'].max()-data[sn]['data']['times'].min())*10**-3
				xtick_interval = work.plots[plotname].get('settings',{}).get('xtick_interval',100.0)
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
					bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
				cbar = plt.colorbar(im,cax=axins,orientation="vertical")

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
		picturesave('fig.contacts_occupancy.residue_lipid_time.%s.cutoff_%.1f'%(mode,cutoff),
			work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: transform incoming data
if 'contact_by_lipid' in routine:

	"""
	Timewise contacts by lipid types, absolute or normed, plotted by simulation or lipid type.
	"""

	ways = {
		'by_lipid':{'outer':'lipid','inner':'system'},
		'by_system':{'outer':'system','inner':'lipid'},
		'by_system.normed':{'outer':'system','inner':'lipid','normed':True},}
	plotspec = {'line_color_strength':0.65,'fs_xlabel':14,'fs_legend':12}

	#---consistent lipid resnames
	lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
	if len(lipid_resnames)!=1: raise Exception('only one composition allowed')
	else: lipid_resnames = list(lipid_resnames[0])

	"""
	debugging mode exlicit, by system, normed
		raw.sum(axis=1).astype(int)
			we have 500 frames and 500 residues which is confusing, but this sum is over frames
			so we are looking at residues
			and at this high cutoff, some of the residues in the middle are e.g. 20000, 30000 and other
			round numbers
			which suggests that they are permanently bound 
				and we are seeing a contact on each frame
				and possibly it is making contact with all of the atoms or something
	debugging mode reduced, by system, not normed
		the reduced counts tell you lipid-protein contacts with no multiplicity
		hence these curves tell us the number of residue-lipid contacts by type
		and if we divided by the total number of residues (500) and looked at the all lipids line
			then it would tell us the probability of any protein residue being bound to a lipid
		but since this is a composite score, we really should norm it by the number of lipids
	"""

	#---consider singletons or multiplicity
	for mode in mode_list:
		postdat,postpostdat = postdat_super[mode],postpostdat_super[mode]

		for plot_sub_name,plot_sub_spec in ways.items():

			#---plot
			outer,inner = [{'system':sns,'lipid':lipid_resnames}[plot_sub_spec[i]] 
				for i in ['outer','inner']]
			axes,fig = square_tiles(len(outer),figsize=(10,10))
			for onum,outside in enumerate(outer):
				ax = axes[onum]
				for inum,inside in enumerate(inner):
					sn,lipid = [{'outer':outside,'inner':inside}[j] for j in 
						[next(key for key,val in plot_sub_spec.items() if val==i) 
						for i in ['system','lipid']]]
					if plot_sub_name=='by_lipid':
						ax.set_title(lipid)
						label = work.meta[sn].get('label',sn)
					elif plot_sub_name in ['by_system','by_system.normed']:
						ax.set_title(work.meta[sn].get('label',sn))
						label = lipid
					raw = postpostdat[sn][lipid].astype(float)
					subject_resnames = data[sn]['data']['subject_residues_resnames']
					residue_kinds = np.array([residue_types[j] for j in subject_resnames])
					kinds = np.sort(np.unique(residue_kinds))
					duration = (data[sn]['data']['times'].max()-data[sn]['data']['times'].min())*10**-3
					xtick_interval = work.plots[plotname].get('settings',{}).get('xtick_interval',100.0)
					ax.set_xticks(np.arange(xtick_interval,duration,xtick_interval))
					ax.set_xlabel('time (ns)',fontsize=plotspec['fs_xlabel'])
					ys = raw.sum(axis=0)
					if plot_sub_spec.get('normed',False):
						ys = ys/((data[sn]['data']['targets_residues_resnames']==lipid).sum() 
							if lipid!='all lipids' else float(
							#---only including lipids that show up in the main list
							np.in1d(data[sn]['data']['targets_residues_resnames'],lipid_resnames).sum()))
					ax.plot(ys,color=colors[inside],label=label)
					#---! this was pretty not good
					if False:
						for kind in kinds:
							color = mpl.cm.__dict__[residue_type_colors[kind]](
								plotspec['line_color_strength'])
							ax.plot(raw[np.where(residue_kinds==kind)].sum(axis=0),color=color)
					if False and mode=='reduced' and plot_sub_name=='by_system.normed' and inside=='POP2':
						import ipdb;ipdb.set_trace()

			legend = ax.legend(loc='upper left',fontsize=plotspec['fs_legend'],
				bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True)
			picturesave('fig.contact_by_lipid.%s.%s.cutoff_%.1f'%(mode,plot_sub_name,cutoff),
				work.plotdir,backup=False,version=True,meta={},extras=[legend])

#---block: transform incoming data
if 'contact_by_lipid_by_restype' in routine:

	"""
	Heat maps of the counts of contacts between residue types and lipid types. 
	Also generates a version which is normed by the quantitity of both partners.
	"""

	#---consider singletons or multiplicity
	for mode in mode_list:
		postdat,postpostdat = postdat_super[mode],postpostdat_super[mode]

		#---settings
		for normed in [True,False]:

			#---consistent lipid resnames
			lipid_resnames = list(set([tuple(postdat[sn]['lipid_resnames']) for sn in sns]))
			if len(lipid_resnames)!=1: raise Exception('only one composition allowed')
			else: lipid_resnames = list(lipid_resnames[0])

			#---precounts
			prepped = {}
			for sn in sns:
				subject_resnames = data[sn]['data']['subject_residues_resnames']
				residue_kinds = np.array([residue_types[j] for j in subject_resnames])
				kinds = np.sort(np.unique(residue_kinds))
				prepped[sn] = dict(kinds=kinds,
					raw=np.array([[postpostdat[sn][lname][np.where(residue_kinds==kind)].sum(axis=0).mean()/
					(np.sum(residue_kinds==kind) if normed else 1.0)/
					(((data[sn]['data']['targets_residues_resnames']==lname).sum() if lname!='all lipids' 
						else float(len(data[sn]['data']['targets_residues_resnames'])))
						if normed else 1.0) 
					for lname in lipid_resnames] for kind in kinds]))

			#---plot
			axes,fig = square_tiles(len(sns),figsize=(8,8),wspace=0.5)
			vmax = max([prepped[sn]['raw'].max() for sn in sns])
			for snum,sn in enumerate(sns):
				ax = axes[snum]
				ax.set_title(work.meta[sn].get('label',sn))
				raw,kinds = [prepped[sn][i] for i in ['raw','kinds']]
				im = ax.imshow(raw.T,origin='lower',interpolation='nearest',
					cmap=mpl.cm.__dict__['Greys'],vmax=vmax)
				ax.set_xticks(np.arange(len(kinds)))
				ax.set_xticklabels(kinds,rotation=90)
				ax.set_yticks(np.arange(len(lipid_resnames)))
				ax.set_yticklabels([work.vars['names']['proper_residue_names_long'].get(l,l) 
					for l in lipid_resnames])
				ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
				ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
				if snum==len(sns)-1:
					axins = inset_axes(ax,width="5%",height="100%",loc=3,
						bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
					cbar = plt.colorbar(im,cax=axins,orientation="vertical")
			picturesave('fig.contact_by_lipid_by_restype.%s.%s.cutoff_%.1f'%(
				mode,'normed' if normed else 'abs',cutoff),work.plotdir,
				backup=False,version=True,meta={},extras=[])

#---block: transform incoming data
if 'lipid_capture' in routine:

	"""
	The above methods consider the total number of bonds or bound lipids (explicit vs reduced) for either 
	specific protein residues (contact maps) or for all residues (contact_by_lipid) or for all residues by
	the protein residue type (contact_by_lipid_by_restype). In this section 
	"""

	specs = {
		'fig.contact_time_distributions.normed':{'normed':True},
		'fig.contact_time_distributions.abs':{'normed':False}}

	def duration_of_stay(series):
		"""
		"""
		zeros = np.concatenate(([False],series==1,[False]))
		diffs = np.diff(zeros*1.0)
		starts, = np.where(diffs==1)
		ends, = np.where(diffs==-1)
		runs = ends - starts
		counts,run = np.histogram(runs,bins=np.arange(1,runs.max()+2))
		return counts,run

	plotspec = {'fs_legend':14,'figsize':(10,6)}

	for specname,spec in specs.items():

		axes,fig = panelplot(
			layout={'out':{'grid':[1,1]},'ins':[
			{'hspace':0.2,'wspace':0.2,'grid':[1,len(sns)]}]},figsize=plotspec['figsize'])

		for snum,sn in enumerate(sns):
			ax = axes[snum]
			dat = data[sn]['data']
			bonds,obs = [dat[i] for i in ['bonds','observations']]
			lipids,idx = np.unique(bonds[:,rowspec.index('target_resid')],return_index=True)
			lipid_resnames = bonds[:,rowspec.index('target_resname')][idx]
			resnames = np.unique(lipid_resnames)
			resnames_index = dict([(i,ii) for ii,i in enumerate(resnames)])
			#---get the rows in the bond list that apply to each distinct lipid that participates in a bond
			rows_by_lipids = [np.where(bonds[:,rowspec.index('target_resid')]==l)[0] for l in lipids]
			traj = np.array([(obs.T[r].sum(axis=0)>0.)*1 for r in rows_by_lipids])
			intermediate = []
			for this_traj,resname in zip(traj,lipid_resnames):
				dist,duration = duration_of_stay(this_traj)
				intermediate.append(dict(dist=dist,duration=duration,resname=resname))
			durations = [i['duration'] for i in intermediate]
			times = np.arange(1,max([i.max() for i in durations])+1)
			counts = [np.zeros((len(times))) for kind in resnames]
			for result in intermediate:
				#---! sloppy
				resname,dist,duration = [result[j] for j in ['resname','dist','duration']]
				counts[resnames_index[resname]][duration[:-1]-1] += dist
			for resname,count in zip(resnames,counts):
				ax.plot(times[:len(count)][np.where(count>0)],count[np.where(count>0)]/count.max(),
					marker='.',lw=0,label=lipid_label(resname),color=colors[resname])
				unbroken = range(0,np.where(count==0)[0][0])
				ax.plot(times[unbroken],count[unbroken]/(count.max() if spec['normed'] else 1.0),
					color=colors[resname])
			if snum==len(sns)-1:
				legend = ax.legend(loc='upper left',fontsize=plotspec['fs_legend'],
					bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True,
					title=r'residence time $\mathrm{(\tau)}$'+'\n'+r'$\mathrm{r\leq%.1f\AA}$'%cutoff)	
			ax.set_title(work.meta[sn].get('label',sn))
			ax.set_yscale('log')
			ax.set_xscale('log')
			ax.set_xlabel(r'$\mathrm{\tau}$ (ns)')
			if snum==0: ax.set_ylabel(r'$\mathrm{{N(\tau)}/{N}_{total}}$')

		picturesave('%s.cutoff_%.1f'%(specname,cutoff),work.plotdir,
			backup=False,version=True,meta={},extras=[legend])
