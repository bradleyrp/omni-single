#!/usr/bin/env python

"""
PLOT hydrogen bonding
Tested by rpb on ptdins on dark on 2017.07.15 at the terminal and in jupyter.
Further development through 2017.08.05 to add subsets of the results, at which point this code is customized
specifically for the PtdIns project. Recall that several weeks ago it was also running on protein-lipid
systems, however 
"""

import copy
import brewer2mpl
if is_live: 
	from ipywidgets import *
	from IPython.display import display

#---block: SETTINGS
show_static,press = True,not is_live
#---aesthetics
mpl.rcParams['hatch.linewidth'] = 1.5
#---hatch color was controllable here in mpl v2.0.0 but they reverted
mpl.rcParams['hatch.color'] = 'w'

#---block: IMPORT THE DATA
if 'data' not in globals(): 
	sns,(data,calc) = work.sns(),plotload(plotname,work)
	data_salt,calc_salt = plotload('salt_bridges',work)
	#---set the plot targets in the YAML for multiuse scripts
	plotspecs = work.plots[plotname]['specs']
	#---subselect interesting bars from the master set
	interesting = plotspecs.get('interesting',{})
	press_routine = plotspecs.get('press_routine',['summary','interesting'])
	collections = plotspecs.get('collection_sets')
	ion_order = plotspecs.get('ion_order')
	legend_mapper = plotspecs.get('legend_mapper')
	relevant_resnames = plotspecs.get('relevant_resnames',None)
	bar_formats_style = plotspecs.get('bar_formats_style','candy')
	legend_maker_specs = plotspecs.get('legend_maker_specs',{})
	figsize_all = plotspecs.get('figsize_all',(16,9))
	sns_explicit_color_names = plotspecs.get('sns_explicit_color_names',{})
	nmols_recount = plotspecs.get('nmols_recount',None)

#---block: UTILITY FUNCTIONS
def reorder_by_ion(sns,ion_order,identity=True):
	"""reorder simulations by ion name"""
	#---identity function means that the order follows the order of the collection
	if identity: return sns
	#---otherwise we order by ion. probably better to use the collection order for more control
	return sorted(sns,key=lambda x:ion_order.index(work.meta[x]['cation']))

def dimmer(bar,child):
	"""Dim a bar in the interactive plots."""
	bar.set_edgecolor(bar.get_facecolor())
	bar.set_facecolor('w')
	child.set_hatch(None)

def count_hydrogen_bonds_by_sn(bonds,obs,nmols):
	"""!!!"""
	#---filter out intralipid hydrogen bonds
	#---some simulations have no salt bridges at all
	try: resids_d = bonds[:,1].astype(int)
	except: return {'result':{},'result_err':{}}
	resids_a = bonds[:,4].astype(int)
	inter_inds = np.where(resids_d!=resids_a)[0]
	inter = bonds[inter_inds]
	#---catalog bonding pairs by residue
	resnames_combos = np.transpose((inter[:,0],inter[:,3]))
	if relevant_resnames: 
		resnames_combos = resnames_combos[np.where(np.all([
			np.in1d(resnames_combos[:,i],relevant_resnames) for i in range(2)],axis=0))]
	inds,counts = uniquify(resnames_combos)
	pairs = resnames_combos[inds]
	#---discard sol
	lipid_only = np.where(np.all(pairs!='SOL',axis=1))[0]
	#---cut cholesterol in half because it is in both leaflets and POPC does not make hbonds
	nmols_leaflet = 400 
	if nmols_recount:
		nmols_recount_f = eval(nmols_recount)
		nmols = nmols_recount_f(nmols)
	#---get the proportions relative to combos
	#---! CHECK THE MEANING OF THIS NUMBER, COMBINATORICALLY-WISE
	#---subsample the obs matrix (frames by combos) to pick out each unique resname combo
	#---summing over the first axis adds up all unique instances of a particular combo
	counts_by_pair = [obs[:,np.where(resnames_combos==p)[0]].sum(axis=1) for p in pairs[lipid_only]]
	#---counts are normalized by the number of each species (note that 400 choose 2 is 79800)
	norm_combos = lambda i,j: (nmols[i]*nmols[j])
	props_mean = np.array([float(counts_by_pair[k].mean())/norm_combos(i,j)
		for k,(i,j) in enumerate(pairs[lipid_only])])
	props_std = np.array([float(counts_by_pair[k].std())/norm_combos(i,j)
		for k,(i,j) in enumerate(pairs[lipid_only])])
	#---convert pairs to PtdIns alias
	pip2_alias = lambda x: 'PtdIns' if x in work.vars['selectors']['resnames_PIP2'] else x
	aliased_names = np.array([[pip2_alias(i) for i in j] for j in pairs[lipid_only]])
	result = dict(zip([tuple(i) for i in aliased_names],props_mean))
	result_err = dict(zip([tuple(i) for i in aliased_names],props_std))
	return {'result':result,'result_err':result_err}

def count_hydrogen_bonds(data,bonds_key='bonds'):
	"""
	Compute the hydrogen bond counts.
	You can send salt bridge data which only has extra columns with the cation information.
	POSSIBLY DEPRECATED BY count_hydrogen_bonds_by_sn?
	"""
	post = {}
	sns = work.sns()
	for sn in sns:
		dat = data[sn]['data']
		bonds,obs = dat[bonds_key],dat['observations']
		#---filter out intralipid hydrogen bonds
		#---some simulations have no salt bridges at all
		try: resids_d = bonds[:,1].astype(int)
		except: 
			post[sn] = {'result':{},'result_err':{}}
			continue
		resids_a = bonds[:,4].astype(int)
		inter_inds = np.where(resids_d!=resids_a)[0]
		inter = bonds[inter_inds]
		#---catalog bonding pairs by residue
		resnames_combos = np.transpose((inter[:,0],inter[:,3]))
		inds,counts = uniquify(resnames_combos)
		pairs = resnames_combos[inds]
		#---discard sol
		lipid_only = np.where(np.all(pairs!='SOL',axis=1))[0]
		nmols = dict(zip(dat['resnames'],dat['nmols']))
		#---cut cholesterol in half because it is in both leaflets and POPC does not make hbonds
		nmols_leaflet = 400 
		if 'CHL1' in nmols: nmols['CHL1'] /= 2.0
		#---get the proportions relative to combos
		#---! CHECK THE MEANING OF THIS NUMBER, COMBINATORICALLY-WISE
		#---subsample the obs matrix (frames by combos) to pick out each unique resname combo
		#---summing over the first axis adds up all unique instances of a particular combo
		counts_by_pair = [obs[:,np.where(resnames_combos==p)[0]].sum(axis=1) for p in pairs[lipid_only]]
		#---counts are normalized by the number of each species (note that 400 choose 2 is 79800)
		norm_combos = lambda i,j: (nmols[i]*nmols[j])
		props_mean = np.array([float(counts_by_pair[k].mean())/norm_combos(i,j)
			for k,(i,j) in enumerate(pairs[lipid_only])])
		props_std = np.array([float(counts_by_pair[k].std())/norm_combos(i,j)
			for k,(i,j) in enumerate(pairs[lipid_only])])
		#---convert pairs to PtdIns alias
		pip2_alias = lambda x: 'PtdIns' if x in work.vars['selectors']['resnames_PIP2'] else x
		aliased_names = np.array([[pip2_alias(i) for i in j] for j in pairs[lipid_only]])
		result = dict(zip([tuple(i) for i in aliased_names],props_mean))
		result_err = dict(zip([tuple(i) for i in aliased_names],props_std))
		post[sn] = {'result':result,'result_err':result_err}
	return post

def make_bar_formats(sns,style='candy'):
	"""Make bar formats, most recently, candy-colored bars with hatches."""
	colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
		for key,val in {
		'red':0,'blue':1,'green':2,'purple':3,'orange':4,
		'yellow':5,'brown':6,'pink':7,'grey':8,}.items()])
	colors['pink'] = mpl.colors.ColorConverter().to_rgb("#f1948a")
	colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")
	#---combine blue and green to denote the dilute Na,Cal simulation
	bgw = 0.3
	colors['bluegreen'] = np.average((colors['blue'],colors['green']),weights=[1-bgw,bgw],axis=0)
	bgw2 = 0.6
	colors['bluegreen2'] = np.average((colors['blue'],colors['green']),weights=[1-bgw2,bgw2],axis=0)
	if style=='original':
		out = dict([(sn,{'c':colorize(work.meta[sn],comparison=comparison)}) for sn in sns])
	elif style=='candy':
		colors_ions = {'NA':'green','Na,Cal':'bluegreen','MG':'pink','Cal':'blue','K':'grey',}
		hatches_lipids = {'PI2P':'//','P35P':'-','PIPU':'xx','PIPP':'++','SAPI':''}
		out = dict([(sn,{
			'c':colors[colors_ions[work.meta[sn]['cation']]],
			'hatch':hatches_lipids[work.meta[sn]['ptdins_resname']]}) for sn in sns])
	elif style=='actinlink':
		out = dict([(sn,{
			#---the actinlink plots have custom colors
			'c':colors[sns_explicit_color_names[sn]],
			'hatch':'//' if work.meta[sn].get('cholesterol',False) else ''}) for sn in sns])
	else: raise Exception('no bar style: %s'%style)
	for k,v in out.items(): v.update(edgecolor=v['c'])
	return out

def legend_maker_stylized(ax,sns_this,title=None,ncol=1,
	bbox=(1.05,0.0,1.,1.),loc='upper left',fs=16,extra_legends=None,
	comparison_spec=None,lipid_resnames=None,sns_explicit=None):
	global bar_formats,comparison
	#---custom items in the legend based on the comparison
	if comparison_spec:
		ion_names,ptdins_names = [comparison_spec[i] for i in ['ion_names','ptdins_names']]
		#---custom simulations must come in through the comparison_spec in the YAML file
		sns_explicit = comparison_spec.get('sns_explicit',None)
	else:
		if comparison not in legend_mapper:
			raise Exception('you must add this comparison (%s) to legend_mapper in the plot specs') 
		ion_names,ptdins_names = [legend_mapper[comparison][k] for k in ['ion_names','ptdins_names']]
	#---we can also add lipid types to the bar plots for the subset plots
	lipid_resnames = [] if not lipid_resnames else lipid_resnames
	rectangle_specs,patches,labels = {},[],[]
	#---loop over relevant cations, PIP2 types, and lipids so they are included in the legend
	for name,group in [('cation',ion_names),('ptdins_resname',ptdins_names),
		('lipid_resname',lipid_resnames)]:
		for item in group:
			if name!='lipid_resname':
				try: sn = next(sn for sn in sns_this if work.meta[sn][name]==item)
				except: 
					import ipdb;ipdb.set_trace()
					raise Exception('cannot get a simulation for %s,%s'%(name,item))
			#---simulation does not matter if we are adding lipids to the legend
			else: sn = sns_this[0]
			if sn not in sns_this: continue
			rectangle_specs[sn] = {}
			if name=='cation':
				rectangle_specs[sn]['fc'] = bar_formats[sn]['edgecolor']
				rectangle_specs[sn]['lw'] = 0
				patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]))
				labels.append(work.meta[sn]['ion_label'])
			elif name=='ptdins_resname':
				patches.append(mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=1.0,fc='w',lw=3,ec='k',
					hatch=bar_formats[sn]['hatch']*(1 if work.meta[sn]['ptdins_resname']=='PI2P' else 2)))
				labels.append(work.meta[sn]['ptdins_label'])
			elif name=='lipid_resname':
				rectangle_specs[sn]['fc'] = colorize(work.meta[sn],resname=item)
				patches.append(mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=1.0,lw=3,ec='k',
					**rectangle_specs[sn]))
				labels.append(item)
	#---additional legend items for explicit simulations that come from sns_explicit
	if sns_explicit:
		for sn in sns_explicit:
			if sn not in rectangle_specs: rectangle_specs[sn] = {}
			rectangle_specs[sn]['ec'] = 'k'
			rectangle_specs[sn]['fc'] = bar_formats[sn]['edgecolor']
			rectangle_specs[sn]['lw'] = 0
			rectangle_specs[sn]['hatch'] = bar_formats[sn].get('hatch','')
			patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]))
			labels.append(work.meta[sn].get('label',sn))
	#---special annotation for cholesterol
	if bar_formats_style=='actinlink':
		for is_chol in [True,False]:
			rectangle_specs[sn]['ec'] = 'k'
			rectangle_specs[sn]['fc'] = 'w'
			rectangle_specs[sn]['lw'] = 2
			rectangle_specs[sn]['hatch'] = {True:'//',False:''}[is_chol]
			patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]))
			labels.append('%s cholesterol'%{True:'with',False:'no'}[is_chol])
	if not patches: patches,labels = [mpl.patches.Rectangle((0,0),1.0,1.0,fc='w')],['empty']
	legend = ax.legend(patches,labels,loc=loc,fontsize=fs,
		ncol=ncol,title=title,bbox_to_anchor=bbox,labelspacing=1.2,
		handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	return legend,patches

def hbonds_bardat(bar_formats,sns,post,format_name=None,ion_order=None):
	"""Generate hydrogen bonding counts bar data."""
	format_name = format_name if format_name else 'color_schemer_original'
	#---! move hatch choices to art_ptdins.py
	hatches_lipids = {'PI2P':'//','P35P':'-','PIPU':'xx','PIPP':'++','SAPI':'',}	
	if ion_order: sns = reorder_by_ion(sns,ion_order)
	#---empty plot uses a random key and zeros it
	empty = not sns
	if not empty: 
		keys = sorted(set(list([tuple(j) for j in 
			np.concatenate([m for m in [post[sn]['result'].keys() for sn in sns] if m])])))
	else: keys = post.keys()[:1]
	bardat,barspec,extras = [],[],[]
	#---sorted gives us a nice order with PtdIns at the end but this might be worth hardcoding
	donors = list(np.unique([key[0] for key in keys]))
	#---group by the first key (the donor)
	for donor in donors:
		group_donor,group_barspec,names = [],[],[]
		#---group by pairs
		for key in [k for k in keys if k[0]==donor]:
			#---group by simulation
			try: group_donor.append([((key,sn),
				(post[sn]['result'].get(key,0.0),post[sn]['result_err'].get(key,0.0)) 
				if not empty else 0.0) for sn in sns])
			except:
				import ipdb;ipdb.set_trace()
			#---include the label with the bar format
			group_barspec.append([dict(tl='%s\m%s'%key,**bar_formats[sn]) for sn in sns])
			new_names = [(key,sn) for sn in sns]
			names.extend(new_names)
			extras.append(('xlabel',{'names':new_names,'label':key[1],'rotation':0,'fs':20}))
		#---pairs are flanked by dividers
		extras.append(('divider',{
			'names':names,'label':key[0],'fs':26,'bbox':dict(boxstyle="round4",fc="w",lw=3)}))
		bardat.append(group_donor)
		barspec.append(group_barspec)
	return dict(bardat=bardat,barspec=barspec,extras=extras)

def postprep(sns):
	"""
	!!!
	"""
	sns_names = dict([(sn,work.meta[sn]['ptdins_resname']+' and '+
		work.meta[sn]['cation']) for sn in sns])
	namer = lambda ((a,b),c): ('%s-%s'%(a,b))
	bar_formats = make_bar_formats(sns,style=bar_formats_style)
	#---premake the data for the outline style
	posts = dict([(key,{}) for key in ['hbonds','salt']])
	for key,(bonds_key,data_name) in [('hbonds',('bonds','data')),('salt',('bonds_salt','data_salt'))]:
		this_data = globals()[data_name]
		posts[key] = dict([(sn,
			count_hydrogen_bonds_by_sn(
				bonds=this_data[sn]['data'][bonds_key],
				obs=this_data[sn]['data']['observations'],
				nmols=dict(zip(this_data[sn]['data']['resnames'],this_data[sn]['data']['nmols'])))) 
			for sn in sns])
	#---merge the two calculations
	posts['both'] = {}
	for sn in sns:
		posts['both'][sn] = {}
		for measure in ['result','result_err']:
			posts['both'][sn][measure] = {}
			pairs = list(set(posts['hbonds'][sn][measure].keys()+posts['salt'][sn][measure].keys()))
			for pair in pairs:
				#---sum the means
				if measure=='result':
					val = posts['hbonds'][sn][measure].get(pair,0.0)+posts['salt'][sn][measure].get(pair,0.0)
				elif measure=='result_err':
					val = np.sqrt(posts['hbonds'][sn][measure].get(pair,0.0)**2+
						posts['salt'][sn][measure].get(pair,0.0)**2)
				else: raise Exception('measure is %s'%measure)
				posts['both'][sn][measure][pair] = val
	#---filter out POPC because we don't care about outer leaflet lipids
	for upkey in ['both','salt']:
		for sn in sns:
			for measure in ['result','result_err']:
				pairs = list(set(posts[upkey][sn][measure].keys()+posts['salt'][sn][measure].keys()))
				pairs = [p for p in pairs if 'POPC' not in p]
				posts[upkey][sn][measure] = dict([(p,posts[upkey][sn][measure][p]) for p in pairs])
	#---select the plotting target
	post,title = posts['both'],'hydrogen bonding + salt bridges'
	#---prepare the data
	bardat = hbonds_bardat(sns=sns,post=post,bar_formats=bar_formats)
	bars_premade = BarGrouper(figsize=figsize_all,spacers={0:0,1:1.0,2:2.0},
		dimmer=dimmer,namer=namer,show_xticks=False,**bardat)
	return dict(namer=namer,bar_formats=bar_formats,posts=posts,sns_names=sns_names)

#---block: plot function and static plot
def hbonds_plotter(bars=None,checks=None,**kwargs):
	"""Plot hydrogen bond counts."""
	title = kwargs.get('title','hydrogen bonds')
	#---prepare the data
	sns = [k for k,v in checks.items() if v]
	ion_order = ['K','NA','Na,Cal','MG','Cal']
	sns = reorder_by_ion(sns,ion_order)
	#---incoming bars creates the outlined-if-unchecked style plot
	#---previously used "outline" method which made bar outlines less transparent when they were selected
	if bars:
		for sn_checked in reorder_by_ion(checks.keys(),ion_order):
			if checks[sn_checked]: bars.on(*[(key,sn) for (key,sn) in bars.names if sn==sn_checked])
			else: bars.off(*[(key,sn) for (key,sn) in bars.names if sn==sn_checked])
	#---no incoming bars tells us to recompute the bar object for the bars that we want
	else:
		global post,bar_formats
		keys = reorder_by_ion([k for k,v in checks.items() if v],ion_order)
		empty = not keys
		bardat = hbonds_bardat(sns=checks.keys()[:1] if empty else keys,post=post,bar_formats=bar_formats)
		bars = BarGrouper(figsize=figsize_all,spacers={0:0,1:1.0,2:2.0},
			dimmer=dimmer,namer=namer,show_xticks=False,empty=empty,**bardat)
	#---additional 
	bars.plot(wait=True)
	bars.ax.tick_params(axis='both',which='major',labelsize=16)
	bars.ax.set_ylabel(r'$\frac{\langle {N}_{bonds} \rangle}{{N}_{donor}{N}_{acceptor}}$',fontsize=26)
	bars.ax.yaxis.grid(which="minor",color='k',linewidth=0.5,alpha=0.5,zorder=0)
	bars.ax.tick_params(axis=u'both',which=u'both',length=0)
	bars.ax.set_title(title,fontsize=40)
	bars.ax.set_ylim((0,bars.ax.get_ylim()[1]))
	legend_kwargs = dict(sns_this=sns)
	#---overrides for the legendmaker
	legend_kwargs.update(**legend_maker_specs)
	legend,patches = legend_maker_stylized(bars.ax,
		sns_this=legend_kwargs.pop('sns'),comparison_spec=legend_kwargs)
	return legend

###---PLOTTING

#---block: static plots and plot press
#---show a static version of the plot
if show_static and is_live:
	checks = dict([(sn,widgets.ToggleButton(value=False,description=sns_names[sn])) for sn in sns])
	hbonds_plotter(bars_premade,extras=bars_premade.extras,checks=checks)
	plt.show()
#---batch plot the summary of all the bonds
if press and 'summary' in press_routine:
	#---which figures to make
	press_specs = {
		'bonding.salt':{'posts_subkey':'salt','title':'salt bridges'},
		'bonding.hbonds':{'posts_subkey':'hbonds','title':'hydrogen bonds'},
		'bonding.both_hbonds_salt':{
			'posts_subkey':'both','title':'hydrogen bonds + salt bridges'},}
	#---which simulations to compare
	sns_sets = dict([(comparison,reorder_by_ion(work.specs['collections'][comparison],ion_order))
		for comparison in collections])
	#---batch of figures with a loop over figure types and simulation sets
	for comparison,sns in sns_sets.items():
		#---put the postprep data into globals
		globals().update(postprep(sns))
		for tag,spec in press_specs.items():
			bars_premade = BarGrouper(figsize=figsize_all,spacers={0:0,1:1.0,2:2.0},
				dimmer=dimmer,namer=namer,show_xticks=False,
				**hbonds_bardat(sns=sns,post=posts[spec['posts_subkey']],bar_formats=bar_formats))
			legend = hbonds_plotter(bars_premade,title=spec['title'],checks=dict([(sn,True) for sn in sns]))
			picturesave('fig.%s.%s'%(tag,comparison),work.plotdir,
				backup=False,version=True,meta={},extras=[legend])
#---batch plot the "interesting" subsets of the main plots
if press and 'interesting' in press_routine:
	#---loop over interesting panel plots, one per figure
	for iname,ispec in interesting.items():
		comparison = ispec['comparison']
		bond_type = ispec['bond_type']
		#---get simulations for this comparison
		sns = work.specs['collections'][comparison]
		globals().update(postprep(sns))
		#---whittle sns by other features
		if ispec.get('cation',False): sns = [sn for sn in sns if 
			ispec['cation']==work.meta[sn]['cation']]
		#---get valid resnames
		pairs = ispec['pairs']
		#---organize the keys first by pairs then by simulation, for the bar grouper
		keys = [[(sn,pair) for sn in sns] for pair in pairs]
		bardat = dict(
			bardat=[[(((r1,r2),sn),tuple([posts[bond_type][sn][r][(r1,r2)] 
				for r in ['result','result_err']])) for sn,(r1,r2) in key] for key in keys],
			barspec=[[dict(tl='%s\m%s'%(r1,r2),**dict(bar_formats[sn],
				c=colorize(work.meta[sn],resname=r1))) for sn,(r1,r2) in key] for key in keys],
			#---the label is the first in the pair since all comparisons in this block are to PtdIns
			extras=[('xlabel',{'label':key[0][1][0],
				'names':[((r1,r2),sn) for sn,(r1,r2) in key]}) for key in keys])
		bars = BarGrouper(figsize=(4,6),show_xticks=False,**bardat)
		bars.plot(wait=True)		
		bars.ax.tick_params(axis=u'both',which=u'both',length=0)
		bars.ax.set_title('%s %s with %s'%({'pip2':'PtdIns','chl1':'cholesterol'}[ispec['pairname']],
			{'both':'hydrogen bonds + salt bridges','hbonds':'hydrogen bonds'}[bond_type],
			work.meta[sns[0]]['ion_label']))
		bars.ax.set_ylabel(r'$\frac{\langle {N}_{bonds} \rangle}{{N}_{donor}{N}_{acceptor}}$',fontsize=14)
		legend,patches = legend_maker_stylized(bars.ax,fs=10,
			sns_this=list(set([sn for sn,pair in key for key in keys])),
			comparison_spec={'ion_names':[],'ptdins_names':['PI2P','P35P']},
			lipid_resnames=sorted(set([r1 for key in keys for sn,(r1,r2) in key])))
		picturesave('fig.bonding_subset.%s.%s'%(bond_type,ispec['pairname']),work.plotdir,
			backup=False,version=True,meta={'interesting':ispec},extras=[legend])
		#---! still needs: legend, ylabel

#---block: interactive plot
if is_live:
	#---! sometimes this throws an annoying "widget javascript not enabled" warning
	#---! ...which you can eliminate by running "jupyter nbextension enable --py widgetsnbextension"
	#---! ...followed by a full shutdown/run of the factory
	def hbonds_plotter_interactive(**checks): 
		hbonds_plotter(bars={'outline':bars_premade,'reveal':None}[live_plot_style],checks=checks)
		plt.show()
	checks = dict([(sn,widgets.ToggleButton(value=sn==sns[0],description=sns_names[sn])) for sn in sns])
	widget = interactive(hbonds_plotter_interactive,checks=checks)
	ws = [widget.children[[w.description for w in widget.children[:-1]].index(sns_names[i])] for i in sns]
	display(HBox(ws),widget.children[-1])
	print('use the toggle buttons above to select different simulations for comparison')
	widget.update()
