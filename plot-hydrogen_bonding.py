#!/usr/bin/env python

"""
PLOT hydrogen bonding
"""

import copy
import brewer2mpl
if is_live: 
	from ipywidgets import *
	from IPython.display import display

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns,(data,calc) = work.sns(),plotload(plotname,work)
	data_salt,calc_salt = plotload('salt_bridges',work)

#---block: prepare data for plotting
def reorder_by_ion(sns,ion_order,identity=True):
	"""reorder simulations by ion name"""
	#---identity function means that the order follows the order of the collection
	if identity: return sns
	#---otherwise we order by ion. probably better to use the collection order for more control
	return sorted(sns,key=lambda x:ion_order.index(work.meta[x]['cation']))

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
	inds,counts = uniquify(resnames_combos)
	pairs = resnames_combos[inds]
	#---discard sol
	lipid_only = np.where(np.all(pairs!='SOL',axis=1))[0]
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
	return {'result':result,'result_err':result_err}

def count_hydrogen_bonds(data,bonds_key='bonds'):
	"""
	Compute the hydrogen bond counts.
	You can send salt bridge data which only has extra columns with the cation information.
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

def make_bar_formats(sns,style):
	"""Make bar formats, most recently, candy-colored bars with hatches."""
	if style=='original':
		out = dict([(sn,{'c':colorize(work.meta[sn],comparison=comparison)}) for sn in sns])
	elif style=='candy':
		colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
			for key,val in {
			'red':0,'blue':1,'green':2,'purple':3,'orange':4,
			'yellow':5,'brown':6,'pink':7,'grey':8,}.items()])
		colors['pink'] = mpl.colors.ColorConverter().to_rgb("#f1948a")
		colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")
		#---combine blue and green to denote the dilute Na,Cal simulation
		bgw = 0.3
		colors['bluegreen'] = np.average((colors['blue'],colors['green']),weights=[1-bgw,bgw],axis=0)
		colors_ions = {'NA':'green','Na,Cal':'bluegreen','MG':'pink','Cal':'blue','K':'grey',}
		hatches_lipids = {'PI2P':'//','P35P':'-','PIPU':'xx','PIPP':'++','SAPI':''}
		out = dict([(sn,{
			'c':colors[colors_ions[work.meta[sn]['cation']]],
			'hatch':hatches_lipids[work.meta[sn]['ptdins_resname']]}) for sn in sns])
	else: raise Exception('no bar style: %s'%style)
	[v.update(edgecolor=v['c']) for k,v in out.items()]
	return out

def legend_maker_stylized(ax,sns_this,title=None,ncol=1,
	bbox=(1.05,0.0,1.,1.),loc='upper left',fs=16,extra_legends=None):
	barspecs = bar_formats
	ion_names = ['K','NA','Na,Cal','MG','Cal']
	ptdins_names = ['PI2P','P35P','SAPI']
	rectangle_specs,patches,labels = {},[],[]
	for name,group in [('cation',ion_names),('ptdins_resname',ptdins_names)]:
		for item in group:
			sn = next(sn for sn in sns if work.meta[sn][name]==item)
			if sn not in sns_this: continue
			rectangle_specs[sn] = {}
			if name=='cation':
				rectangle_specs[sn]['fc'] = barspecs[sn]['edgecolor']
				rectangle_specs[sn]['lw'] = 0
				patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,**rectangle_specs[sn]))
				labels.append(work.meta[sn]['ion_label'])
			elif name=='ptdins_resname':
				patches.append(mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=1.0,fc='w',lw=3,ec='k',
					hatch=barspecs[sn]['hatch']*(1 if work.meta[sn]['ptdins_resname']=='PI2P' else 2)))
				labels.append(work.meta[sn]['ptdins_label'])
	if not patches: patches,labels = [mpl.patches.Rectangle((0,0),1.0,1.0,fc='w')],['empty']
	legend = ax.legend(patches,labels,loc=loc,fontsize=fs,
		ncol=ncol,title=title,bbox_to_anchor=bbox,labelspacing=1.2,
		handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)
	frame = legend.get_frame()
	frame.set_edgecolor('black')
	frame.set_facecolor('white')
	return legend,patches

def hbonds_bardat(format_name=None,ion_order=None,**kwargs):
	sns,post = kwargs['sns'],kwargs['post']
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
	#---group by the first key (the donor)
	for donor in list(set(zip(*keys)[0])):
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

#---block: plot settings
bar_style = ['original','candy'][-1]
live_plot_style = ['outline','reveal'][-1]
show_static,press = True,is_live
#---choices
comparison = ['asymmetric_all','symmetric_all'][0]
sns = work.specs['collections'][comparison]
ion_order = ['K','NA','Na,Cal','MG','Cal']
sns = reorder_by_ion(sns,ion_order)
#---aesthetics
mpl.rcParams['hatch.linewidth'] = 1.5
#---hatch color was controllable here in mpl v2.0.0 but they reverted
mpl.rcParams['hatch.color'] = 'w'
def dimmer(bar,child):
	"""Dim a bar."""
	bar.set_edgecolor(bar.get_facecolor())
	bar.set_facecolor('w')
	child.set_hatch(None)
sns_names = dict([(sn,work.meta[sn]['ptdins_resname']+' and '+
	work.meta[sn]['cation']) for sn in sns])
namer = lambda ((a,b),c): ('%s-%s'%(a,b))
bar_formats = make_bar_formats(sns,style=bar_style)
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
post = posts['both']
#---prepare the data
bardat = hbonds_bardat(sns=sns,post=post)
bars_premade = BarGrouper(figsize=(16,8),spacers={0:0,1:1.0,2:2.0},
	dimmer=dimmer,namer=namer,show_xticks=False,**bardat)

#---block: plot function and static plot
def hbonds_plotter(bars=None,**checks):
	"""Plot hydrogen bond counts."""
	#---prepare the data
	sns = [k for k,v in checks.items() if v]
	ion_order = ['K','NA','Na,Cal','MG','Cal']
	sns = reorder_by_ion(sns,ion_order)
	#---incoming bars creates the outlined-if-unchecked style plot
	#---! outline method is deprecated for hatching scheme
	if bars:
		for sn_checked in reorder_by_ion(checks.keys(),ion_order):
			if checks[sn_checked]: bars.on(*[(key,sn) for (key,sn) in bars.names if sn==sn_checked])
			else: bars.off(*[(key,sn) for (key,sn) in bars.names if sn==sn_checked])
	#---no incoming bars tells us to recompute the bar object for the bars that we want
	else:
		keys = reorder_by_ion([k for k,v in checks.items() if v],ion_order)
		empty = not keys
		bardat = hbonds_bardat(sns=checks.keys()[:1] if empty else keys)
		bars = BarGrouper(figsize=(16,8),spacers={0:0,1:1.0,2:2.0},
			dimmer=dimmer,namer=namer,show_xticks=False,empty=empty,**bardat)
	#---additional 
	bars.plot(wait=True)
	bars.ax.tick_params(axis='both', which='major',labelsize=16)
	bars.ax.set_ylabel(r'$\frac{\langle {N}_{bonds} \rangle}{{N}_{donor}{N}_{acceptor}}$',fontsize=26)
	bars.ax.yaxis.grid(which="minor",color='k',linewidth=0.5,alpha=0.5,zorder=0)
	bars.ax.tick_params(axis=u'both', which=u'both',length=0)
	bars.ax.set_title('hydrogen bonding',fontsize=40)
	legend,patches = legend_maker_stylized(bars.ax,sns_this=sns)
	return legend

#---initialize check marks for both live/static plotting
if is_live: checks = dict([(sn,widgets.ToggleButton(value=False,description=sns_names[sn])) for sn in sns])
else: checks = dict([(sn,True) for sn in sns])
#---show a static version of the plot
if show_static and is_live:
	hbonds_plotter(bars_premade,extras=bars_premade['extras'],**checks)
	plt.show()
#---batch plot several version of the data
if press:
	#---which figures to make
	press_specs = {
		'hydrogen_bonding.salt':{'posts_subkey':'salt'},
		'hydrogen_bonding.hbonds':{'posts_subkey':'hbonds'},
		'hydrogen_bonding.both_hbonds_salt':{'posts_subkey':'both'},}
	#---batch of figures
	for tag,spec in press_specs.items():
		bars_premade = BarGrouper(figsize=(16,8),spacers={0:0,1:1.0,2:2.0},
			dimmer=dimmer,namer=namer,show_xticks=False,
			**hbonds_bardat(sns=sns,post=posts[spec['posts_subkey']]))
		legend = hbonds_plotter(bars_premade,**checks)
		picturesave('fig.%s'%(tag),work.plotdir,
			backup=False,version=True,meta={},extras=[legend])

#---block: interactive plot
if is_live:
	def hbonds_plotter_interactive(**checks): 
		hbonds_plotter(bars={'outline':bars_premade,'reveal':None}[live_plot_style],**checks)
		plt.show()
	checks = dict([(sn,widgets.ToggleButton(value=sn==sns[0],description=sns_names[sn])) for sn in sns])
	widget = interactive(hbonds_plotter_interactive,**checks)
	ws = [widget.children[[w.description for w in widget.children[:-1]].index(sns_names[i])] for i in sns]
	display(HBox(ws),widget.children[-1])
	print('use the toggle buttons above to select different simulations for comparison')
	widget.update()
