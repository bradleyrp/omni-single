#!/usr/bin/env python

"""
ART: PtdIns project
Note: you can run the ipynb header to refresh this.
"""

#---forgo the use of the colors/labels since this project has complicated rendering
colors = None
labels = None

#---CUSTOM additions to globals
_art_words = ['colorize','barmaker']
_art_words += ['uniquify','catalog','delve','delveset','delveset_list',
	'BarGrouper','barspecs_stylized','colorscale',
	'make_bar_formats','legend_maker_stylized','sidestack_images','get_blank_border',
	'ptdins_manuscript_settings','residue_codes']
#! extras for actinlink
_art_words += ['actinlink_monolayer_indexer','actinlink_replicate_mapping','actinlink_color_by_simulation',
	'actinlink_extra_labels','color_by_simulation','actinlink_sns_mdia2','actinlink_color_by_simulation']

#---canonical colors for this project from brewer2mpl
import brewer2mpl
import re

palette_colors = dict([(val,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[key]) for key,val in enumerate('red blue green purple orange yellow brown pink grey'.split())])
	
def colorize(metadat,comparison='',resname='',named=False,overrides=None):
	"""
	Master listing of colors for different PtdIns comparisons.
	"""
	if resname != '':
		#---switched CHL1 from green to gray for the hydrogen bonding plots
		colordict = {
			'DOPC':'grey',
			'POPC':'grey',
			'DOPS':'red',
			'DOPE':'blue',
			'PIP2':'purple',
			'PI2P':'purple',
			'P35P':'purple',
			'PIPU':'purple',
			'PIPP':'purple',
			'SAPI':'purple',
			'PtdIns':'purple',
			'CHL1':'grey'}
		if overrides: colordict.update(**overrides)
		if named: return colordict[resname]
		else: return palette_colors[colordict[resname]]
	else:
		colordict = {
			'ions':{
				('NA',):'green',
				('MG',):'red',
				('Cal',):'blue',
				('K',):'grey',
				('Na,Cal',):'blue',},
			'phosphate_position':{
				('NA','PI2P'):'green',
				('MG','PI2P'):'red',
				('Cal','PI2P'):'blue',
				('MG','P35P'):'purple',
				('Cal','P35P'):'orange',},
			'protonation':{
				('NA',):'green',
				('MG',):'red',
				('Cal',):'blue',
				('K',):'grey',
				('Na,Cal',):'blue',
				('NA','PI2P'):'green',
				('Na,Cal','PI2P',):'blue',
				('K','PI2P'):'grey',
				('MG','PI2P'):'red',
				('Cal','PI2P'):'blue',
				('MG','P35P'):'purple',
				('Cal','P35P'):'orange',
				('NA','PI2P'):'green',
				('NA','PIPU'):'brown',
				#---hackish
				('NA','SAPI'):'brown',
				('NA','PIPP'):'grey',},}
		compare_to_value = {
			'symmetric':['cation'],	
			'asymmetric':['cation'],
			'symmetric_all':['cation','ptdins_resname'],	
			'asymmetric_all':['cation','ptdins_resname'],
			'phosphate_position':['cation','ptdins_resname'],
			'protonation':['cation','ptdins_resname']}		
		#---aliases
		colordict['position'] = colordict['phosphate_position']
		compare_to_value['position'] = compare_to_value['phosphate_position']
		colordict['charge'] = colordict['protonation']
		compare_to_value['charge'] = compare_to_value['protonation']
		colordict['symmetric'] = dict(colordict['ions'])
		colordict['asymmetric'] = dict(colordict['ions'])
		colordict['symmetric_all'] = dict(colordict['protonation'])
		colordict['asymmetric_all'] = dict(colordict['protonation'])
		descriptor = tuple([metadat[i] for i in compare_to_value[comparison]])
		if named: return colordict[resname]
		else: return palette_colors[colordict[comparison][descriptor]]

def barmaker(ax,yvals,**kwargs):

	"""
	standardized bar plotter with groupings
	bars are plotted as a list of lists
	"""

	#---defaults have a space between bar groups which are otherwise flush
	width = kwargs.pop('width',1.0)
	gap_out = kwargs.pop('gap',0.1)
	gap_in = kwargs.pop('gap_in',0.0)
	zero_threshold = kwargs.pop('zero_threshold',0.0)
	lw_zero = kwargs.pop('lw_zero',3.0)
	barspecs_all = kwargs.get('barspecs_all',{})
	barspecs = kwargs.get('barspecs',{})
	barspecs_stack = kwargs.get('barspecs_stack',{})
	xoffset = kwargs.get('xoffset',0)

	#---if there are no bar groups and/or stacks we reformat the objects for the triple loop
	if all([type(i)!=list for i in yvals]): yvals = [[[i]] for i in yvals]
	elif all([all([type(j)!=list for j in i]) for i in yvals]):
		yvals = [[[j] for j in i] for i in yvals]
	assert all([all([type(j)==list for j in i]) for i in yvals]),'yvals have weird dimensions'

	ymax,xpos,xvals = 0,xoffset,[]
	#---loop over groups of bars
	for ii,i in enumerate(yvals):
		xpos += gap_out
		xvals_grp = []
		#---loop over bars within the group
		for jj,j in enumerate(i):
			xvals_grp.append(xpos)
			bottom = 0
			#---loop over stacks of bars
			for kk,k in enumerate(j):
				#---assemble the right specs for the bar (by group) and stack position
				this_barspecs = dict(barspecs_all)
				if barspecs: this_barspecs.update(**barspecs[ii][jj])
				if barspecs_stack: this_barspecs.update(**barspecs_stack[ii][jj][kk])
				#---zero thresholds get a nice line to remind you they are present
				if k <= zero_threshold and bottom == 0:
					ax.plot([xpos,xpos+width],[k,k],color=this_barspecs.get('color','k'),
						lw=2,zorder=2,clip_on=False)
					ax.plot([xpos,xpos+width],[k,k],color='k',lw=4,zorder=1,clip_on=False)
				ax.bar(xpos,k,width=width,bottom=bottom,**this_barspecs)
				ymax = ymax if ymax>k+bottom else k+bottom
				bottom += k
			xpos += width
			xpos += gap_in
		xvals.append(xvals_grp)

	return {'ymax':ymax,'xvals':xvals,'width':width}

import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
import copy

class BarGrouper:
	"""
	Process arbitrarily nested lists of 
	A bar plot is just a linear object with bars and spaces and possibly dividers.
	Each bar object has a left and right object.
	We wish to retain sequences from the incoming data.
	"""
	def __init__(self,**kwargs): 
		self.dat = kwargs.get('dat',kwargs.get('bardat',None))
		self.width = kwargs.get('width',1.0)
		self.spacers = kwargs.get('spacers',{0:0,1:0.5,2:1.0,3:2.0})
		self.dividers = kwargs.get('dividers',{})
		self.figsize = kwargs.get('figsize',None)
		self.bars,self.features = [],[]
		self.names,self.namesdict = [],{}
		self.cursor,self.level = 0.0,None
		self.extras = kwargs.get('extras',{})
		self.routes = list(catalog(self.dat))
		self.error_bar_lw = kwargs.get('error_bar_lw',3)
		self.dimmer = kwargs.get('dimmer',None)
		self.namer = kwargs.get('namer',lambda x:x)
		self.show_xticks = kwargs.get('show_xticks',True)
		self.empty = kwargs.get('empty',False)
		self.ax = kwargs.get('ax',None)
		if self.ax: self.fig = kwargs['fig']
		#---mimic the incoming data structure and store bar properties
		self.barsave = copy.deepcopy(self.dat)
		#---plot specs for each bar
		self.specs = kwargs.get('specs',kwargs.get('barspec',None))
		self.proc()
	def proc(self):
		"""Process the bars."""
		#---unpack the dat recursively and make one bar for each child
		for path,(name,val) in self.routes:
			if not self.level: self.level,self.last_path = len(path),None
			#---change in depth interpreted before we add the bar
			self.next(self.level,len(path),path)
			#---make the bar
			bar = {'l':self.cursor,'v':val,'w':self.width,'r':self.cursor+self.width}
			#---apply special formatting via specs (which should mimic the structure of the data)
			if self.specs:
				#---the first item in the tuple is the name 
				try: spec = delve(self.specs,*path)
				except: spec = {}
				#---valid features includes: colors, hatches
				bar['c'] = spec.get('c','k')
				if 'hatch' in spec: bar['hatch'] = spec['hatch']
			#---bars go to a list and a structure with the same topology as the incoming
			self.bars.append(bar)
			self.names.append(name)
			self.namesdict[name] = True
			delveset_list(self.barsave,bar,*path)
			self.bar_done()
			self.last_path = path
	def next(self,before,after,path):
		"""When the group level changes we apply spacers and dividers."""
		change = abs(after-before)
		#---update the spacer for the features
		space = self.spacers[change]
		#---if there is not a level change there might be a level break
		if not space:
			#---check for a level change one above the current level
			if self.last_path and path[:-1]!=self.last_path[:-1]:
				space = self.spacers[len(path[:-1])]
		self.cursor += space
		self.level = after
	def bar_done(self):
		"""When a bar is completed we advance the cursor."""
		self.cursor += self.width
	def plot(self,wait=False):
		if not self.ax:
			fig = self.fig = plt.figure(figsize=self.figsize)
			ax = self.ax = plt.subplot(111)
		else: ax,fig = self.ax,self.fig
		#---! here instead of set_ ???
		xs,ws = [np.array([b[k] for bb,b in enumerate(self.bars)]) for k in 'lw']
		#---incoming values are tuples of mean,std
		yvals = np.array([b['v'] for bb,b in enumerate(self.bars)])
		ys,yerrs = yvals[:,0],yvals[:,1]
		alpha = 1.0 if not self.empty else 0.0
		#---! ytf did they change default alignment?
		#---no outline around the bars. the hatch colors are set below, and match the edge color
		rendered = ax.bar(xs,ys,width=ws,align='edge',zorder=3,alpha=alpha,lw=0)
		ax.errorbar(xs+self.width/2.,ys,yerr=yerrs,zorder=5,alpha=1.0,lw=self.error_bar_lw,c='k',ls='none')
		if self.show_xticks:
			xticks = np.mean(np.array([xs,ws+xs]),axis=0)
			ax.set_xticks(xticks)
			ax.set_xticklabels([self.namer(i) for i in self.names],rotation=90)
		else: ax.set_xticks([])
		if self.specs:
			#---post processing works on both the rendered objects and on the children
			for bb,(bar,child) in enumerate(zip(self.bars,rendered.get_children())):
				#---! between mpl v2.0.0 and v2.0.2 they must have changed behavior so hatches are 
				#---! ...using the edge color. note that you can set the color for each bar specifically
				#---! ...in the ax.bar command above but we retain full control down here
				rendered[bb].set_hatch(self.bars[bb].get('hatch',None))
				rendered[bb].set_color(bar.get('c','k'))
				#---by mpl v2.0.2 (or earlier?) the hatches have the edge color
				#---hatch color is the edge color set here
				rendered[bb].set_edgecolor('k')
				#---apply the dimmer if this bar is "off"
				if not self.namesdict[self.names[bb]] and self.dimmer: self.dimmer(rendered[bb])
		self.plot_extras()
		if not wait: plt.show()
		#---! mimics the divider code
		if self.bars:
			lims = zip(*[(i['l'],i['r']) for i in self.bars])
			left,right = min(lims[0]),max(lims[1])
			ax.set_xlim(left-self.width,right+self.width)
		fig.canvas.draw()
	def name_to_path(self,name,repeats=False):
		"""Get the paths associated with a name."""
		#---linear search through routes
		matches = [ii for ii,(path,(key,val)) in enumerate(self.routes) if key==name]
		if not repeats and len(matches)>1: 
			#---! add feature that allows redundant labels and somehow still groups correctly
			print(name)
			print(matches)
			raise Exception('cannot provide a unique path for name %s. use unique names'%str(name))
		elif len(matches)==0: raise Exception('no bar named "%s"'%name)
		else: return self.routes[matches[0]]
	def bars_by_names(self,names):
		"""Look up bars by names"""
		#---get the paths for each name
		paths = [self.name_to_path(name)[0] for name in names]
		#---find the minimum common path
		for depth in range(max([len(i) for i in paths]))[::-1]:
			if len(list(set([tuple(p[:depth]) for p in paths])))==1: break
		#---look up the bars in the group
		bars = [delve(self.barsave,*p) for p in [tuple(p) for p in paths]]
		return bars
	def plot_extras(self):
		"""Handle extra annotations for different groups."""
		ax = self.ax
		for name,extra in self.extras:
			if name=='divider':
				bars = self.bars_by_names(extra['names'])
				lims = zip(*[(i['l'],i['r']) for i in bars])
				left,right = min(lims[0]),max(lims[1])
				for spot in [left-self.spacers[1],right+self.spacers[1]]: 
					ax.axvline(spot,c=extra.get('c','k'))
				label = extra.get('label',None)
				if label:
					xpos = (left+right)/2.
					tform = mpl.transforms.blended_transform_factory(ax.transData,ax.transAxes)
					ax.annotate(label,xy=(xpos,0),xycoords=tform, 
						xytext=(xpos,0.9),textcoords=tform,fontsize=extra.get('fs',None),
						ha='center',va='center',arrowprops={},
						bbox=extra.get('bbox',None))
			elif name=='xlabel':
				bars = self.bars_by_names(extra['names'])
				lims = zip(*[(i['l'],i['r']) for i in bars])
				left,right = min(lims[0]),max(lims[1])
				ax.set_xticks(list(ax.get_xticks()) + [(left+right)/2.])
				ax.set_xticklabels([x.get_text() for x in ax.get_xticklabels()][:-1]+[extra['label']],
					rotation=extra.get('rotation',0),fontsize=extra.get('fs',None))
			else: raise Exception('no extra named %s'%name)
	def off(self,*args): self._change_state(False,*args)
	def on(self,*args): self._change_state(True,*args)
	def _change_state(self,val,*args):
		for arg in args: self.namesdict[arg] = val

def legend_maker_stylized(ax,sns_this,work,title=None,ncol=1,
	bbox=(1.05,0.0,1.,1.),loc='upper left',fs=16,extra_legends=None,fancy=True,
	comparison_spec=None,lipid_resnames=None,sns_explicit=None,bar_formats=None):
	"""
	############## copied from hydrogen bonding!
	needs work compared to the original
	"""
	if not comparison_spec: global comparison
	if not bar_formats: bar_formats = globals()['bar_formats']
	#---custom items in the legend based on the comparison
	if comparison_spec:
		ion_names = comparison_spec.get('ion_names',[])
		ptdins_names = comparison_spec.get('ptdins_names',[])
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
	###---!!!!!!!!!!!!!!!!!!!!
	if 'bar_formats_style' in globals() and bar_formats_style=='actinlink':
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
		handleheight=2.0,markerscale=0.5,shadow=fancy,fancybox=fancy)
	frame = legend.get_frame()
	frame.set_edgecolor('black' if fancy else 'white')
	frame.set_facecolor('white')
	return legend,patches

def make_bar_formats(sns,work,style='candy'):
	"""
	Make bar formats, most recently, candy-colored bars with hatches.
	"""
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

def figspec_to_panelplot(layout):
    """
    Convert a list of explicit "plot these sns here" into the layout expected by panelplot.
    """
    #---! no nesting
    out,ins = {'grid':[1,1]},[]
    #---infer the dimensions of the inner grid
    rcs = [(i['col'],i['row']) for i in layout]
    rows,cols = [sorted(i) for i in zip(*rcs)]
    rmin,rmax,cmin,cmax = min(rows),max(rows),min(cols),max(cols)
    ins = [{'grid':[cmax-cmin+1,rmax-rmin+1]}]
    return {'out':out,'ins':ins}

def uniquify(array):
    """Get unique rows in an array."""
    #---contiguous array trick
    alt = np.ascontiguousarray(array).view(
        np.dtype((np.void,array.dtype.itemsize*array.shape[1])))
    unique,idx,counts = np.unique(alt,return_index=True,return_counts=True)
    #---sort by count, descending
    idx_sorted = np.argsort(counts)[::-1]
    return idx[idx_sorted],counts[idx_sorted]

def catalog(base,path=None):
    """Traverse all paths in a nested dictionary."""
    if not path: path=[]
    if isinstance(base,list):
        for i,x in enumerate(base):
            local_path = path[:]+[i]
            for b in catalog(base[i],local_path): yield b
    else: yield path,base
def delve(o,*k): 
    """Return items from a nested dict."""
    return delve(o[k[0]],*k[1:]) if len(k)>1 else o[k[0]]
def delveset_list(o,value,*k): 
    """Utility function for adding a path to a nested list. Copied from delveset for dict."""
    if len(k)==0: raise Exception('deepset needs a path')
    elif len(k)==1: o[k[0]] = value
    else:
        if k[0] not in range(len(o)): raise Exception('invalid delveset_list path: o="%s", k="%s"'%(o,k))
        delveset_list(o[k[0]],value,*k[1:])
        
"""
Common figures specifications.
"""

_art_words += ['figplacer','figlayout','figplace','blank_unused_axes']

def figplacer(sn,figplace):
	"""
	Return the outer and inner axes indices for a particular simulation given a figplace.
	"""
	rown = next(row for row,i in enumerate(figplace) if sn in i)
	coln = next(i.index(sn) for row,i in enumerate(figplace) if sn in i)
	return rown,coln

figlayout = {}
figlayout['4x4'] = {'out':{'grid':[1,1]},'ins':[{'grid':[2,2]}]}
figlayout['summary1a'] = {'out':{'grid':[4,1]},'ins':[
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},]}

figlayout['summary1'] = {'out':{'grid':[5,1]},'ins':[
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},
	{'grid':[1,3],'wspace':0.5},]}
	
#---FIGURE PLACEMENTS	

figplace = {}
figplace['summary1a'] = [
	['membrane-v509','membrane-v510','membrane-v511'],
	['membrane-v530','membrane-v531','membrane-v532'],
	['membrane-v536','membrane-v533','membrane-v534'],
	['membrane-v538','membrane-v514','membrane-v515'],]	

figplace['summary1'] = [
	['membrane-v509','membrane-v510','membrane-v511'],
	['membrane-v530','membrane-v531','membrane-v532'],
	['membrane-v536','membrane-v533','membrane-v534'],
	['membrane-v538','membrane-v514','membrane-v515'],
	['membrane-v543','membrane-v542',None],]	

def blank_unused_axes(axes,fig,figplace):
	"""Remove blank figures."""
	for rr,cc in [(ii,jj) for ii,i in enumerate(figplace) 
		for jj,j in enumerate(i) if j==None]:
		fig.delaxes(axes[rr][cc])

def colorscale(name=None,count=256,cmap=None,bands=None,
	reverse=False,sharp=False,zero=False,return_cmap=False,nsegs=None):
	"""
	Divide a matplotlib color map into discrete colors.
	Copied from diffusion_ions_zoned.
	"""
	from matplotlib.colors import LinearSegmentedColormap
	if cmap != None: thiscmap = cmap
	elif name != None: thiscmap = plt.cm.get_cmap(name)
	elif bands != None:
		vmax = len(bands)-1.
		color_list = [(ii/vmax,i) for ii,i in enumerate(bands)]
		segs = (int(vmax)+1 if sharp else (256 if not zero else 255)) if nsegs==None else nsegs
		thiscmap = LinearSegmentedColormap.from_list('',color_list,N=segs)
	else: raise Exception('unclear inputs')
	grads = [thiscmap(i) for i in np.array(range(0,count)[::(-1 if reverse else 1)])/float(count)]
	if return_cmap: return grads,thiscmap
	else: return grads

import matplotlib.image as mpimg

def get_blank_border(img):
	"""Return the border limits for removing whitespace."""
	nonblank = np.any(img!=1.0,axis=2)*1
	lims = [map(lambda x:(x[0]-1,x[-1]+1),
		[np.where(np.any(nonblank!=0,axis=j))[0]])[0] for j in range(2)]
	return lims

def sidestack_images(fns):
	"""Combine snapshots into one image with same zoom and minimal whitespace."""
	#---preassemble images for a particular row
	imgs = [mpimg.imread(fn) for fn in fns]
	borders = np.array([get_blank_border(img) for img in imgs])
	lims = np.array([[i.min() for i in borders.T[0]]]+[[i.max() for i in borders.T[1]]]).T
	#---stack horizontally
	buffer_width = 100
	imgs_zoomed = [img[slice(*lims[1]),slice(*lims[0])] for img in imgs]
	imgs_cropped = [i[:,slice(*(get_blank_border(i)[0]))] for i in imgs_zoomed]
	buffer_strip = np.ones((imgs_cropped[0].shape[0],buffer_width,3))
	#---add a vertical buffer strip
	imgs_combo = np.concatenate([np.concatenate((i,buffer_strip),axis=1) 
		for i in imgs_cropped[:-1]]+[imgs_cropped[-1]],axis=1)
	return imgs_combo

def ptdins_manuscript_settings():
	"""Alternate colors for late-stage ptdins manuscript preparation."""
	colors = dict([(key,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[val])
		for key,val in {
		'red':0,'blue':1,'green':2,'purple':3,'orange':4,
		'yellow':5,'brown':6,'pink':7,'grey':8,}.items()])
	colors['white'] = mpl.colors.ColorConverter().to_rgb("#000000")
	colors['pink'] = mpl.colors.ColorConverter().to_rgb("#f1948a")
	colors['beige'] = mpl.colors.ColorConverter().to_rgb("#C3C3AA")
	bgw = 0.3
	colors['bluegreen'] = np.average((colors['blue'],colors['green']),weights=[1-bgw,bgw],axis=0)
	bgw2 = 0.6
	colors['bluegreen2'] = np.average((colors['blue'],colors['green']),weights=[1-bgw2,bgw2],axis=0)
	scale = 0.5
	colors['light_blue'] = np.average((colors['blue'],colors['white']),weights=[1-scale,scale],axis=0)
	colors_ions = {'NA':'green','Na,Cal':'bluegreen','MG':'pink','Cal':'blue','K':'grey',}
	hatches_lipids = {'PI2P':'//','P35P':'-','PIPU':'xx','PIPP':'++','SAPI':''}
	return dict(colors_ions=colors_ions,colors=colors,hatches_lipids=hatches_lipids)

residue_codes = {'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E',
	'SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C','SEL':'U','GLY':'G','PRO':'P',
	'ALA':'A','ILE':'I','LEU':'L','MET':'M','PHE':'F','TRP':'W','TYR':'Y','VAL':'V'}

#! actinlink settings

def actinlink_monolayer_indexer(sn,abstractor):
	if abstractor=='lipid_chol_com':
		#! had to add 'mdia2bilayer10' for this to get lipid_rdfs to work
		if sn in ['mdia2bilayer10','mdia2bilayer10_2','mdia2bilayer30_2']: return 1
		else: return 0
	elif abstractor=='lipid_com':
		if sn in ['mdia2bilayerphys2','mdia2bilayer30_2']: return 1
		else: return 0
	else: raise Exception

actinlink_replicate_mapping = [('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
	('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
	('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
	('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])]

actinlink_extra_labels = {
	'pip2_20_no_chol':r'mDia2, 20% $PIP_2$, no CHOL ($\times2$)',
	'pip2_20':r'mDia2, 20% $PIP_2$ ($\times2$)',
	'pip2_30':r'mDia2, 30% $PIP_2$ ($\times2$)',
	'pip2_10':r'mDia2, 10% $PIP_2$ ($\times2$)',}

def actinlink_color_by_simulation(sn):
	"""
	Choose colors for each simulation across many plots.
	"""
	replicate_mapping = actinlink_replicate_mapping
	#! remove from actinlink_bonds_analysis
	#! move to art? this is used elsewhere, in plot-ptdins_partners.py
	# tetradic colors via https://www.sessions.edu/color-calculator/
	colors = ['#ff6c28','#28c6ff','#a928ff','#ffd828']
	# tetradic colors via http://paletton.com/#uid=7030Z0kqbujggFvlnx6qcY+wDkl
	#! these colors must exactly match the ordering in the replicate matching
	colors = ['#f2b52c','#2e47a4','#22b93c','#f23a2c']
	refs = ['^mdia2bilayer_nochl','^mdia2bilayer10','^mdia2bilayer(?!(_nochl|[0-9]))','^mdia2bilayer30']
	ref_to_color = dict(zip(refs,colors))
	matches = [color for ref,color in ref_to_color.items() if re.match(ref,sn)]
	if len(matches)>1: raise Exception
	elif len(matches)==0:
		#! color order must match the replicate mapping order
		return dict(zip(zip(*replicate_mapping)[0],colors))[sn]
	else: return matches[0]

actinlink_sns_mdia2 = ['mdia2bilayer_nochl2','mdia2bilayer_nochl3','mdia2bilayer10','mdia2bilayer10_2',
	'mdia2bilayerphys','mdia2bilayerphys2','mdia2bilayer30','mdia2bilayer30_2']

def make_legend(ax,keys=None):
	#! copied from plot-lipid_rdfs.pf
	# legend at last tile
	legendspec = []
	for sn_general,sns in replicate_mapping:
		if not keys or sn_general in keys:
			legendspec.append(dict(name=extra_labels[sn_general],
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.))
	frame = legend.get_frame()
	frame.set_edgecolor('k')
	frame.set_facecolor('w')
	return legend
