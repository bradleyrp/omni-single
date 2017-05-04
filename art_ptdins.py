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
_art_words += ['uniquify','catalog','delve','delveset','delveset_list','BarGrouper','barspecs_stylized']

#---canonical colors for this project from brewer2mpl
import brewer2mpl

palette_colors = dict([(val,brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[key]) for key,val in enumerate('red blue green purple orange yellow brown pink grey'.split())])
	
def colorize(metadat,comparison='',resname=''):
	"""
	Master listing of colors for different PtdIns comparisons.
	"""
	if resname != '':
		colordict = {
			'DOPC':'grey',
			'DOPS':'red',
			'PIP2':'purple',}
		return palette_colors[colordict[resname]]
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
		return palette_colors[colordict[comparison][descriptor]]

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
	def __init__(self,dat,**kwargs): 
		self.dat = dat
		self.width = kwargs.get('width',1.0)
		self.spacers = kwargs.get('spacers',{0:0,1:0.5,2:1.0,3:2.0})
		self.dividers = kwargs.get('dividers',{})
		self.figsize = kwargs.get('figsize',None)
		self.bars,self.features = [],[]
		self.names,self.namesdict = [],{}
		self.cursor,self.level = 0.0,None
		self.extras = kwargs.get('extras',{})
		self.routes = list(catalog(dat))
		self.dimmer = kwargs.get('dimmer',None)
		self.namer = kwargs.get('namer',lambda x:x)
		self.show_xticks = kwargs.get('show_xticks',True)
		self.empty = kwargs.get('empty',False)
		#---mimic the incoming data structure and store bar properties
		self.barsave = copy.deepcopy(self.dat)
		#---plot specs for each bar
		self.specs = kwargs.get('specs',None)
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
		fig = plt.figure(figsize=self.figsize)
		ax = plt.subplot(111)
		#---! here instead of set_ ???
		xs,ws = [np.array([b[k] for bb,b in enumerate(self.bars)]) for k in 'lw']
		#---incoming values are tuples of mean,std
		yvals = np.array([b['v'] for bb,b in enumerate(self.bars)])
		ys,yerrs = yvals[:,0],yvals[:,1]
		alpha = 1.0 if not self.empty else 0.0
		#---! ytf did they change default alignment?
		rendered = ax.bar(xs,ys,width=ws,align='edge',zorder=3,alpha=alpha)
		# ax.errorbar(xs+self.width/2.,ys,yerr=yerrs,zorder=4,alpha=1.0,lw=5,c='w',ls='none')
		ax.errorbar(xs+self.width/2.,ys,yerr=yerrs,zorder=5,alpha=1.0,lw=3,c='k',ls='none')
		if self.show_xticks:
			xticks = np.mean(np.array([xs,ws+xs]),axis=0)
			ax.set_xticks(xticks)
			ax.set_xticklabels([self.namer(i) for i in self.names],rotation=90)
		else: ax.set_xticks([])
		if self.specs:
			#---post processing works on both the rendered objects and on the children
			for bb,(bar,child) in enumerate(zip(self.bars,rendered.get_children())):
				rendered[bb].set_color(bar.get('c','k'))
				hatch = self.bars[bb].get('hatch',None)
				child.set_hatch(hatch)
				#---apply the dimmer if this bar is "off"
				if not self.namesdict[self.names[bb]] and self.dimmer or self.empty: 
					#---dimmer takes the bar and the child
					#self.dimmer(bar=rendered[bb],child=child)
					child.set_hatch('||||')
					child.set_visible(False)
					#bar.set_visible(False)
			if False:
				#---! SEPARATE HATCH ROUTINE WHYYYYYY???
				children = rendered.get_children()
				for cc,child in enumerate(children):
					hatch = self.bars[cc].get('hatch',None)
					child.set_hatch(hatch)
					#---apply the dimmer if this bar is "off"
					if not self.namesdict[self.names[bb]] and self.dimmer: self.dimmer(rendered[bb])
		self.plot_extras(ax)
		if not wait: plt.show()
		self.ax = ax
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
	def plot_extras(self,ax):
		"""Handle extra annotations for different groups."""
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