#!/usr/bin/env python

"""
Plot hydration distributions.

Note that this is the third plot script. in the new style
neat: https://stackoverflow.com/questions/12200580/numpy-function-for-simultaneous-max-and-min
also neat: https://stackoverflow.com/questions/29321835/is-it-possible-to-get-color-gradients-under-curve-in-matplotlib

todo:
	1. finalize new plot style format
	2. better color gradient
	3. better way of specifying the zones, plus a good legend for them
	4. breakdown with inset
	5. no baselines please
	6. individual normalizers
"""

import scipy
import scipy.interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#---settings
binsize = 0.01
farcut = 1.05
#---we must declare variables here. this should match the top of the reload function
variables = 'sns,scanrange,distributions,distances,data,calc,normalizers,middles'.split(',')

from codes.looptools import basic_compute_loop

###---STANDARD

#---declare standard variables
required_variables = 'printers,routine'.split(',')
for v in variables+required_variables:
	if v not in globals(): globals()[v] = None

def register_printer(func):
	"""Add decorated functions to a list of "printers" which are the default routine."""
	global printers
	if printers is None: printers = []
	printers.append(func.__name__)
	return func

###---PLOTS

def histogram_stack(sn,index):
	"""Parallel function for computing histograms."""
	global distances,scanrange
	return np.histogram(distances[sn][index],bins=scanrange)[0]

def reload():
	"""Load everything for the plot only once."""
	#---canonical globals list
	#---!? can this be made programmatic?
	global sns,scanrange,distributions,distances,data,calc,normalizers,middles
	#---reload sequence goes here
	data,calc = plotload(plotname)
	sns = work.specs['collections']['position']+['membrane-v538']
	#---compute distance distributions
	cutoff = max([data[sn]['data']['water_distances'].max() for sn in sns])
	#---globals for parallel functions
	scanrange = np.arange(0,cutoff,binsize)
	#---distances are indexed by concatenated frames, then ions
	distances = dict([(sn,np.concatenate(data[sn]['data']['water_distances'])) for sn in sns])
	#---parallel compute
	looper = [dict(index=i,sn=sn) for sn in sns for i in range(len(distances[sn]))]
	incoming = np.array(basic_compute_loop(histogram_stack,looper=looper))
	distributions = dict([(sn,np.array([incoming[ii] for ii,i in enumerate(looper) 
		if i['sn']==sn])) for sn in sns])
	#---normalization factors
	middles = (scanrange[1:]+scanrange[:-1])/2
	areas = np.array([4*np.pi*binsize*middles[i]**2 for i in range(len(middles))])
	#---atoms in a nm3: (1000g/18g is mol/L) * 1L/1000ml * 1ml/cm3 / (10**7 nm/cm)**3 = 33.46
	water_density = 6.023*10**23*(1000.0/18)/1000/(10**(9-2))**3
	#---window to estimate bulk, lower than the dropoff, higher than the first two shells
	bulk_window_raw = (0.75,1.0)
	#---normalize the ion-water density for all RDF measurements
	#---note that this is repeated below for zone-specific normalizations
	normalizers = {}
	for sn in sns:
		bulk_window = np.where(np.all((scanrange>=bulk_window_raw[0],
			scanrange<=bulk_window_raw[1]),axis=0))[0]
		#---reestimate ion-water pseudo-density at supposed bulk-like distances
		water_density = (distributions[sn][slice(None,None),bulk_window]/areas[bulk_window]).mean()
		normalizers[sn] = areas*water_density

def water_distribution_sweep_plot(incoming=None,fs=None):
	"""Plot the RDFs for water near ions using a gradient that indicates proximity to lipids."""
	global colormaps
	if not fs: fs = {}
	extras = []
	def normalize_window(source,bulk_window_raw=(0.75,1.0),sliced=False):
		"""Compute far densities to normalize the RDF for particular zones."""
		areas = np.array([4*np.pi*binsize*middles[i]**2 for i in range(len(middles))])
		bulk_window = np.where(np.all((scanrange>=bulk_window_raw[0],
			scanrange<=bulk_window_raw[1]),axis=0))[0]
		#---reestimate ion-water pseudo-density at supposed bulk-like distances
		water_density = (source[bulk_window]/areas[bulk_window]).mean()
		return areas*water_density
	#---other functions can run this plot on an incoming axis
	#---choose the order for readability (must match keys above)
	sns_ordered = ['membrane-v538','membrane-v531','membrane-v532']
	#---defining the bins
	fenceposts = np.arange(0.16,1.0,step=0.4)
	fenceposts = np.linspace(0.16,1.0,num=10)
	zones_subsel_gradient = np.array([0,1,2,3,4,5,6,7,8])
	zones_subsubsel_gradient = [0]
	def make_colormaps(fade_color='white'):
		"""Colormaps for RDF curves, possibly with a gradient."""
		colormaps = {
			'membrane-v531':colorscale(bands=[
				colorize(work.meta['membrane-v531'],comparison='protonation'),
				fade_color],return_cmap=True)[1],
			'membrane-v532':colorscale(bands=[
				colorize(work.meta['membrane-v532'],comparison='protonation'),
				fade_color],return_cmap=True)[1],
			'membrane-v538':colorscale(bands=[
				colorize(work.meta['membrane-v538'],comparison='protonation'),
				fade_color],return_cmap=True)[1],}
		return colormaps
	#---set colormaps and sns for the primary comparison. settings placed near the plotter
	colormaps = make_colormaps()
	#---insets benefit from repeated calls to a single plot function
	def plot(ax,xlims=None,ylims=None,one=False,plot_line_style='many',
		do_opacity=False,axlabels=False,normed_windows=True,is_inset=False):
		global colormaps
		def minmax(a): return a.min(),a.max()
		lipid_distance_lims = minmax(np.array([minmax(data[sn]['data']['lipid_distances']) 
			for sn in sns]).reshape(-1))
		#---exclude some zones that have low populations for t he standard method
		if plot_line_style=='many': zones_subsel = np.array([0,1,6,7,8])
		elif plot_line_style=='single': zones_subsel = np.array([-1])
		elif plot_line_style in ['gradient','simple']: 
			zones_subsel,zones_subsubsel = zones_subsel_gradient,zones_subsubsel_gradient
		else: raise Exception('invalid plot_line_style %s'%plot_line_style)
		zones = dict([(sn,[np.all((data[sn]['data']['lipid_distances']>=i,
			data[sn]['data']['lipid_distances']<=j),axis=0) 
			for i,j in np.array(zip(fenceposts[:-1],fenceposts[1:]))[zones_subsel]]) for sn in sns])
		view_window = scanrange[:-1]<=farcut
		sns_ordered_core = ['membrane-v531','membrane-v532']
		#---plot the lipid distance-dependant g(r) for one simulation
		for snum,sn in enumerate(sns_ordered if is_inset else sns_ordered_core):
			if plot_line_style in ['many','single']:
				for znum,zone in enumerate(zones[sn]):
					curve = distributions[sn][np.concatenate(zone)].mean(axis=0)
					curve = curve/(normalizers[sn] if not normed_windows else normalize_window(curve))
					ax.plot(10*middles[view_window],curve[view_window],'-',lw=1,zorder=3,
						color=colormaps[sn](znum/float(len(zones[sn])-1)) if plot_line_style=='many' else 
							colormaps[sn](1))
			elif plot_line_style in ['gradient','fill_between','simple']:
				offset_color = 1.0
				#---! this method is currently set for adjacent zones only
				#---! ...the zones_subsubsel chooses the marker between adjacent zones
				#---! ...and the zones_subsel can be adjusted for coarse/fine zones
				#---! obviously this is clumsy but we are trying to be expedient
				#---! it might be worth replacing this wholesale since it's an important calculation
				for inum,index in enumerate(zones_subsubsel):
					zonepair = [zones[sn][zones_subsel[index+i]] for i in range(2)]
					curves = [distributions[sn][np.concatenate(z)].mean(axis=0) for z in zonepair]
					#---normalize between the two limits, effectively normalizing the sweep
					average_normalizer = np.array([normalize_window(c) for c in curves]).mean(axis=0)
					curves = [c/(normalizers[sn] if not normed_windows else average_normalizer) 
						for c in curves]
					#---fill between works pretty well but we have discarded it
					if plot_line_style=='fill_between':
						ax.fill_between(10*middles[view_window],y1=curves[1][view_window],
							y2=curves[0][view_window],
							lw=1,zorder=3,color=colormaps[sn]((inum)/
								float(len(zones_subsubsel)-1+offset_color)))
					elif plot_line_style=='gradient':
						ngrad = float(100)
						for s in range(int(ngrad))[::-1]:
							ys = ((curves[1][view_window]-curves[0][view_window])*float(s)/
								ngrad+curves[0][view_window])
							kwargs = dict(lw=1,zorder=3,
								color=colormaps[sn](float(s)/ngrad) 
									if not do_opacity else colormaps[sn](0.0),
								#---power level below tells you how fast the fade is
								alpha=((1.-s/ngrad)**2) if do_opacity else 1.0)
							ax.plot(10*middles[view_window],ys,**kwargs)
						ys = ((curves[1][view_window]-curves[0][view_window])*float(ngrad-1)/
							ngrad+curves[0][view_window])
						#---outline the mountain peaks
						kwargs.update(alpha=0.75,lw=0.5)
						ax.plot(10*middles[view_window],ys,**kwargs)
					elif plot_line_style=='simple':
						colormaps = make_colormaps(fade_color='black')
						lw = 3 if not is_inset else 2
						posts = [fenceposts[zones_subsel[index+i]] for i in range(2)]
						#---! need to check if labeling is accurate
						label = {'label':'%s $\mathrm{d_{c,l}=%.1f-%.1f\AA}$'%(work.meta[sn]['ion_label'],
							10*fenceposts[zones_subsel][index],10*fenceposts[zones_subsel][index+1])}
						ax.plot(10*middles[view_window],curves[0][view_window],
							lw=lw,zorder=3,color=colormaps[sn](0.0),**(label if not is_inset else {}))
						label = {'label':'%s $\mathrm{d_{c,l}=%.1f-%.1f\AA}$'%(work.meta[sn]['ion_label'],
							10*fenceposts[zones_subsel][index+1],10*fenceposts[zones_subsel][index+2])}
						if not is_inset: ax.plot(10*middles[view_window],curves[1][view_window],
							lw=2,zorder=3,color=colormaps[sn](0.3),alpha=0.5,**label)
					else: raise Exception('invalid plot_line_style %s'%plot_line_style)
			else: raise Exception('invalid plot_line_style %s'%plot_line_style)
			if is_inset and plot_line_style in ['many','single','simple']:
				if plot_line_style=='simple': curve = curves[0]
				peak_arg = np.argmax(curve[view_window])
				peak_x,peak_y = 10*middles[view_window][peak_arg],curve[view_window][peak_arg]
				bbox_props = dict(boxstyle="round",fc="w",ec="k",alpha=1.0)
				ann = ax.text(peak_x+0.08*10,peak_y-1.0,work.meta[sn]['ion_label'],ha="center",va="center",size=9,
					rotation=0,bbox=bbox_props)
				extras.append(ann)
		ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
		#---move the label so it doesn't overlap with tagboxes
		if not is_inset: ax.xaxis.set_label_coords(0.3,0.08)
		#---???
		if not is_inset: 
			ann = ax.annotate('',xy=(10*0.46,1.85),xycoords='data',xytext=(10*0.46,1.4),
				textcoords='data',arrowprops=dict(arrowstyle="<->",
					connectionstyle="bar,fraction=1.0",
				ec="k",shrinkA=5, shrinkB=5))
			extras.append(ann)
			bbox_props = dict(boxstyle="round",fc="w",ec="k",alpha=1.0)
			ann = ax.text(10*0.55,1.625,"specificity",ha="center",va="center",size=10,
				rotation=-90,bbox=bbox_props)
			extras.append(ann)
			ann = ax.annotate('',xy=(10*0.42,2.1),xycoords='data',xytext=(10*0.47,2.1),
				textcoords='data',arrowprops=dict(arrowstyle="<->",
					connectionstyle="bar,fraction=0.4",
				ec="k",shrinkA=5, shrinkB=5))
			extras.append(ann)
			ann = ax.text(10*0.445,2.35,"distance",ha="center",va="center",size=10,bbox=bbox_props)
			extras.append(ann)
		if xlims: ax.set_xlim(xlims)
		if ylims: ax.set_ylim(ylims)
		if one: ax.axhline(1.0,c='k',lw=0.5,zorder=1,alpha=0.35)
		ax.tick_params(axis='both',which='both',length=0)
		if axlabels:
			ax.set_xlabel('r ($\mathrm{\AA}$)',
				fontsize=fs.get('xlabel%s'%('_inset' if is_inset else''),14))
			ax.set_ylabel('g(r)',fontsize=fs.get('ylabel%s'%('_inset' if is_inset else''),14))
	#---assemble the plot
	if incoming==None: fig,ax = plt.subplots()
	else: fig,ax = [incoming[i] for i in ['fig','ax']]
	xlims,ylims = (3.0,10.0),(0.,2.8)
	ax.set_aspect((xlims[1]-xlims[0])/(ylims[1]-ylims[0])*1.0)
	#---key settings
	#---note that we replaced gradient with single which required a different colormap and a legend instead of colorbar
	plot_line_style = ['gradient','simple'][-1]
	plot(ax,xlims=xlims,ylims=ylims,plot_line_style=plot_line_style,one=False,do_opacity=False,axlabels=True)
	axins = inset_axes(ax,width="45%",height="45%",loc=1)
	plot(axins,xlims=(1,7),ylims=None,plot_line_style=plot_line_style,axlabels=True,is_inset=True)
	#---main colorbars
	norm = mpl.colors.Normalize(vmin=5, vmax=10)
	for snum,sn in enumerate(sns_ordered):
		#---inset colorbars to explain the gradient
		if plot_line_style=='gradient':
			sm = plt.cm.ScalarMappable(cmap=colormaps[sn],norm=plt.Normalize(vmin=0,vmax=1))
			sm._A = []
			axins = inset_axes(ax,width="3%",height="15%",loc=2,
				bbox_to_anchor=(0.08*(snum+1),0.,1.,0.90),bbox_transform=ax.transAxes,borderpad=0)
			if snum<len(sns_ordered)-1: cbar = plt.colorbar(sm,cax=axins,ticks=[])
			else: 
				ticks = [0.,1.]
				#---labels are set for one zone only
				if len(zones_subsubsel_gradient)!=1: raise Exception('labels need only one sub-sub-selection')
				lims = [fenceposts[zones_subsel_gradient[zones_subsubsel_gradient[0]+i]] for i in range(2)]
				cbar = plt.colorbar(sm,cax=axins,ticks=ticks)
				axins.set_yticklabels(['%.2f'%i for i in lims],fontsize=fs.get('colorbar_label',8))
				axins.set_xlabel('cation-lipid\ndistance ($\mathrm{\AA}$)',
					fontsize=fs.get('colorbar_label',8),labelpad=10)
			axins.set_title(work.meta[sn]['ion_label'],fontsize=fs.get('colorbar_label',8))
			cbar.ax.tick_params(axis='both',which='both',length=0)
		elif plot_line_style=='simple':
			ax.legend(loc='lower right',frameon=False)
		else: raise Exception('invalid plot_line_style %s'%plot_line_style)
	#---only save this figure if no incoming figure
	if incoming==None:
		picturesave('fig.hydration_distribution',work.plotdir,backup=False,version=True,
			meta={},extras=extras)

@register_printer
def water_distribution_sweep_plot_with_snapshots():
	fig = plt.figure(figsize=(12,12))
	axes = []
	axes.append(plt.subplot2grid((3,3),(0,0),colspan=2,rowspan=2))
	for pos in [(0,2),(1,2),(2,2),(2,1),(2,0)]: axes.append(plt.subplot2grid((3,3),pos))
	#---bigger fonts
	ax = axes[0]
	ax.tick_params(axis='both',labelsize=18)
	ax.tick_params(axis='both',labelsize=18)

	#---!!! needs centralized
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

	#---hardcoded paths for now (switched to v15 from v12)
	snaps = {'folder':'fig.hydrogen_bonding.v15_ptdins_solvated','files':[
		'fig.snapshot.membrane-v531.fr703.795_796_o0.png',
		'fig.snapshot.membrane-v531.fr50.780_782_o6.png',
		#'fig.snapshot.membrane-v532.fr1.793_790_o1.png',
		'fig.snapshot.membrane-v532.fr11.779_781_o0.png',
		'fig.snapshot.membrane-v532.fr363.798_794_o6.png',
		'fig.snapshot.membrane-v532.fr3.773_775_o7.png',],}

	extras = []
	letter_reord = [0,1,4,3,2]
	#---ensure the order matches below
	sns_order = ['membrane-v531' for i in range(2)]+['membrane-v532' for i in range(3)]
	tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
	for inum,snap_fn in enumerate(snaps['files']):
		image = sidestack_images([os.path.join(work.plotdir,snaps['folder'],snap_fn)])
		#---make the image square
		this_shape = image.shape[:2]
		extra_buffer = (max(this_shape)-min(this_shape))/2
		white_buffer = np.ones((extra_buffer,max(this_shape),3))
		if np.argmin(this_shape)==1: reorder = (0,1,2),(1,0,2),(1,0,2)
		elif np.argmin(this_shape)==0: reorder = (0,1,2),(0,1,2),(0,1,2)
		else: raise Exception('incorrect dimensions of image %s'%image.shape)
		image = np.concatenate((white_buffer.transpose(reorder[0]),
			image.transpose(reorder[1]),white_buffer.transpose(reorder[0]))).transpose(reorder[2])
		ax = axes[1+inum]
		ax.imshow(image)
		ax.axis('off')
		tb = ax.text(0.0,1.0,chr(ord('A')+letter_reord[inum]),fontsize=14,
			bbox=tagbox,rotation=0,ha="center",va="top",color='k',
			transform=ax.transAxes)
		extras.append(tb)
		tb = ax.text(0.2,1.0,work.meta[sns_order[inum]]['ion_label'],fontsize=14,
			bbox=tagbox,rotation=0,ha="center",va="top",color='k',
			transform=ax.transAxes)
		extras.append(tb)
	water_distribution_sweep_plot(incoming=dict(fig=fig,ax=axes[0]),
		fs=dict(xlabel=20,ylabel=20,ylabel_inset=16,xlabel_inset=16,colorbar_label=12))
	picturesave('fig.hydration_distribution.snapshots',
		work.plotdir,backup=False,version=True,meta={},extras=extras)

###---STANDARD

def printer():
	"""Load once per plot session."""
	global variables,routine,printers
	#---reload if not all of the globals in the variables
	if any([v not in globals() or globals()[v] is None for v in variables]): reload()
	#---after loading we run the printers
	printers = list(set(printers if printers!=None else []))
	if routine is None: routine = list(printers)	
	#---routine items are function names
	for key in routine: 
		status('running routine %s'%key,tag='printer')
		globals()[key]()

if __name__=='__main__': printer()
