#!/usr/bin/env python

"""
Plot hydration distributions.
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
required_variables = 'printers routine'.split()
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
	#---canonical globals list is loaded systematically 
	#---...but you have to load it into globals manually below
	global sns,scanrange,distributions,distances,data,calc,normalizers,middles
	#---custom reload sequence goes here
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
	"""
	Plot the RDFs for water near ions versus distance waters.
	Note that we copied this plot script to start over and simplify the calculation of the curves.
	"""
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
	colormaps = make_colormaps(fade_color='white')
	def make_rdf(sn,lims=None):
		"""Make the RDFs for lipids within a particular range."""
		if lims==None: return distributions[sn].mean(axis=0)
		else: 
			zone = np.all((data[sn]['data']['lipid_distances']>=lims[0],
				data[sn]['data']['lipid_distances']<=lims[1]),axis=0)
			#---! critical concatenation step. CHECK THAT THIS IS THE RIGHT AXIS ORDER!
			return distributions[sn][np.concatenate(zone)].mean(axis=0)
	def plot(ax,sns,zonelist=None,xlims=None,ylims=None,
		do_move_xlabel=False,do_peak_annotations=False,do_annotations=False,
		do_axlabels=True,is_inset=False,lw=2,do_legend=False,zone_alpha=None,zone_color=None,
		zone_lw=None,zone_zorder=None):
		"""Repeatedly call the plot for the main axis and the inset."""
		view_window = scanrange[:-1]<=farcut
		if not zonelist: zonelist = [None]
		#---plot the lipid distance-dependant g(r) for one simulation
		for snum,sn in enumerate(sns):
			#---index the zones for the colormap
			for znum,zone in enumerate(zonelist):
				curve = make_rdf(sn,lims=zone)
				curve = curve/normalize_window(curve)
				#---prepare settings
				kwargs = dict()
				if do_legend: kwargs.update(label='$\mathrm{%.1f-%.1f\AA}$%s'%(
					zone[0]*10,zone[1]*10,' (%s)'%work.meta[sn]['ion_label'] if znum==0 else ''))
				if zone_alpha: kwargs.update(alpha=zone_alpha(znum,len(zonelist),colormaps[sn]))
				if zone_color: kwargs.update(color=zone_color(znum,len(zonelist),colormaps[sn]))
				else: kwargs.update(color=colormaps[sn](0.0))
				if zone_lw: kwargs.update(lw=zone_lw(znum,len(zonelist)))
				else: kwargs.update(lw=lw)
				if zone_zorder: kwargs.update(zorder=zone_zorder(znum,len(zonelist)))
				ax.plot(10*middles[view_window],curve[view_window],'-',**kwargs)
			#---annotate the large peaks in the inset for extra clarity
			if do_peak_annotations:
				peak_arg = np.argmax(curve[view_window])
				peak_x,peak_y = 10*middles[view_window][peak_arg],curve[view_window][peak_arg]
				bbox_props = dict(boxstyle="round",fc="w",ec="k",alpha=1.0)
				ann = ax.text(peak_x+0.08*10,peak_y-1.0,work.meta[sn]['ion_label'],ha="center",va="center",size=9,
					rotation=0,bbox=bbox_props)
				extras.append(ann)
		ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower'))
		#---move the label so it doesn't overlap with tagboxes
		if do_move_xlabel: ax.xaxis.set_label_coords(0.3,0.08)
		#---annotations
		if do_annotations: 
			bbox_props = dict(boxstyle="round",fc="w",ec="k",alpha=1.0)
			#---draw horizontal and vertical lines for emphasis
			spot_l,spot_r,spot_b,spot_t,spot_m = [10*0.42,10*0.455,1.54,2.2,1.9]
			#for i in [spot_t,spot_b,spot_m]: ax.axhline(i,c='k',lw=0.5,alpha=0.5,zorder=2)
			for i in [spot_l,spot_r]: ax.axvline(i,c='k',lw=0.5,alpha=0.5,zorder=2)
			#---note that the bar fraction does not count the little offset between the arrow and the point
			for name,(i,j,k,l),rot,frac,label_move in zip(['distance','specificity','dehydration'],
				#---position the arrows
				[[spot_l,spot_t,spot_r,spot_t],[spot_l,spot_b,spot_l,spot_t],[spot_r,spot_m,spot_r,spot_b]],
				#---rotations and bar offsets
				[0,90,-90],[0.85,0.45,0.6],
				#---move the labels
				[[0,0.17],[-0.85,0],[0.6,0]]):
				extras.append(ax.annotate('',xy=(i,j),xycoords='data',xytext=(k,l),
					textcoords='data',arrowprops=dict(arrowstyle="<->",
						connectionstyle="bar,fraction=%.2f"%frac,ec="k",shrinkA=5,shrinkB=5)))
				extras.append(ax.text((i+k)/2.+label_move[0],(j+l)/2.+label_move[1],name,
					ha="center",va="center",size=10,rotation=rot,bbox=bbox_props))
			if False:
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
		ax.tick_params(axis='both',which='both',length=0)
		if do_axlabels:
			ax.set_xlabel('r ($\mathrm{\AA}$)',
				fontsize=fs.get('xlabel%s'%('_inset' if is_inset else''),14))
			ax.set_ylabel('g(r)',fontsize=fs.get('ylabel%s'%('_inset' if is_inset else''),14))
		if do_legend: 
			legend = ax.legend(title='cation-lipid distance $\mathrm{({d}_{c,l})}$',loc='lower right',
				frameon=False,fontsize=12)
			legend.get_title().set_fontsize('12')
	#---assemble the plot
	if incoming==None: fig,ax = plt.subplots()
	else: fig,ax = [incoming[i] for i in ['fig','ax']]
	xlims,ylims = (3.0,10.0),(0.,2.5)
	ax.set_aspect((xlims[1]-xlims[0])/(ylims[1]-ylims[0])*1.0)
	def zone_color(i,n,colormap): return colormap(float(i)/n)
	def zone_lw(i,n): return 3.5 if i==0 else 2
	def zone_zorder(i,n): return 4+n-i
	plot(ax,xlims=xlims,ylims=ylims,lw=3,
		zonelist=[(0.,0.22,),(0.22,0.46),(0.46,1)],
		sns=['membrane-v531','membrane-v532'],
		do_move_xlabel=True,do_legend=True,
		do_annotations=True,zone_color=zone_color,zone_lw=zone_lw,
		zone_zorder=zone_zorder)
	axins = inset_axes(ax,width="45%",height="45%",loc=1)
	plot(axins,xlims=(1,7),ylims=None,
		sns=['membrane-v538','membrane-v531','membrane-v532'],
		do_peak_annotations=True,is_inset=True)
	#---only save this figure if no incoming figure
	if incoming==None:
		picturesave('fig.hydration_distribution',work.plotdir,backup=False,version=True,
			meta={},extras=extras)

@register_printer
def water_distribution_sweep_plot_with_snapshots():
	"""
	Arrange snapshots alongside the RDF plot.
	"""
	#---hardcoded paths for now (switched to v15 from v12)
	snaps = {'folder':'fig.hydrogen_bonding.v15_ptdins_solvated','files':[
		'fig.snapshot.membrane-v531.fr703.795_796_o0.png',
		'fig.snapshot.membrane-v531.fr50.780_782_o6.png',
		#'fig.snapshot.membrane-v532.fr1.793_790_o1.png',
		'fig.snapshot.membrane-v532.fr11.779_781_o0.png',
		'fig.snapshot.membrane-v532.fr363.798_794_o6.png',
		'fig.snapshot.membrane-v532.fr3.773_775_o7.png',],}
	fig = plt.figure(figsize=(12,12))
	axes,extras = [],[]
	axes.append(plt.subplot2grid((3,3),(0,0),colspan=2,rowspan=2))
	for pos in [(0,2),(1,2),(2,2),(2,1),(2,0)]: axes.append(plt.subplot2grid((3,3),pos))
	#---bigger fonts
	ax = axes[0]
	ax.tick_params(axis='both',labelsize=18)
	ax.tick_params(axis='both',labelsize=18)
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
	#---add the RDF plot
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
	printers = list(set(printers if printers else []))
	if routine is None: routine = list(printers)	
	#---routine items are function names
	for key in routine: 
		status('running routine %s'%key,tag='printer')
		globals()[key]()

if __name__=='__main__': printer()
