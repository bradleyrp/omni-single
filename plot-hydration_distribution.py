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
	#---separate normalizer for each simulation because different ions
	#---note that these normalizers are for all zones, so it might not approach unity on subsets
	normalizers = {}
	for sn in sns:
		bulk_window = np.where(np.all((scanrange>=bulk_window_raw[0],
			scanrange<=bulk_window_raw[1]),axis=0))[0]
		#---reestimate ion-water pseudo-density at supposed bulk-like distances
		water_density = (distributions[sn][:,bulk_window]/areas[bulk_window]).mean()
		normalizers[sn] = areas*water_density

@register_printer
def water_distribution_sweep_plot():
	"""Plot the RDFs for water near ions using a gradient that indicates proximity to lipids."""
	#---we can use the normalizers for each simulation computed above or compute them for a specific window
	normed_windows = True
	def normalize_window(source,bulk_window_raw=(0.75,1.0)):
		"""..."""
		areas = np.array([4*np.pi*binsize*middles[i]**2 for i in range(len(middles))])
		bulk_window = np.where(np.all((scanrange>=bulk_window_raw[0],
			scanrange<=bulk_window_raw[1]),axis=0))[0]
		#---reestimate ion-water pseudo-density at supposed bulk-like distances
		water_density = (source[bulk_window]/areas[bulk_window]).mean()
		return areas*water_density
	#---set colormaps and sns for the primary comparison. settings placed near the plotter
	fade_color = 'gray'
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
	sns_ordered = ['membrane-v538','membrane-v531','membrane-v532']
	#---insets benefit from repeated calls to a single plot function
	def plot(ax,xlims=None,ylims=None,one=False,plot_line_style='many',do_opacity=False,axlabels=False):
		#---compute filters
		n_zones = 10
		def minmax(a): return a.min(),a.max()
		lipid_distance_lims = minmax(np.array([minmax(data[sn]['data']['lipid_distances']) 
			for sn in sns]).reshape(-1))
		lipid_distance_lims = (0.16,1.0)
		fenceposts = np.linspace(*lipid_distance_lims,num=n_zones)
		#---exclude some zones that have low populations for t he standard method
		if plot_line_style=='many': zones_subsel = np.array([0,1,6,7,8])
		elif plot_line_style=='single': zones_subsel = np.array([-1])
		elif plot_line_style=='gradient': 
			zones_subsel = np.array([0,1,2,3,4,5,6,7,8])
			zones_subsel = np.array([0,1,8])
			#---best to only subselect the closest zone since it already captures a lot
			zones_subsubsel = [0]
		else: raise Exception('invalid plot_line_style %s'%plot_line_style)
		zones = dict([(sn,[np.all((data[sn]['data']['lipid_distances']>=i,
			data[sn]['data']['lipid_distances']<=j),axis=0) 
			for i,j in np.array(zip(fenceposts[:-1],fenceposts[1:]))[zones_subsel]]) for sn in sns])
		view_window = scanrange[:-1]<=farcut
		#---plot the lipid distance-dependant g(r) for one simulation
		for snum,sn in enumerate(sns_ordered):
			if plot_line_style in ['many','single']:
				for znum,zone in enumerate(zones[sn]):
					curve = distributions[sn][np.concatenate(zone)].mean(axis=0)
					curve = curve/(normalizers[sn] if not normed_windows else normalize_window(curve))
					ax.plot(middles[view_window],curve[view_window],'-',lw=1,zorder=3,
						color=colormaps[sn](znum/float(len(zones[sn])-1)) if plot_line_style=='many' else 
							colormaps[sn](1))
			elif plot_line_style=='gradient':
				offset_color = 1.0
				for inum,index in enumerate(zones_subsubsel):
					zonepair = [zones[sn][zones_subsel[index+i]] for i in range(2)]
					curves = [distributions[sn][np.concatenate(z)].mean(axis=0) for z in zonepair]
					curves = [c/(normalizers[sn] if not normed_windows else normalize_window(c))
						for c in curves]
					#---fill between works pretty well
					if False:
						ax.fill_between(middles[view_window],y1=curves[1][view_window],
							y2=curves[0][view_window],
							lw=1,zorder=3,color=colormaps[sn]((inum)/
								float(len(zones_subsubsel)-1+offset_color)))
					else:
						ngrad = float(100)
						for s in range(int(ngrad))[::-1]:
							ax.plot(middles[view_window],
								(curves[1][view_window]-curves[0][view_window])*float(s)/ngrad+
								curves[0][view_window],lw=1,zorder=3,
								color=colormaps[sn](float(s)/ngrad) 
									if not do_opacity else colormaps[sn](0.0),
								alpha=((1.-s/ngrad)**4.) if do_opacity else 1.0)
			else: raise Exception('invalid plot_line_style %s'%plot_line_style)
		if xlims: ax.set_xlim(xlims)
		if ylims: ax.set_ylim(ylims)
		if one: ax.axhline(1.0,c='k',lw=0.5,zorder=1,alpha=0.35)
		ax.tick_params(axis='both',which='both',length=0)
		if axlabels:
			ax.set_xlabel('r ($\mathrm{\AA}$)')
			ax.set_ylabel('g(r)')
	#---assemble the plot
	ax = plt.subplot(111)
	xlims,ylims = (0.3,1.0),(0.,4.)
	ax.set_aspect((xlims[1]-xlims[0])/(ylims[1]-ylims[0])*1.0)
	plot(ax,xlims=xlims,ylims=ylims,plot_line_style='gradient',one=True,do_opacity=False,axlabels=True)
	axins = inset_axes(ax,width="40%",height="40%",loc=1)
	xlims,ylims = (0.3,1.0),(0.,4.)
	plot(axins,xlims=(0.1,0.7),ylims=None,plot_line_style='single',axlabels=True)
	picturesave('fig.test12224',work.plotdir,backup=False,version=True,meta={},extras=[])

###---STANDARD

def printer():
	"""Load once per plot session."""
	global variables,routine,printers
	#---reload if not all of the globals in the variables
	if any([v not in globals() or globals()[v] is None for v in variables]): reload()
	#---after loading we run the printers
	printers = list(set(printers))
	if routine is None: routine = list(printers)	
	#---routine items are function names
	for key in routine: 
		status('running routine %s'%key,tag='printer')
		globals()[key]()

if __name__=='__main__': printer()
