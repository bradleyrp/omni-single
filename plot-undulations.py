#!/usr/bin/env python

"""
Plot bilayer areas for the ptdins project.
"""

from codes.undulate import calculate_undulations
from codes.undulate_plot import undulation_panel,add_undulation_labels,add_axgrid,add_std_legend

@autoload(plotrun)
def load():
	"""Load the data."""
	#---condition for loading the data
	if 'data' not in globals(): 
		#---load once
		global data,calc,data_prot
		data,calc = plotload('undulations',work)
		data_prot,_ = plotload('protein_abstractor',work)

@autoplot(plotrun)
def plot_height_profiles():
	"""
	Plot the bilayer height profile with protein positions.
	"""
	#---one plot per simulation
	for sn in work.sns():
		mesh = data[sn]['data']['mesh']
		surf = mesh.mean(axis=0).mean(axis=0)
		surf -= surf.mean()
		hmax = np.abs(surf).max()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		vecs = data[sn]['data']['vecs'].mean(axis=0)
		im = ax.imshow(surf.T,origin='lower',
			interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
			extent=[0,vecs[0],0,vecs[1]],vmax=hmax,vmin=-1*hmax)
		if sn in data_prot:
			from render.wavevids import plothull
			points_all = data_prot[sn]['data']['points_all'] 
			points_all_mean_time = points_all.mean(axis=0)
			plothull(ax,points_all_mean_time[...,:2],griddims=surf.shape,vecs=vecs,c='k',lw=0)	
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		cbar.set_label(r"$\mathrm{\langle z \rangle (nm)}$",labelpad=-40,y=1.05,rotation=0)
		axins.tick_params(axis='y',which='both',left='off',right='off',labelright='on')
		ax.set_title('average bilayer height')
		ax.set_xlabel('x (nm)')
		ax.set_ylabel('y (nm)')
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		picturesave('fig.height_average.%s'%sn,work.plotdir,backup=False,version=True,meta={})

@autoplot(plotrun)
def undulation_spectra(style='all_on_one'):
	"""
	Manage different plot combinations.
	This is currently under development, and needs to be resolved with plot-analyze_dextran_undulations.py.
	"""
	#---external settings
	plotspecs = work.plots[plotname].get('specs',{})
	wavevector_limits = plotspecs.get('wavevector_limits',[1.0])
	#---top level sweep over plot settings and styles
	sweep_specs = {
		'wavevector_limit':wavevector_limits,
		'style':['all_on_one'],
		#---! need custom heights for average_normal
		'midplane_method':['flat','average','average_normal'][:-1],}
	#---get colors first from art which always
	colors_reqs = ['binned','fitted','line']
	if colors!=None: 
		#---must match the receiver in the plot function
		sweep_specs['color'] = [dict([(sn,dict([(key,colors[sn]) 
			for key in colors_reqs])) for sn in colors])]
	#---we could also get the colors systematically from the metadata here?
	#---fallback randomly generates random colors
	else: 
		def random_colormaker(x,n,name='jet'): 
			return mpl.cm.__dict__[name](np.linspace(0.,1.0,len(work.sns()))[x]/n)
		colors_random = dict([(sn,dict([(key,random_colormaker(snum,len(work.sns())))
			for key in colors_reqs])) for snum,sn in enumerate(work.sns())])
		sweep_specs['color'] = [colors_random]
	plots_this = sweeper(**sweep_specs)
	#---settings
	art = {'fs':{'legend':12}}
	#---loop over all plot settings and styles
	for pnum,plotspec in enumerate(plots_this):
		#---unpack settings
		style = plotspec['style']
		#---router over plot layouts
		if style=='all_on_one':
			figsize = (5,5)
			layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
			axes,fig = panelplot(layout,figsize=figsize)
			ax = axes[0]
			for snum,sn in enumerate(work.sns()):
				plot_undulation_spectrum(ax,sn,**plotspec)
			decorate_undulation_plot(ax=ax,art=art)
		else: raise Exception('invalid plot style %s'%style)
		#---save the plot, with relevant keys
		meta_reqs = ['midplane_method','wavevector_limit','style']
		meta = dict([(key,plotspec[key]) for key in meta_reqs])
		picturesave('fig.undulations',work.plotdir,backup=False,version=True,meta=meta)

def plot_undulation_spectrum(ax,sn,**kwargs):
	"""
	Plot a single undulation spectrum.
	"""
	mesh = data[sn]['data']['mesh']
	surf = mesh.mean(axis=0)
	surf = np.array([s-surf.reshape(len(surf),-1).mean(axis=1)[ss] for ss,s in enumerate(surf)])
	vecs = data[sn]['data']['vecs']
	fit_style = kwargs.get('fit_style','band,perfect,curvefit')
	lims = kwargs.get('lims',[0.,kwargs.get('wavevector_limit',1.0)])
	colors = kwargs.get('color','k')
	midplane_method = kwargs.get('midplane_method','flat')
	colors_reqs = ['binned','fitted','line']
	if type(colors) in str_types: 
		colors = dict([(key,colors) for key in colors_reqs])
	elif not all([key in colors.get(sn,{}) for key in colors_reqs]): 
		raise Exception('need keys in colors %s'%colors_reqs)
	uspec = calculate_undulations(surf,vecs,fit_style=fit_style,lims=lims,
		midplane_method=midplane_method,fit_tension=kwargs.get('fit_tension',False))
	label = work.meta[sn]['label']+'\n'+r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
	q_binned,energy_binned = uspec['q_binned'][1:],uspec['energy_binned'][1:]
	ax.plot(q_binned,energy_binned,'.',lw=0,markersize=10,markeredgewidth=0,
		label=None,alpha=0.2,color=colors[sn]['binned'])
	q_fit,energy_fit = np.transpose(uspec['points'])
	ax.plot(q_fit,energy_fit,'.',lw=0,markersize=4,markeredgewidth=0,
		label=label,alpha=1.,zorder=4,color=colors[sn]['fitted'])
	def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
		return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
	ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
		area=uspec['area']),lw=1,zorder=3,color=colors[sn]['line'])

def decorate_undulation_plot(ax,art):
	"""
	Uniform appearance for all undulation spectra.
	"""
	add_undulation_labels(ax,art=art)
	add_std_legend(ax,loc='upper right',art=art)
	add_axgrid(ax,art=art)
