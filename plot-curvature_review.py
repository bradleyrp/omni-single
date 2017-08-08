#!/usr/bin/env python

"""
INSPECT THE CURVATURE COUPLING RESULTS
"""

import scipy
import scipy.interpolate

deg2rad = lambda d: d/180.*np.pi
if is_live: 
	from ipywidgets import *
	from IPython.display import display

#---block: get the curvature investigation class
#---! this is pretty much redundant!
import codes.curvature_coupling.InvestigateCurvature
codes.curvature_coupling.InvestigateCurvature.plotload = plotload
codes.curvature_coupling.InvestigateCurvature.plotname = plotname
codes.curvature_coupling.InvestigateCurvature.work = work
if 'ic' not in globals(): ic = codes.curvature_coupling.InvestigateCurvature.InvestigateCurvature(mode='plot')

#---block: organize the post-processing data
class PostDat:
	"""Organize the post-processing data for brevity in the plotters."""
	def __init__(self):
		global work
		self.sns = work.sns()
		self.isotropy_sweep = ic.spec.get('isotropy',[1.0])
		self.theta_sweep = ic.spec.get('theta',[1.0])
		self.data,self.toc,self.range = [],[],[]
		#---double loop gets C_0,sigma_a landscapes over isotropy and theta
		for isotropy in self.isotropy_sweep:
			for theta in self.theta_sweep:
				spec = dict(isotropy=isotropy,theta=deg2rad(theta))
				self.data.append(ic.generate_standard_manual_landscapes(**spec))
				self.range.append([f([f(self.data[-1][sn].reshape(-1)) for i in self.isotropy_sweep 
					for t in self.theta_sweep for sn in self.sns]) for f in [min,max]])
				self.toc.append(dict(spec,theta=theta))
		#---get the master range and histograms
		self.range_master = [getattr(np.array(self.range),k)() for k in ['min','max']]
		self.counts = [dict([(sn,np.histogram(data[sn].reshape(-1),
			bins=np.linspace(self.range_master[0],self.range_master[1],100)))
			for sn in self.sns]) for data in self.data]
		self.maxcount = max([counts[sn][0].max() for counts in self.counts for sn in self.sns])
		#---prepare contour plots
		contour_interp_pts = 100
		self.curvature_sweep = ic.spec.get('C_0',[1.0])
		self.extent_sweep = ic.spec.get('sigma_a',[1.0])
		c0,c1 = min(self.curvature_sweep),max(self.curvature_sweep)
		e0,e1 = min(self.extent_sweep),max(self.extent_sweep)
		finex,finey = np.linspace(c0,c1,contour_interp_pts),np.linspace(e0,e1,contour_interp_pts)
		X,Y = np.meshgrid(finex,finey)
		self.contour_xy = X,Y
		reworked_lands = [dict([(sn,np.array([(c,e,land[sn][cc,ee]) 
			for cc,c in enumerate(self.curvature_sweep) for ee,e in enumerate(self.extent_sweep)])) 
			for sn in self.sns]) for land in self.data]
		self.contours = [dict([(sn,scipy.interpolate.griddata(
			land[sn][:,:2],land[sn][:,2],(X,Y),method='cubic')) for sn in self.sns]) 
			for land in reworked_lands]
	def get(self,name='data',**kwargs):
		these = [ii for ii,i in enumerate(self.toc) if kwargs==i]
		if len(these)!=1: raise Exception('lookup failure: %s'%([self.data[ii] for ii in these]))
		else: index = these[0]
		return self.__dict__[name][these[0]]
if 'post' in globals(): del post
if 'post' not in globals(): post = PostDat()

#---block: plot the landscapes
def landscaper(isotropy=1.0,theta=0.0):
	"""Plot the landscapes."""
	nrows = 3
	axes,fig = panelplot(figsize=(16,10),layout={
		'out':{'grid':[nrows,1]},'ins':[{'grid':[1,len(post.sns)]} for r in range(nrows)]})
	lands,(vmin,vmax),counts,contours = [post.get(name=name,isotropy=isotropy,theta=theta) 
		for name in ['data','range','counts','contours']]
	for snum,sn in enumerate(post.sns):
		#---raw error landscape
		ax = axes[0][snum]
		ax.imshow(lands[sn].T,
			interpolation='nearest',origin='lower',
			cmap=mpl.cm.__dict__['gist_stern_r'],vmin=vmin,vmax=vmin+(vmax-vmin)*0.1)
		ax.set_title(work.meta[sn].get('name',sn))
		#---histogram of the errors
		ax = axes[1][snum]
		ax.plot(counts[sn][1][:-1],counts[sn][0],lw=0,marker='.',ms=10,color='r')
		for count in post.counts: 
			ax.plot(count[sn][1][:-1],count[sn][0],lw=0,marker='.',color='k',alpha=0.1)
		ax.set_yscale('log')
		ax.set_ylim((1,post.maxcount*1.05))
		#---contour plots
		ax = axes[2][snum]
		ax.contourf(post.contour_xy[0],post.contour_xy[1],contours[sn])
	plt.show()

#---block: interactive plot
if is_live:
	slider_isotropy = widgets.SelectionSlider(options=post.isotropy_sweep,
		description='isotropy',disabled=False,continuous_update=False,readout=True)
	slider_theta = widgets.SelectionSlider(options=post.theta_sweep,
		description='theta',disabled=False,continuous_update=False,readout=True)
	interact(landscaper,isotropy=slider_isotropy,theta=slider_theta);
>>>>>>> 7467280eceaa8ba989bf424a168b4681d5645ee7
