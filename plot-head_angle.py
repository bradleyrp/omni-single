#!/usr/bin/env python

from scipy.stats import kde

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns,(data,calc) = work.sns(),plotload(plotname,work)

#---block: settings
fsbase = 20
nbins = 50
nlevels = 25
xmin,xmax,ymin,ymax = -180,180,-180,180
color_map_name = 'gnuplot2 pink jet'.split()[-1]
plot_simple = True
angle_ticks = np.arange(-150,200,50)
routine = ['simple']
cmap = mpl.cm.get_cmap(color_map_name)

#---block: utility functions
def angle_image_format(ax,sn=None):
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	plt.setp(ax.xaxis.get_majorticklabels(),rotation=90,fontsize=fsbase-2)
	plt.setp(ax.yaxis.get_majorticklabels(),rotation=0,fontsize=fsbase-2)
	ax.axhline(0,lw=1.5,c='w',zorder=4)
	ax.axvline(0,lw=1.5,c='w',zorder=4)
	if sn:
		ax_side = ax.twinx()
		ax_side.set_yticks([])
		ax_side.set_ylabel(work.meta[sn]['composition_name'],fontsize=fsbase,rotation=270,labelpad=25)

#---block: add the kernel density estimates
if 'postdat' not in globals():
	postdat = {}
	for pnum,sn in enumerate(sns):
		status('kernel density estimates %s'%sn,looplen=len(sns),i=pnum,tag='plot')
		theta,phi = [data[sn]['data'][i] for i in ['theta','phi']]
		#---HACK for backwards PHI in symmetric systems
		#---! deprecated: factor = -1.0 if work.meta[sn]['composition_name']=='symmetric' else 1.0
		#---! ...however this needs checked!
		factor = 1.0
		phi = phi*factor
		raw = np.array([theta,phi]).T
		#---no averaging: everything goes into one bin
		catted = np.concatenate(raw)
		#---KDE
		x,y = catted.T
		k = kde.gaussian_kde(catted.T)
		xi, yi = np.mgrid[xmin:xmax:nbins*1j,ymin:ymax:nbins*1j]
		zi = k(np.vstack([xi.flatten(), yi.flatten()]))
		postdat[sn] = dict(xi=xi,yi=yi,zi=zi,x=x,y=y)

#---block: make a thumbnail plot
if 'simple' in routine:

	zmax = max([postdat[sn]['zi'].max() for sn in sns])
	axes,fig = square_tiles(len(sns),figsize=(10,10))

	for pnum,sn in enumerate(sns):
		status('plotting %s'%sn,looplen=len(sns),i=pnum,tag='plot')
		#---! previously used a custom function for getting the right axis
		ax = axes[pnum]
		xi,yi,zi = [postdat[sn][key] for key in ['xi','yi','zi']]
		ax.pcolormesh(xi,yi,zi.reshape(xi.shape),vmin=0,vmax=zmax,cmap=cmap)
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)
		cs = ax.contourf(xi,yi,zi.reshape(xi.shape),vmax=zmax,vmin=0,levels=np.linspace(0,zmax,nlevels),
			extend='both',origin='lower',lw=2,zorder=3,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
		#im = ax.imshow(xi,yi,zi.reshape(xi.shape),vmax=zmax,vmin=0,origin='lower',lw=2,zorder=2)
		ax.set_aspect('equal')
		ax.set_title('%s, %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label']))
		angle_image_format(ax,sn)
		ax.set_xlabel(r'$\mathrm{\theta}$',fontsize=fsbase-2,labelpad=-15)
		ax.set_ylabel(r'$\mathrm{\phi}$',fontsize=fsbase-2,labelpad=-15)

	picturesave('fig.head_angle_thumbs',work.plotdir,backup=False,version=False,meta={})
	plt.close()
