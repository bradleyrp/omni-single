#!/usr/bin/env python

import scipy.ndimage
import scipy.spatial
import tempfile
import shutil
import joblib

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from plotter import *

#---update for new omnicalc
from plotter.panels import panelplot

import numpy as np

def print_birdseye_snapshot_render(surfs,protpts,mvecs,nprots,handle='',
	outdir='',fi=0,pbc_expand=1.0,smooth=1.,extrema=None,metadat=None,fs=None,
	titles=None,panelspecs=None,cmap_mpl_name='RdBu_r',fn=None,**kwargs):
	"""
	Outer loop for rendering protein videos.
	"""
	if not all([len(np.shape(s))==2 for s in surfs]):
		raise Exception('incoming data must be simulation by x by y')
	npanels = len(surfs)
	axes,fig = panelplot(**panelspecs)
	#---loop over panels
	for ni in range(npanels):
		m,n = griddims = np.shape(surfs[ni])
		grid_spacing = mvecs[ni][:2]/griddims
		#---smoothing parameter is in nm
		grid_smooth = smooth/np.mean(grid_spacing,axis=0)
		border = [int(float(pbc_expand)*i) for i in (m,n)]
		griddat = np.tile(surfs[ni],(3,3))[m-border[0]:2*m+border[0],n-border[1]:2*n+border[1]]
		extent = [-(border*grid_spacing)[0],mvecs[ni][0]+(border*grid_spacing)[0],
			-(border*grid_spacing)[1],mvecs[ni][1]+(border*grid_spacing)[1]]
		ax = axes[ni]
		#---plot the surface
		im = ax.imshow((scipy.ndimage.filters.gaussian_filter(griddat,grid_smooth)).T,
			interpolation='nearest',origin='lower',cmap=mpl.cm.__dict__[cmap_mpl_name],
			vmax=extrema[1],vmin=extrema[0],extent=extent)
		if pbc_expand != 0.:
			ax.plot([0,mvecs[ni][0],mvecs[ni][0],0,0],[0,0,mvecs[ni][1],mvecs[ni][1],0],'-',lw=2,c='k')
		if type(protpts)!=type(None):
			plothull(ax,protpts[ni],vecs=None,griddims=griddims,fill=False,lw=2)
		ax.set_xlabel(r'$\mathrm{x\,(nm)}$',fontsize=fs['axlabel'])
		ax.set_ylabel(r'$\mathrm{y\,(nm)}$',fontsize=fs['axlabel'])
		ax.tick_params(axis='x',which='both',bottom='off',top='off',
			labelbottom='on',labelsize=fs['axlabel'])
		ax.tick_params(axis='y',which='both',left='off',right='off',
			labelbottom='off',labelsize=fs['axlabel'])
		ax.set_xlim(extent[0],extent[1])
		ax.set_ylim(extent[2],extent[3])
		if titles: ax.set_title(titles[ni])
	for i in range(len(axes)-npanels): fig.delaxes(axes[-1])
	#---colorbar for the last axis	
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.06,0.,1.,1.),
		bbox_transform=ax.transAxes,
		borderpad=0)
	cbar_kwargs = {}
	nzticks = kwargs.get('nzticks',None)
	if nzticks: 
		zticks_precision = kwargs.get('zticks_precision',1)
		ins_ticks = np.arange(0,extrema[1],extrema[1]/nzticks).round(zticks_precision)
		ins_ticks = np.concatenate((ins_ticks[1:][::-1]*-1,ins_ticks))
		cbar_kwargs['ticks'] = ins_ticks
	cbar = plt.colorbar(im,cax=axins,orientation="vertical",**cbar_kwargs)
	cbar.ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
	zlabel = kwargs.get('zlabel',r'$\mathrm{z\,(nm)}$')
	cbar_label_specs = kwargs.get('cbar_label_specs',dict(rotation=270,labelpad=20))
	cbar.ax.set_ylabel(zlabel,**cbar_label_specs)
	#---print
	if outdir != '' or fn:
		plt.savefig(outdir+'/'+fn+'.png',
			dpi=200,bbox_inches='tight')
		plt.close(fig)
	else: plt.show()

def plothull(ax,points,griddims=None,vecs=None,c=None,mset=None,
	alpha=None,nine=False,fill=None,radius=0.5,ec=None,lw=1):
	"""
	Plot an outline around each protein.
	This function was copied from elsewhere.
	Copied in from simuluxe and reworked for arbitrary protein points dimensions.
	"""
	if alpha == None: alpha = 0.65
	if fill == None: fill = True
	if type(vecs)==type(None): vecs = griddims
	m,n = griddims
	#import ipdb;ipdb.set_trace()
	#pts = np.array([[i[0]*m/vecs[0],i[1]*n/vecs[1]] for i in points])
	patches = []
	#---you must send a vector with three dimensions: protein index, atom, xy
	for pts in points:
		hull = scipy.spatial.ConvexHull(pts)
		p = ax.add_patch(mpl.patches.Polygon(pts[hull.vertices],
			closed=True,fill=(True if fill else False),
			facecolor=c,lw=lw,alpha=alpha,edgecolor=(c if ec == None else ec)))
		patches.append(p)
	return patches
