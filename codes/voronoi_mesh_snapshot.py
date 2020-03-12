#!/usr/bin/env python

import numpy as np
import scipy
import scipy.spatial
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_snapshots(sn,data,ax=None,mn=None):

	custom = False
	if (ax!=None and type(mn)==type(None)) or (ax==None and type(mn)!=type(None)): 
		raise Exception('define ax and monolayer together')
	elif ax!=None: custom = True
	# these are VMD colors but the only difference was DOPE is blue2
	lipid_colors = {'ptdins':'purple','PI2P':'purple','SAPI':'purple',
		'P35P':'purple','CHL1':'gray','DOPS':'red','DOPE':'blue','POPC':'green'}
	# colors above clash too much
	pip2_color = '#ff0700'
	lipid_colors = {'ptdins':pip2_color,'PI2P':pip2_color,'SAPI':pip2_color,
		'P35P':pip2_color,'CHL1':'green','DOPS':'#64aa2b','DOPE':'#2c17b1','POPC':'gray'}

	show_coms = False
	fr_specific = 500
	show_voronoi_lines = False
	vor_edge_color = 'w'
	vor_edge_lw = 0.4
	vor_alpha = 0.7
	vor_alpha_image = 0.45
	margin = 0.2
	delaunay_color = 'k'
	show_delaunay_lines = False
	do_simple_voronoi = False
	do_voronoi_pbc = True
	n_dims = 2

	if not ax:
		figsize = (8,8)
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,1],'wspace':0.3},'ins':[{'grid':[2,1],'wspace':0.0,'hspace':0.0}]})
		extras = []
	else: axes,mn = [ax],[mn]
	for ax,mn in zip(axes,[0,1]):
		ax = axes[mn]
		fr = fr_specific
		try: dat = data.this[sn]
		except: 
			print('???')
			import ipdb;ipdb.set_trace()
		vec = dat['%d.%d.%s'%(mn,fr,'vec')]
		lims = np.concatenate([(v*-1*margin,v*(1.+margin)) for v in vec][:n_dims])
		npts = sum(data.this[sn]['monolayer_indices']==mn)
		pts = np.array(dat['%d.%d.%s'%(mn,fr,'points')])
		# get simplices
		simps = np.array(dat['%d.%d.%s'%(mn,fr,'simplices')])
		simps = simps[:,np.array([2,1,0])]
		# turn simplices into lines by index
		pairs = np.concatenate([np.transpose([simps[:,i%3],simps[:,(i+1)%3]]) for i in range(3)])
		# sort the indices to ensure that they are uniform (i.e. direction is irrelevant) and take uniques
		# note that this step requires np version 1.13 or higher
		pairs = np.array(np.unique(np.sort(pairs),axis=0))
		lines_ready = lines = np.array(pts[...,:2][pairs])
		# get colors
		colors = np.array([lipid_colors[j] for j in dat['resnames'][np.where(dat['monolayer_indices']==mn)]])
		if show_coms: ax.scatter(pts[:,0],pts[:,1],
			color=colors[np.array(dat['%d.%d.%s'%(mn,fr,'ghost_ids')])],s=5,zorder=5,alpha=1.0)
		if show_delaunay_lines:
			# ensure the Delaunay mesh expandes beyond PBCs for aesthetic reasons
			lines_1d = np.array(lines.reshape((-1,4)))
			# get PBC shifts
			shifts = np.array([(i,j) for i in range(-1,2) for j in range(-1,2)])*vec[:2]
			# double the shifts to account for two points per line
			shifts_dub = np.tile(shifts,(1,2))
			shifts_expanded = np.tile(shifts_dub,(len(lines_1d),1,1)).reshape((-1,4))
			lines_expanded = np.tile(lines_1d,(9,1))
			lines_pbc = np.array(shifts_expanded + lines_expanded)
			# the format here is x,y,x,y and we filter so they are inside the box
			#! ...!!!
			ins = np.where(~np.any((
				np.all(lines_pbc[:,np.array([0,2])]<vec[0]*(-1.*margin),axis=1),
				np.all(lines_pbc[:,np.array([0,2])]>vec[0]*(1.+margin),axis=1),
				np.all(lines_pbc[:,np.array([1,3])]<vec[1]*(-1.*margin),axis=1),
				np.all(lines_pbc[:,np.array([1,3])]>vec[1]*(1.+margin),axis=1),
				),axis=0))[0]
			lines_ready = lines_pbc[ins].reshape((-1,2,2))
			#! something is severely wrong in this sequence for some frames e.g. 101: 
			#! ... see too small: print(np.unique(lines_ready,axis=0).shape)
			# collection of Delaunay lines
			lc = mpl.collections.LineCollection(lines_ready,colors=delaunay_color,linewidths=0.5,alpha=0.5)
			ax.add_collection(lc)		
		# compute voronoi
		lines_voronoi = []
		vmap = scipy.spatial.Voronoi(pts[:,:2])
		# preserving the simple method because it was so much easier
		if do_simple_voronoi:
			# use npts to avoid ghost points
			for p in range(npts):
				vertices = [vmap.vertices[i] for i in vmap.regions[vmap.point_region[p]]]
				#! fc = mpl.cm.__dict__['jet'](random.random())
				p = mpl.patches.Polygon(vertices,alpha=vor_alpha,
					facecolor=colors[p],#('w' if r >= len(colorcodes[row]) else colorcodes[row][r]),
					lw=vor_edge_lw,edgecolor=vor_edge_color)
				ax.add_patch(p)
				# save lines for later
				pairs = zip(vertices, vertices[1:]+vertices[0:1])
				lines_voronoi.extend(pairs)
			if show_voronoi_lines:
				lcv = mpl.collections.LineCollection(lines_voronoi,colors='k',linewidths=0.5)
				ax.add_collection(lcv)
		elif do_voronoi_pbc:
			# create a regular data structure containing the points for each cell
			vertices_core = [np.array([vmap.vertices[i] 
				for i in vmap.regions[vmap.point_region[p]]]) 
				for p in range(npts)]
			n_objs = len(vertices_core)
			# get the maximum constituents (number of points) in each object in the list (a cell)
			n_cons = max([len(i) for i in vertices_core]) +2
			# construct a flat n_objects by (dimension * max constituents) object
			vertices_raw = np.ones((n_objs,n_dims*n_cons))*np.inf
			# load the vertices into the flat object, reshaping over dimensions
			for ii,i in enumerate(vertices_core): vertices_raw[ii][:len(i)*2] = i.reshape(-1)
			# recall that reshape automatically does rows then columns
			# get a single set of vector shifts
			shifts = np.array([(i,j) for i in range(-1,2) for j in range(-1,2)])*vec[:2]
			# duplicate each row so it matches the shifts
			vertices_raw_pbc = np.repeat(vertices_raw,np.ones(n_objs).astype(int)*len(shifts),axis=0)
			# at this point each row is a list of constituents that need a single shift
			# so we tile the shifts in columns to match the constituents, and in rows to match the objects
			shifts_pbc = np.tile(shifts,(n_objs,n_cons))
			vertices_shifted = vertices_raw_pbc + shifts_pbc
			# now we must filter the cells by taking all that are not totally outside the box
			# totally outside the box means that any dimension is outside the box
			# this would mean that all points have any outside
			# outside means any high or low
			# lims is x low, x high, y low y high
			# dim 0 (x) and less: np.__dict__['less'](row[np.arange(0,n_dims*n_cons,n_dims)],lims[0])
			# automatic the command above: np.array([np.__dict__[k](row[
			# ... np.arange(d,n_dims*n_cons,n_dims)],lims[kk]) for d in range(2) 
			# ... for kk,k in enumerate(['less','greater'])])
			# shape above is len(lims) by n_obs by n_cons
			# since the points are paired for 2D we can add a condition that checks if the first item
			# ... in each dimension is infinity
			# check valid: np.isnan(vertices_shifted[:,np.arange(0,n_cons*n_dims,n_dims)])
			conditions = np.array([
				np.__dict__[k](
					vertices_shifted[:,np.arange(d,n_dims*n_cons,n_dims)],
					lims[kk]) 
				for d in range(2) 
				for kk,k in enumerate(['less','greater'])])
			condition_not_inside = (~np.all([np.any([conditions[kk+d*n_dims] for kk,k in enumerate(['less','greater'])],axis=0) for d in range(n_dims)],axis=0))
			condition_exists = ~np.isinf(vertices_shifted[:,np.arange(0,n_cons*n_dims,n_dims)])
			valid = (~np.all(~np.all((condition_not_inside,condition_exists),axis=0),axis=1))
			# split into core regions and image regions
			is_core = np.all(shifts_pbc==0,axis=1)
			regions_inds = np.where(valid)[0]
			regions = vertices_shifted[regions_inds]
			# break out colors by shifts
			#colors = np.tile(colors,(len(shifts),1)).reshape(-1)
			colors_expanded = np.repeat(colors,np.ones(len(colors)).astype(int)*len(shifts))
			# now we are ready to plot
			for rr,region in zip(regions_inds,regions):
				r = region[np.where(~np.isinf(region))]
				p = mpl.patches.Polygon(r.reshape((len(r)/2,2)),alpha=vor_alpha if is_core[rr] else vor_alpha_image,
					facecolor=colors_expanded[rr],#('w' if r >= len(colorcodes[row]) else colorcodes[row][r]),
					lw=vor_edge_lw,edgecolor=vor_edge_color)
				ax.add_patch(p)
		# set limits
		ax.set_xlim(-margin*vec[0],(1.+margin)*vec[0])
		ax.set_ylim(-margin*vec[1],(1.+margin)*vec[1])
		ax.set_aspect('equal')		
		# draw the pbc box
		corners = [(0,0),(0,vec[1]),(vec[0],vec[1]),(vec[0],0)]
		ax.add_collection(
			mpl.collections.LineCollection(
				[(corners[i%4],corners[(i+1)%4]) for i in range(4)],
				colors='k',linewidths=1.0))
		# aesthetics
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')		
		# vertical orientation
		ax.set_xlabel('x (nm)',fontsize=14)
		ax.set_ylabel('y (nm)',fontsize=14)

	if not custom: picturesave('fig.mesh_voronoi.%s'%sn,work.plotdir,backup=False,extras=extras)
