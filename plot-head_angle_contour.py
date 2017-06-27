#!/usr/bin/env python

if False:

	#---plot prep
	if 'plotload' not in globals(): execfile('/etc/pythonstart')
	execfile('./omni/base/header.py')
	from plotter import *
	from base.store import plotload
	execfile('./calcs/specs/figures.py')
	execfile('./calcs/specs/colors.py')
	import matplotlib.patches as mpatches
	from codes.diffusion import msd_fit,msd_fit_deprecated
	from copy import deepcopy
	import numpy as np
	ns = sys.argv[0]

	from scipy.stats import kde
	from codes import vmdmake

import scipy
import scipy.stats
try: from codes import vmdmake
except: raise Exception(
	'cannot import vmdmake. try `git clone http://github.com/bradleyrp/amx-vmd %s`'%
	os.path.join(os.getcwd(),'calcs','codes','vmdmake'))

#---! hard-coded paths
kwargs_vmdmake = {'CGBONDSPATH':'~/libs/cg_bonds.tcl','GMXDUMP':'~/libs/gmxdump'}

#---settings
plotname = 'head_angle'
routine = ['simple','head_angle_snaps','simple_detailed'][-2:-1]

if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

sns = work.vars['orders']['canon']['symmetric']+work.vars['orders']['canon']['asymmetric']
sns_special = ['membrane-v531','membrane-v532']
sns_top = [s for s in sns if s not in sns_special]

fsbase = 20
nbins = 50
nlevels = 25
xmin,xmax,ymin,ymax = -180,180,-180,180
color_map_name = 'gnuplot2 pink jet'.split()[-1]
plot_simple = False
angle_ticks = np.arange(-150,200,50)

if False:
	bands = [mpl.colors.ColorConverter().to_rgb(i) for i in [
		"gray",
		colors['blue'],
		]]
	grads,cmap = colorscale(bands=bands,return_cmap=True)
cmap = mpl.cm.get_cmap(color_map_name)

#---snapshotting positions
positions = ['modal','modal_low']

#---snapshot naming tag
tag = 'v4'

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

if 'postdat' not in globals():
	postdat = {}
	for pnum,sn in enumerate(sns):
		status('kernel density estimates %s'%sn,looplen=len(sns),i=pnum,tag='plot')
		theta,phi = [data[sn]['data'][i] for i in ['theta','phi']]
		#---HACK for backwards PHI in symmetric systems
		factor = -1.0 if work.meta[sn]['composition_name']=='symmetric' else 1.0
		phi = phi*factor
		raw = np.array([theta,phi]).T
		#---no averaging: everything goes into one bin
		catted = np.concatenate(raw)
		#---KDE
		x,y = catted.T
		k = scipy.stats.kde.gaussian_kde(catted.T)
		xi, yi = np.mgrid[xmin:xmax:nbins*1j,ymin:ymax:nbins*1j]
		zi = k(np.vstack([xi.flatten(), yi.flatten()]))
		postdat[sn] = dict(xi=xi,yi=yi,zi=zi,x=x,y=y)

if 'simple' in routine:

	zmax = max([postdat[sn]['zi'].max() for sn in sns])
	if plot_simple: axes,fig = square_tiles(len(sns),figsize=(10,10))
	else:
		#---discarding this and repeating for a better layout
		if False:
			#---logic for a sequential square of tiles
			ntiles = len(sns_top)
			favor_rows = False
			nrows = int(np.sqrt(ntiles))+1
			ncols = float(ntiles)/nrows
			ncols = int(ncols) if float(ncols)==int(ncols) else int(ncols+1)
			if not favor_rows: nrows,ncols = ncols,nrows
			axes,fig = panelplot(figsize=(10,10),
				layout={'out':{'grid':[1,1]},
				'ins':[{'grid':[nrows,ncols],'hspace':0.7,'wspace':0.7},
				{'grid':[1,1]}]})
			for i in range(nrows*ncols-len(sns_top)): fig.delaxes(axes[0][-1*(i+1)])
		#---side-by-side symmetric vs asymmetric
		axes,fig = panelplot(figsize=(16,10),
			layout={'out':{'grid':[1,2],'wspace':0.2},
			'ins':[
				{'grid':[3,2],'hspace':0.7,'wspace':0.7},
				{'grid':[3,2],'hspace':0.7,'wspace':0.7},
				]})
		def get_which_plot(sn):
			compname = work.meta[sn]['composition_name']
			return axes[
				['symmetric','asymmetric'].index(compname)][
				[s for s in [i for i in sns_top if work.meta[i]['composition_name']==compname]].index(sn)]

	for pnum,sn in enumerate(sns_top):
		status('plotting %s'%sn,looplen=len(sns),i=pnum,tag='plot')
		ax = get_which_plot(sn)
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

def get_blank_border(img):
	"""Return the border limits for removing whitespace."""
	nonblank = np.any(img!=1.0,axis=2)*1
	lims = [map(lambda x:(x[0]-1,x[-1]+1),
		[np.where(np.any(nonblank!=0,axis=j))[0]])[0] for j in range(2)]
	return lims

good_side = """
draw delete all
set me [atomselect top "all"]
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
$me moveby [vecscale -1.0 [expr [$me_c11 get {x y z}]]]
display resetview
mouse stoprotation
rotate x to -90
if { $reverse_view == 1 } { rotate y by 90 } else { rotate y by -90 }
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
set vec_ring [vecsub [expr [$me_c14 get {x y z}]] [expr [$me_c11 get {x y z}]]]
set vec_up {0.0 0.0 1.0}
set vec_look [veccross $vec_ring $vec_up]
set relevant_angle [expr [tcl::mathfunc::acos [vecdot [vecnorm $vec_ring] $vec_up]]/3.1415926*180.0]
set tmat [transvecinv $vec_look]
$me move $tmat
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
# REMOVED draw cylinder [expr [$me_c11 get {x y z}]] [expr [$me_c14 get {x y z}]] radius 0.1
"""

def snapshot_namer(**kwargs):
	"""Name the snapshots in the plot folder."""
	canon_order = ['tag','sn','pos']
	mark_out,mark_in = '.','_'
	suffix = mark_out.join(
		'%s%s%s'%(str(key),mark_in,str(val))
		for key,val in [(k,kwargs[k]) for k in canon_order])
	return 'fig.head_angle.snap.'+suffix

#---store the snapshots in the post_plot_spot according to the tag
tempdir = os.path.join(work.paths['post_plot_spot'],'fig.head_angle.%s'%tag)
if not os.path.isdir(tempdir): os.mkdir(tempdir)
status('snapshots dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')

if 'head_angle_snaps' in routine:

	"""
	The following code was cannibalized from plot-head_angle.
	"""

	def exemplar_lipid(theta,phi,x,y,nframes,nlipids):
		ind = np.argmin(np.sum((np.array((x,y)).T-np.array([theta,phi]))**2,axis=1))
		return np.unravel_index(ind,(nframes,nlipids))
	
	"""
	checking that the angles are correct per david's observation that the detail plot was wrong
	data[sn]['data']['theta'][frame][resid_rel] should equal theta but it doesn't
	this means that the exemplar_lipid function is failing
	recall that x,y are the theta,phi angles (I THINK -- CHECK THIS)
	error is obv due to mishandling at the unravel_index stage
	only four possibilities so tried all four
	data[sn]['data']['theta'][reps[(sn,pos)]['frame']][reps[(sn,pos)]['resid_rel']],reps[(sn,pos)]['theta']
	"""
	
	exemplar_rank = 1
	
	def exemplar_lipid(theta,phi,x,y,nframes,nlipids,rank=0):
		resid_rel,frame = np.unravel_index(np.argsort(
			np.sum((np.array((x,y)).T-np.array([theta,phi]))**2,axis=1))[rank],(nlipids,nframes))
		return frame,resid_rel
	
	#---pick select postions for plotting
	reps = {}
	for sn in sns_special:
		xi,yi,zi,x,y = [postdat[sn][key] for key in ['xi','yi','zi','x','y']]
		nframes,nlipids = data[sn]['data']['theta'].shape
		for pos in positions:
			if pos=='modal':
				ind = np.unravel_index(np.argmax(zi.reshape(xi.shape)),xi.shape)
				theta,phi = xi[ind],yi[ind]
			#---somewhat hackish
			elif pos=='modal_low':
				#---actually ultra hackish
				xyz = np.transpose((xi.reshape(-1),yi.reshape(-1),zi));
				theta,phi,_ = next(i for i in sorted(xyz,key=lambda a_entry:a_entry[2])[::-1]  
					if i[0]<0 and i[1]<0)
			frame,resid_rel = exemplar_lipid(theta,phi,x,y,nframes,nlipids,rank=exemplar_rank)
			reps[(sn,pos)] = {'frame':frame,'resid_rel':resid_rel,'theta':theta,'phi':phi}
			
	rendered_fns = {}
	for sn in sns_special[:]:

		#---move the cursor back
		work.cursor = ('ptdins','xtc')
		#---! previously in legacy: gro = work.slice(sn)['current']['all']['gro']
		#---modified the calculation output to include the slices for use-cases like this one
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(calc['extras'][sn]['slice_path'],suf))
			for suf in ['gro','xtc']]
		#---get the tpr from the raw data
		tpr = work.raw.get_last(sn,subtype='tpr')

		#---precompute snapshot filenames and plot if any are missing
		snapshot_fns = [snapshot_namer(sn=work.raw.prefixer(sn_this),tag=tag,pos=pos) 
			for (sn_this,pos),val in reps.items() if sn_this==sn]
		if any([not os.path.isfile(os.path.join(tempdir,i+'.png')) for i in snapshot_fns]):

			#---! CHANGED THE LOOP TO SLOW VERSION WHEN WE ADDED IONS

			#---loop over desired lipid/frame pairs for rendering
			for key in [r for r in reps if r[0]==sn]: 

				#---run VMD to make the snapshots via vmdmake
				view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
					frames='',xres=2000,yres=2000,**kwargs_vmdmake)
				view.do('load_dynamic','standard','bonder')

				val = reps[key]
				frame,resid_rel,angle = [val[i] for i in ['frame','resid_rel','theta']]
				sn,pos = key
				view['snapshot_filename'] = snapshot_namer(sn=work.raw.prefixer(sn),pos=pos,tag=tag)
				view.command("set lipids_select [atomselect top \"resname %s\"]"%
					work.meta[sn]['ptdins_resname'])
				view.command("set my_index [lindex [lsort -unique [$lipids_select get resid]] %d]"%
					resid_rel)
				headsel = '(%s)'%' or '.join(['name %s'%i for i in work.vars['head_atoms'].split()])
				view.select(**{'me':"resname %s and resid $my_index and not hydrogen and %s"%(
					work.meta[sn]['ptdins_resname'],headsel),
					'smooth':False,'style':'licorice','goodsell':True})
				view.command("animate goto %d"%frame)
				view.do('reset','xview')
				view.command('set reverse_view %d'%(0 if angle>0 else 1))
				#---! EXTREMELY IMPORTANT NOTE: DON'T MESS WITH SELECTIONS WITHOUT REMEMBERING
				#---! ...that the GOODSIDE code rotates the molecules to get the right view
				#---! ...since we cannot rotate the camera
				view.command(re.sub('PIP2RESNAME',work.meta[sn]['ptdins_resname'],good_side))
				view.set_color_cursor({'MG':'pink','Cal':'blue'}[work.meta[sn]['cation']])
				view.select(**{'near_ions':
					"name %s and within 4.6 of (resname %s and resid $my_index)"%(
					work.meta[sn]['cation'],
					work.meta[sn]['ptdins_resname']),
					'smooth':False,'style':'Beads 0.6 12.0','goodsell':True,'color_specific':True})
				#---skipping water because it wasn't drawing bonds after transformations, other problems
				if False: view.select(**{'hydration_shell':
					"water and (same residue as (within 3.01 of (name %s "%work.meta[sn]['cation']+
					"and within 4.6 of (resname %s and resid $my_index))))"%(
					work.meta[sn]['ptdins_resname']),
					'smooth':False,'style':'cpk_water','goodsell':True})
				view.command('scale by 0.9')
				view.do('snapshot')
				view.command("mol delrep 1 top")
				view.command("mol delrep 0 top")
				#---YOU MUST TERMINATE THE REPS IN THIS LOOP FOR ALL SELECTIONS!
				view.show(quit=True)

#---SEPARATE DETAILED VIEW
if 'simple_detailed' in routine:

	#---rows and columns on the left side
	nrl,ncl = 4,2
	axes,fig = panelplot(figsize=(16,12),
		layout={'out':{'grid':[2,2],'wspace':0.0,'hspace':0.3,'wratios':[1,2]},
		'ins':[
			{'grid':[2,2],'hratios':[4,1][::-1],'wratios':[1,4][::-1],'wspace':0.1,'hspace':0.1},
			{'grid':[1,2],'wspace':0},
			{'grid':[2,2],'hratios':[4,1][::-1],'wratios':[1,4][::-1],'wspace':0.1,'hspace':0.1},
			{'grid':[1,2],'wspace':0},
			]})
	axn = {'contour':1,'y':0,'x':3,'del':2}
	axn = {'contour':2,'y':3,'x':0,'del':1}

	zmax = max([postdat[sn]['zi'].max() for sn in sns_special])
	for ss,sn in enumerate(sns_special):
		pnum = [0,2][ss]
		ax = axes[pnum][axn['contour']]
		xi,yi,zi,x,y = [postdat[sn][key] for key in ['xi','yi','zi','x','y']]
		ax.pcolormesh(xi,yi,zi.reshape(xi.shape),vmin=0,vmax=zmax,cmap=cmap)
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)

		labels = {'theta':r'tilt angle $\mathrm{\theta}$',
			'phi':r'rotation angle $\mathrm{\phi}$'}
		ax.set_xlabel(labels['theta'],fontsize=fsbase)
		ax.set_ylabel(labels['phi'],fontsize=fsbase)

		cs = ax.contourf(xi,yi,zi.reshape(xi.shape),vmax=zmax,vmin=0,levels=np.linspace(0,zmax,nlevels),
			extend='both',origin='lower',lw=2,zorder=3,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
		xbins = np.arange(xmin,xmax,(xmax-xmin)/nbins)
		ybins = np.arange(ymin,ymax,(ymax-ymin)/nbins)

		angle_image_format(ax)

		ax = axes[pnum][axn['x']]
		ax.hist(x,bins=xbins,color='gray',lw=1,edgecolor='w')

		ax.spines['left'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])

		ax = axes[pnum][axn['y']]
		ax.hist(y,bins=ybins,orientation='horizontal',color='gray',lw=1,edgecolor='w')

		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.yaxis.set_ticks_position('none')
		ax.xaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])

		fig.delaxes(axes[pnum][axn['del']])
		for ax in [axes[pnum][axn[i]] for i in ['x','y','contour'][-1:]]:
			ax.set_xticks(angle_ticks)
			ax.set_yticks(angle_ticks)

	patches,counter = [],0
	tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")

	#---MUST CREATE SNAPSHOTS FIRST
	for ss,sn in enumerate(sns_special):
		pnum = [1,3][ss]
		for posnum,pos in enumerate(positions):
			fn = os.path.join(tempdir,snapshot_namer(sn=work.prefixer(sn),tag=tag,pos=pos)+'.png')
			im = mpl.image.imread(fn)
			ax = axes[pnum][posnum]
			ax.imshow(im)
			ax.axis('off')
			tb = ax.text(0.5,0,chr(ord('A')+counter),
				bbox=tagbox,rotation=0,ha="center",va="top",color='k',
				fontsize=fsbase,transform=ax.transAxes)
			patches.append(tb)
			if posnum == 0:
				tb = ax.text(0.1,0.9,work.meta[sn]['ion_label'],
					bbox=tagbox,rotation=0,ha="center",va="center",color='k',
					fontsize=fsbase+4,transform=ax.transAxes,zorder=10)
				patches.append(tb)
			#---add the marker to the histogram
			ax = axes[[0,2][ss]][axn['contour']]
			theta,phi = [reps[(sn,pos)][i] for i in ['theta','phi']]
			tb = ax.text(theta,phi,chr(ord('A')+counter),
				rotation=0,ha="center",va="center",color='k',
				fontsize=fsbase-4)
			tb.set_path_effects([path_effects.Stroke(linewidth=4,foreground='w'),
				path_effects.Normal()])
			#bb = tb.get_bbox_patch()
			#bb.set_boxstyle("round,pad=0.1")
			counter += 1

	#---saving the snapshot tag here so we can keep track of the fix to exemplar_lipid above, fixed on v4
	picturesave('fig.head_angle_detail',work.plotdir,backup=False,version=True,
		meta={'tag':tag,'exemplar_rank':exemplar_rank},extras=patches)
	plt.close()

