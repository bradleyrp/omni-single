#!/usr/bin/env python

import scipy
import scipy.ndimage

def clone_vmdmake():
	"""Clone a copy of VMDMAKE."""
	from config import bash
	vmdmake_dn = 'calcs/codes/vmdmake'
	vmdmake_source = 'https://github.com/bradleyrp/amx-vmd'
	try: from codes import vmdmake
	except: 
		if os.path.isdir(vmdmake_dn): 
			raise Exception('cannot import vmdmake from codes but %s exists'%vmdmake_dn)
		bash('git clone %s %s'%(vmdmake_source,vmdmake_dn))
		from codes import vmdmake
	globals()['vmdmake'] = vmdmake

def get_snapshots_folder(supername,overwrite_snaps=False,**kwargs):
	"""Generate a snapshot folder from a method name."""
	# store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],supername)
	# only make figures if the folder is empty
	if (os.path.isdir(tempdir) and not overwrite_snaps): 
		status('found %s and refusing to overwrite'%tempdir,tag='warning')
		return None
	try: os.mkdir(tempdir)
	except: status('directory exists so we are overwriting',tag='warning')
	status('snapshots dropping to %s'%tempdir,tag='note')
	return tempdir

@autoload(plotrun)
def load():
	sns_contacts,(data_contacts,calc_contacts) = work.sns(),plotload('contacts',work)
	clone_vmdmake()
	do_goodsell = True
	licorice_thick = ['licorice','Licorice 2.0 12.0 12.0'][-1]
	lipid_colors = {'ptdins':'purple','PI2P':'purple','SAPI':'purple',
		'P35P':'purple','CHL1':'orange','DOPS':'green3','DOPE':'blue2'}

def snapshot_namer(group,sn,frame): 
	return 'snap.overhead.%s.%s.frame_%d'%(group_name,sn,frame)

def make_overhead_snapshot(group_name,**kwargs):
	"""Make a snapshot of the overhead view."""
	sns_this = kwargs.pop('sns')
	snaps_drop_dn = kwargs.pop('snapshot_drop')
	if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
	tempdir = get_snapshots_folder(snaps_drop_dn,overwrite_snaps=True)
	for sn_this in sns_this:
		slice_path = calc_contacts['extras'][sn_this]['slice_path']
		if not hasattr(work,'source'): work.parse_sources()
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(slice_path,i)) for i in ['gro','xtc']]
		tpr = work.source.get_last(sn_this,subtype='tpr')
		kwargs_vmdmake = {'CGBONDSPATH':'/home/share/libs/cg_bonds.tcl',
			'GMXDUMP':'/usr/local/gromacs/bin/gmx'}
		frame = 500
		view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
			frames='',res=(2000,2000),**kwargs_vmdmake)
		view.do('load_dynamic','standard','bonder')
		view['snapshot_filename'] = snapshot_namer(sn=sn_this,frame=frame,
			group=group_name)
		view.set_color_cursor('black')
		view.select(**{'protein':'noh and protein','style':'cartoon','structure_color':False,
			'goodsell':do_goodsell,'smooth':False,'xy':True,'color_specific':True})
		for resname in lipid_colors.keys():
			view.set_color_cursor(lipid_colors[resname])
			view.select(**{'residue_lipid_%s'%resname:'noh and resname %s'%resname,
				'goodsell':do_goodsell,'smooth':False,'style':licorice_thick,
				'color_specific':True,'xy':True})
		view.command("animate goto %d"%frame)
		view.do('reset','zview')
		view.command('scale by 0.8')
		view.do('snapshot')
		view.show(quit=True)

@autoplot(plotrun)
def snapshots():
	"""Make many snapshots."""

	# ongoing
	snapshot_batches = {
		'v1':{'sns':sns_contacts,'snapshot_drop':'snapshots_overhead'},}

	# render snapshots and summaries
	if 0:
		for key,val in snapshot_batches.items(): 
			make_overhead_snapshot(group_name=key,**val)

	# combine overhead snapshots
	if 0:

		replicate_mapping = [('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
			('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
			('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
			('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])]

		#! hardcoded
		sns = sns_contacts
		ncols = 1
		group_name = 'v1'
		#! frame is fixed
		frame = 500
		snaps_drop_dn = 'snapshots_overhead'

		# plot panels together
		figsize = (6,12)
		nrows,ncols = 4,2
		axes,fig = panelplot(
			layout={'out':{'hspace':0.5,'wspace':0.5,'grid':[1,ncols]},
			'ins':[{'hspace':0.0,'wspace':0.0,'grid':[nrows,1]}
			for i in range(len(sns))]},figsize=figsize)
		for rnum,(name,sns_those) in enumerate(replicate_mapping):
			for snum,sn in enumerate(sns_those):
				ax = axes[snum][rnum]
				tempdir = get_snapshots_folder(snaps_drop_dn,overwrite_snaps=True)
				fn = os.path.join(work.plotdir,tempdir,
					snapshot_namer(sn=sn,group=group_name,frame=frame)+'.png')
				image = scipy.ndimage.imread(fn)
				ax.imshow(image)
				ax.axis('off')
				ax.set_title(work.meta[sn]['label'])
		picturesave('fig.snapshots_overhead.%s'%group_name,directory=work.plotdir,meta={})

	# nice side view of one simulation
	if 1:

		# one target for the nice rendering
		sn_this = 'mdia2bilayerphys'
		snaps_drop_dn = 'snapshots_side'
		tempdir = get_snapshots_folder(snaps_drop_dn,overwrite_snaps=True)
		licorice_thick = 'Licorice 0.5 12.0 12.0'
		lipid_colors_full = dict(POPC='gray',**lipid_colors)
		alt_views = 'transparent brushedmetal diffuse ghost glass1 glass2 glass3 glossy hardplastic metallicpastel steel translucent edgy edgyshiny edgyglass aoshiny aochalky aoedgy blownglass glassbubble rtchrome'.split()
		# after doing the survey it looks best with edgy on the outside and maybe 
		standard_view = ['goodsell','aochalky'][-1]
		alt_views = ['aochalky']
		xres,yres = 2000,1400
		frame = 300
		zoom = 2.0

		slice_path = calc_contacts['extras'][sn_this]['slice_path']
		if not hasattr(work,'source'): work.parse_sources()
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(slice_path,i)) for i in ['gro','xtc']]
		tpr = work.source.get_last(sn_this,subtype='tpr')
		kwargs_vmdmake = {'CGBONDSPATH':'/home/share/libs/cg_bonds.tcl',
			'GMXDUMP':'/usr/local/gromacs/bin/gmx'}

		# loop over alternative view on the sides for comparison representations
		for aa,alt_view in enumerate(alt_views):
			view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
				frames='',res=(xres,yres),**kwargs_vmdmake)
			view.do('load_dynamic','standard','bonder')
			extra_tag = '.rep_side_%s'%alt_views if len(alt_views)>1 else ''
			view['snapshot_filename'] = 'snap.side.v1.%s.frame_%d%s'%(sn_this,frame,extra_tag)
			view.set_color_cursor('black')
			#! note that we do not plot the proteins in the images so be careful to zoom all the way
			view.select(**{'protein':'noh and protein','style':licorice_thick,'structure_color':False,
				'smooth':False,'ystrip':False,'goodsell':True})
			view.set_color_cursor('black')		
			view.select(**{'protein_cartoon':'noh and protein','style':'cartoon','structure_color':False,
				'smooth':False,'ystrip':False,'color_specific':True,'goodsell':True})
			for resname in lipid_colors_full.keys():
				view.set_color_cursor(lipid_colors_full[resname])
				view.select(**{'residue_lipid_%s'%resname:'noh and resname %s'%resname,
					'smooth':False,'style':licorice_thick,standard_view:True,
					'color_specific':True})
				view.select(**{'residue_lipid_%s_glass'%resname:'noh and resname %s'%resname,
					'style':licorice_thick,'smooth':False,alt_view:True,
					'color_specific':True,'ystrip':True})
			view.command("animate goto %d"%frame)
			view.command('package require pbctools')
			view.command('pbc box -color black')
			view.do('reset','xview')
			view.command('scale by %.2f'%zoom)
			view.do('snapshot')
			view.show(quit=True)
			print('done %d/%d'%(aa,len(alt_views)))
