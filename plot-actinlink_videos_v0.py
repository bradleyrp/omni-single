#!/usr/bin/env python

"""
"""

#---block: imports
try: from codes import vmdmake
except: raise Exception(
	'cannot import vmdmake. try `git clone http://github.com/bradleyrp/amx-vmd %s`'%
	os.path.join(os.getcwd(),'calcs','codes','vmdmake'))
#---! hard-coded paths
kwargs_vmdmake = {'CGBONDSPATH':'~/libs/cg_bonds.tcl','GMXDUMP':'~/libs/gmxdump'}

#---block: load
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)
	sns = work.sns()

#---block: make videos
drop_dn = 'vmdmake_videos'
#---store the snapshots in the post_plot_spot according to the tag
tempdir = os.path.join(work.paths['post_plot_spot'],drop_dn)
if not os.path.isdir(tempdir): os.mkdir(tempdir)
status('snapshots dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')

#---block: render videos
do_video = True
do_smooth = True
sns_ordered = sns#['gelbilayerphys']
lipid_colors = {'PI2P':'cyan','DOPS':'pink','DOPE':'blue','POPC':'gray','CHL1':'green'}
film_cuts = {
	'protein_bilayer.binds':{'lipid_style':'phosphate_beads','protein_residue_highlight':['LYS','ARG']},
	'protein_bilayer.dynamics':{'lipid_style':'lines','show_near_ptdins':False},
	}

#---loop over video styles
for cut_name,cut_spec in film_cuts.items():
	#---loop over simulations
	for sn in sns_ordered:

		#---settings
		ptdins_cutoff = 10.0
		non_ptdins_atom_names = 'name P or name P3'
		ptdins_resname = work.meta[sn]['ptdins_resname']
		#---make a video of the nearby PIP2
		slice_path = calc['extras'][sn]['slice_path']
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(slice_path,j)) for j in ['gro','xtc']]
		tpr = work.raw.get_last(sn,subtype='tpr')
		view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
			frames='',xres=4000,yres=4000,**kwargs_vmdmake)
		view.do('load_dynamic','standard','bonder')
		view.select(protein='protein',style='cartoon',structure_color=True,smooth=do_smooth,goodsell=True)
		if cut_spec.get('protein_residue_highlight',None):
			view.select(lys='protein and (%s)'%' or '.join(['resname %s'%i for i in 
				cut_spec['protein_residue_highlight']]),
				style='licorice',smooth=do_smooth,goodsell=True,resname_color=True)
		if cut_spec.get('show_near_ptdins',True):
			select_ptdins = ('not hydrogen and (resname %s and '
				'same residue as (resname %s and within %.1f of protein))'%(
				ptdins_resname,ptdins_resname,ptdins_cutoff))
			view.select(near_pip2=select_ptdins,smooth=do_smooth,goodsell=True,style='licorice',update=True)
		#---loop over lipids
		for lipid_name,lipid_color in lipid_colors.items():
			view.set_color_cursor(lipid_color)
			if cut_spec['lipid_style']=='phosphate_beads':
				view.select(**{'lipid_%s'%lipid_name:
					'(%s) and not hydrogen and (%s) and not (%s)'%(
					non_ptdins_atom_names,'resname %s'%lipid_name,
					select_ptdins),'smooth':do_smooth,'style':'beads','color_specific':True,'goodsell':True})
			elif cut_spec['lipid_style']=='lines':
				view.select(**{'lipid_%s_lines'%lipid_name:'not hydrogen and (%s) and not (%s)'%(
					'resname %s'%lipid_name,select_ptdins),
					'style':'Licorice 0.400000 12.000000 12.000000',
					'goodsell':True,'smooth':do_smooth,'color_specific':True})
			else: raise Exception('unclear lipid style: %s'%cut_spec['lipid_style'])
		view.do('reset','xview')
		view.command('scale by 1.2')
		if do_video: view.video()
		view.command('animate goto last')
		view.command('')
		view['snapshot_filename'] = 'snap.%s.%s'%(cut_name,sn)
		view.do('snapshot')
		view.show(quit=True)
		if do_video: view.render(name='vid.%s.%s'%(cut_name,sn))
