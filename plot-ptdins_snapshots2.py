#!/usr/bin/env python

routine = ['full_bilayer_snapshot']

from base.tools import dictsum
try: from codes import vmdmake
except: raise Exception(
	'cannot import vmdmake. try `git clone http://github.com/bradleyrp/amx-vmd %s`'%
	os.path.join(os.getcwd(),'calcs','codes','vmdmake'))

if 'data' not in globals():
	data,calc = plotload(plotname)

if 'full_bilayer_snapshot' in routine:

	frame = 56
	sn = 'membrane-v532'
	vmd_fn = 'fig.snapshot.vmd.bilayer'
	do_ptdins_color = True
	lipid_style = {'goodsell':True,'style':'Licorice 1.3 12.0 12.0'}
	lipid_colors = {'PI2P':'purple','SAPI':'purple',
		'P35P':'purple','CHL1':'gray','DOPS':'red','DOPE':'blue2','POPC':'white','DOPC':'white'}
	#---trying different colors
	lipid_colors = {'PI2P':'purple','SAPI':'purple',
		'P35P':'purple','CHL1':'gray','DOPS':'mauve','DOPE':'iceblue','POPC':'white','DOPC':'white'}
	for key in ['ptdins','P35P','PIPP','PIPU','SAPI','PI2P']: lipid_colors[key] = 'magenta'
	cation_color = 'green'
	anion_color = 'gray'

	gro,xtc = [os.path.join(work.postdir,'%s.%s'%(calc['extras'][sn]['slice_path'],suf))
		for suf in ['gro','xtc']]
	#---get the tpr from the raw data
	tpr = work.raw.get_last(sn,subtype='tpr')
	view = vmdmake.VMDWrap(site=work.paths['post_plot_spot'],gro=gro,xtc=xtc,tpr=tpr,
		frames='',res=(4000,4000))
	view['snapshot_filename'] = vmd_fn
	view.do('load_dynamic','standard','bonder')
	view.command("animate goto %d"%(frame))
	resname = work.meta[sn]['ptdins_resname']
	select = {'ptdins':'resname %s and not hydrogen'%resname}
	if do_ptdins_color: 
		view.set_color_cursor(lipid_colors[resname])
		view.select(**dictsum(select,lipid_style,{'color_specific':True}))
	else: view.select(**dictsum(select,lipid_style))
	#---show other lipids in their respective colors
	for resname in [i for i in lipid_colors if i not in work.vars['selectors']['resnames_PIP2']]:
		select = {'res_%s'%resname:'resname %s and not hydrogen'%resname}
		view.set_color_cursor(lipid_colors[resname])
		view.select(**dictsum(select,lipid_style,{'color_specific':True}))
	#---show the ions
	all_lipids = ' or '.join(['resname %s'%r for r in lipid_colors.keys()])
	view.set_color_cursor(cation_color)
	select = {'cation':'name %s and within 20 of %s'%(work.meta[sn]['cation'],all_lipids)}
	view.select(**dictsum(select,{'goodsell':True,'style':'vdw'},{'color_specific':True}))
	view.set_color_cursor(anion_color)
	select = {'anion':'name %s and within 20 of %s'%(work.meta[sn]['anion'],all_lipids)}
	view.select(**dictsum(select,{'goodsell':True,'style':'vdw'},{'color_specific':True}))
	view.do('reset','xview')
	view.command('scale by 1.5')
	view.do('snapshot')
	view.show(quit=True)

	#---trim the image
	fig = plt.figure(figsize=(16,16))
	ax = fig.add_subplot(111)
	image = sidestack_images([os.path.join(work.plotdir,vmd_fn+'.png')])
	ax.imshow(image)
	ax.axis('off')
	picturesave('fig.snapshot.bilayer',work.plotdir,
		backup=False,version=True,meta={},extras=[])
