#!/usr/bin/env python

"""
Rendering instructions for a protein video.
"""

import os

drop_dn = 'vmdmake_videos'
do_smooth = True
lipid_material = 'diffuse'
film_cuts = {
	'bilayer.side':{'debug':False,'zoom':1.8,'does':['bonder'],'nframes_max':500,
		'kwargs':{'cgbondspath':os.path.join(os.getcwd(),'amx-martini/bin/cg_bonds.tcl')},
		'selections':[dict(subject='protein and noh',style='Licorice 0.6 12.0 12.0',
				smooth=do_smooth,goodsell=True),
			dict(subject_cartoon='protein and noh',style='cartoon',diffuse=True,
				smooth=do_smooth,goodsell=True),]},}

def interpreter(): 
	# return film cuts for dynamics_videos.py
	global film_cuts
	return film_cuts
