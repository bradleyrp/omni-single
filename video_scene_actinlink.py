#!/usr/bin/env python

"""
The "interpreter" function returns instructions for making dynamics videos.
This scene must be called via e.g. `exec(this_script_text,globals(),local_dict)` where globals contains
the workspace variable (work) since it may use parts of the workspace to decide what to render. 
The function which sends these instructions to the vmdmake interface at codes/dynamics_videos.py should
run the interpreter from the local_dict populated after the exec function. If you think this scheme is too
complicated, ask Ryan! Note that vmdmake automates VMD and dynamics_videos.py connects it to omnicalc.
"""

import os

drop_dn = 'vmdmake_videos'
do_smooth = True
lipid_material = ['goodsell','glass1','edgyglass','diffuse'][-1]
film_cuts = {
	'bilayer.side':{'debug':False,'zoom':1.8,'does':['bonder'],'nframes_max':300,
		'kwargs':{'cgbondspath':os.path.join(os.getcwd(),'amx-martini/bin/cg_bonds.tcl')},
		'selections':[
			{'lipids_r%d'%rnum:'noh and resname %s'%resname,
				'style':'Licorice 0.3 12.0 12.0','smooth':do_smooth,lipid_material:True,
					'color_specific':{'eval':
					'colorize(work.meta[sn],resname=\"'+resname+
						'\",named=True,overrides={"CHL1":"white"})'}}
			for rnum,resname in enumerate(work.vars['selectors']['resnames_lipid']+['CHL1'])]+[
			dict(subject='protein and noh',style='Licorice 0.6 12.0 12.0',
				smooth=do_smooth,goodsell=True),
			dict(subject_cartoon='protein and noh',style='cartoon',diffuse=True,
				smooth=do_smooth,goodsell=True,color_specific='black'),
			]},}

def interpreter(): 
	# return film cuts for dynamics_videos.py
	global film_cuts
	return film_cuts
