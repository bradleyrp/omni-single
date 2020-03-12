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


lipid_colors = {'PI2P':'purple','SAPI':'purple',
	'P35P':'purple','CHL1':'gray','DOPS':'mauve','DOPE':'iceblue','POPC':'white','DOPC':'white'}
for key in ['ptdins','P35P','PIPP','PIPU','SAPI','PI2P']: 
	lipid_colors[key] = 'magenta'
colors = lipid_colors
cation_color = 'green'
for key in ['Cal','Mg']: colors[key] = cation_color
anion_color = 'gray'
for key in ['CL']: colors[key] = anion_color
all_lipids = ' or '.join(['resname %s'%r for r in lipid_colors.keys()])

#! should not be necessary!
sn = 'membrane-v532'
#! want to recode? # frate = 24.0 # if not size: bitrate = None  
#! elif duration==0.0: bitrate = float(size)*8192/(nframes/frate)


drop_dn = 'vmdmake_videos_v1'
do_smooth = True
lipid_material = ['goodsell','glass1','edgyglass','diffuse'][0]
film_cuts = {
	'bilayer.side':{'debug':False,'zoom':2.6,'does':['bonder'],
		'nframes_max':300,'resolution':(3000,2000),'video_size':25.0,
		'kwargs':{'cgbondspath':os.path.join(os.getcwd(),'amx/bin/cg_bonds.tcl')},
		'selections':[
			{'lipids_r%d'%rnum:'noh and resname %s'%resname,
				'style':'Licorice 1.3 12.0 12.0','smooth':do_smooth,lipid_material:True,
					'color_specific':colors[resname]}
			for rnum,resname in enumerate(work.vars['selectors']['resnames_lipid']+['CHL1'])]+[
			dict(anion='name %s and within 20 of (%s)'%(work.meta[sn]['anion'],all_lipids),
				style='vdw',smooth=do_smooth,goodsell=True,
				color_specific=colors[work.meta[sn]['anion']]),
			dict(cation='name %s and within 20 of (%s)'%(work.meta[sn]['cation'],all_lipids),
				style='vdw',smooth=do_smooth,goodsell=True,
				color_specific=colors[work.meta[sn]['cation']]),
			]},}

def interpreter(): 
	# return film cuts for dynamics_videos.py
	global film_cuts
	return film_cuts

#! iterated on nframes_max 5 until I liked what I saw
