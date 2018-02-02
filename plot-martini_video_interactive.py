#!/usr/bin/env python

__doc__ = """
Automatic videomaker adapted for protein demo.
This script and its companion were adapted from the protein version.
"""

#---block: imports
if not is_live or not os.path.basename(os.getcwd())=='calcs':
	raise Exception('you must run this from jupyter in the factory')
import json
from config import bash
try: import vmdmake
except: 
	#---clone vmdmake codes if they are absent
	vmdmake_spot = os.path.join('vmdmake')
	if os.path.isdir(vmdmake_spot): raise Exception('could not import vmdmake but %s exists'%vmdmake_spot)
	else: bash('git clone http://github.com/bradleyrp/amx-vmd vmdmake')
#---clone the martini repo for the bonder code
if not os.path.isdir('amx-martini'): 
    bash('git clone http://github.com/bradleyrp/amx-martini')
#---use a calculation to get the trajectory files, set by the martini_video_interactive entry in plots
if 'data' not in globals(): data,calc = plotload(plotname)

#---block: video requests for MARTINI bilayer simulations
drop_dn = 'vmdmake_videos'
do_smooth = True
film_cuts = {
	'bilayer.side':{'does':['bonder'],'kwargs':{'cgbondspath':
		os.path.join(os.getcwd(),'amx-martini/bin/cg_bonds.tcl')},
		'selections':[dict(
			lipids=' or '.join(['resname %s'%i for i in work.vars['selectors']['resnames_lipid']]),
			style='Licorice 1.000000 12.000000 12.000000',smooth=do_smooth,goodsell=True)]},}

#---block: make videos
#---store the snapshots in the post_plot_spot according to the tag
tempdir = os.path.join(work.paths['post_plot_spot'],drop_dn)
if not os.path.isdir(tempdir): os.mkdir(tempdir)
status('snapshots are dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')
sns = work.sns()
#---! Jupyter has a bug currently forbidding Popen so we have to do this all day
with open('../video_requests.json','w') as fp: json.dump(film_cuts,fp)

#---block: make the videos directly in bash
!cd .. && make plot martini_video render_from_json

#---block: show the videos
if is_live:
	with open('../video_catalog.json') as fp: vids = json.loads(fp.read())
	from IPython.display import HTML
	#---loop over generated videos
	for sn in vids:
		for cut_key in vids[sn]:
			#---link in the video
			os.system('ln -s %s'%vids[sn][cut_key])
			#---! present a link for download here? webm might not be the preferred format
			#---display
			display(HTML(('<h3>%s, %s</h2><video width="50%%" '
				'controls><source src="%s" type="video/mp4"></video>')%(
				sn,cut_key,os.path.basename(vids[sn][cut_key]))))
