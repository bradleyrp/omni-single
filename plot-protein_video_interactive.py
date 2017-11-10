#!/usr/bin/env python

__doc__ = """
Automatic videomaker adapted for protein demo.
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
#---load the data for protein_rmsd to get trajectory files
if 'data' not in globals(): data,calc = plotload('protein_rmsd')

#---block: video requests
drop_dn = 'vmdmake_videos'
do_smooth = True
film_cuts = {
	'protein.centered':{'does':['align_backbone'],'kwargs':{'backbonename':'name CA'},
		'selections':[dict(protein='protein',style='cartoon',structure_color=True,
			smooth=do_smooth,goodsell=True)]},
	'protein.basic_residues':{'selections':[
		dict(protein='protein',style='cartoon',structure_color=True,
			smooth=do_smooth,goodsell=True),
		dict(basic='resname ARG or resname LYS',style='licorice',
			smooth=do_smooth,goodsell=True),]},}

#---block: make videos
#---store the snapshots in the post_plot_spot according to the tag
tempdir = os.path.join(work.paths['post_plot_spot'],drop_dn)
if not os.path.isdir(tempdir): os.mkdir(tempdir)
status('snapshots are dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')
sns = work.sns()
#---! Jupyter has a bug currently forbidding Popen so we have to do this all day
with open('../video_requests.json','w') as fp: json.dump(film_cuts,fp)

#---block: make the videos directly in bash
!cd .. && make plot protein_video render_from_json

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
