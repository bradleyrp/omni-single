#!/usr/bin/env python

"""
Non-interactive codes for the automatic videomaker adapted for protein demo.
Note that this script must be called from plot-protein_video_interactive.py for now.
"""

import json
from config import bash

@autoload(plotrun)
def load():
	"""Standard load sequence."""
	global vmdmake,data,sns,calc,tempdir
	plotname = 'protein_rmsd'
	try: import vmdmake
	except: 
		vmdmake_spot = os.path.join('calcs','vmdmake')
		if os.path.isdir(vmdmake_spot): raise Exception('could not import vmdmake but %s exists'%vmdmake_spot)
		else: bash('git clone http://github.com/bradleyrp/amx-vmd vmdmake',cwd='calcs')
	if 'data' not in globals(): data,calc = plotload(plotname)
	drop_dn = 'vmdmake_videos'
	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],drop_dn)
	if not os.path.isdir(tempdir): os.mkdir(tempdir)
	status('snapshots are dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')
	do_video,do_smooth = True,True
	sns = work.sns()

@autoplot(plotrun)
def render_from_json():
	"""
	Read video_catalog.json and make videos from it if they do not already exist.
	Video viles are written back to video_catalog.json.
	"""
	#---import settings due to problems with 
	with open('video_requests.json') as fp: film_cuts = json.loads(fp.read())
	vids = dict([(sn,{}) for sn in sns])
	for cut_name,cut_spec in film_cuts.items():
		#---loop over simulations
		for sn in sns:
			#---settings
			slice_path = calc['extras'][sn]['slice_path']
			gro,xtc = [os.path.join(work.postdir,'%s.%s'%(slice_path,j)) for j in ['gro','xtc']]
			tpr = work.raw.get_last(sn,subtype='tpr')
			view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
				frames='',xres=4000,yres=4000,**cut_spec.get('kwargs',{}))
			view.do('load_dynamic','standard',*cut_spec.get('does',[]))
			for sel in cut_spec.get('selections',[]): view.select(**sel)
			view.do('reset','xview')
			view.command('scale by 1.2')
			view.video()
			view.command('animate goto last')
			view.command('')
			view['snapshot_filename'] = 'snap.%s.%s'%(cut_name,sn)
			view.do('snapshot')
			#---no repeats 
			#---! currently the video suffix is hard-coded
			vid_fn = os.path.join('vid.%s.%s.%s'%(cut_name,sn,'mp4'))
			is_complete = os.path.isfile(os.path.join(tempdir,vid_fn))
			if not is_complete: 
				view.show(quit=True)
				view.render(name='vid.%s.%s'%(cut_name,sn))
			vid_fn_mp4 = os.path.join(tempdir,vid_fn)
			#---intervene to convert this to webm
			vid_fn = re.sub('\.mp4$','.webm',vid_fn_mp4)
			if not is_complete: bash('ffmpeg -i %s %s'%(vid_fn_mp4,vid_fn))
			#---save the video file and cut details
			#---! should you resave more cut information?
			vids[sn][cut_name] = vid_fn
	#---files go back out to the interactive python
	with open('video_catalog.json','w') as fp: json.dump(vids,fp)
