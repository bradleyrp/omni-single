#!/usr/bin/env

import os,json,copy,re
import numpy as np
# get the bash function from omni.config
from config import bash

def prepare():
	"""
	Run once to prepare for videos.
	"""
	global vmdmake,tempdir
	try: import vmdmake
	except: 
		vmdmake_spot = os.path.join('calcs','vmdmake')
		if os.path.isdir(vmdmake_spot): 
			raise Exception('could not import vmdmake but %s exists'%vmdmake_spot)
		else: bash('git clone http://github.com/bradleyrp/amx-vmd vmdmake',cwd='calcs')
	if 'data' not in globals(): 
		global data,calc
		data,calc = plotload(plotname)
	drop_dn = 'vmdmake_videos'
	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],drop_dn)
	if not os.path.isdir(tempdir): os.mkdir(tempdir)
	status('snapshots are dropping to temporary folders in %s. it is your responsibility to delete them.'%
		tempdir+' after the video is created we will not remake it unless you delete it.',tag='note')
	do_video,do_smooth = True,True

def render_from_json(request_fn='video_requests.json',catalog_fn='video_catalog.json'):
	"""
	Read video_catalog.json and make videos from it if they do not already exist.
	Video viles are written back to video_catalog.json.
	"""
	sns = work.sns()
	#---import settings due to problems with 
	if not os.path.isfile(request_fn): raise Exception('cannot find %s'%request_fn)
	with open(request_fn) as fp: film_cuts = json.loads(fp.read())
	vids = dict([(sn,{}) for sn in sns])
	for cut_name,cut_spec_loop in film_cuts.items():
		#---loop over simulations
		for sn in sns:
			#---make a copy because we modify some arguments
			cut_spec = copy.deepcopy(cut_spec_loop)
			#---settings
			slice_path = calc['extras'][sn]['slice_path']
			gro,xtc = [os.path.join(work.postdir,'%s.%s'%(slice_path,j)) for j in ['gro','xtc']]
			#! legacy plot mode. remove this once omnicalc development branch is completed.
			#! ... it is fairly inefficient to parse the entire source just to get a single TPR 
			#! ... so consider using meta to specify it explicitly, or some other method?
			try: tpr = work.raw.get_last(sn,subtype='tpr')
			except:
				work.parse_sources()
				tpr = work.source.get_last(sn,subtype='tpr')
			#! limit the number of frames if it is excessive. overridden by step in kwargs in cut_spec
			nframes_max = cut_spec.get('nframes_max',None)
			if nframes_max:
				nframes = ((calc['extras'][sn]['end']-calc['extras'][sn]['start'])/calc['extras'][sn]['skip'])
				# get the largest integer step size that will keep the number of frames below the max
				step = int(np.ceil(float(nframes)/nframes_max))
				if step<1: raise Exception('negative step size')
			else: step = 1
			view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
				frames='',res=cut_spec.pop('resolution',(1000,1000)),step=step,**cut_spec.get('kwargs',{}))
			view.do('load_dynamic','standard',*cut_spec.get('does',[]))
			for sel in cut_spec.get('selections',[]): 
				#! elaborate color handling here
				color_specific = sel.get('color_specific',None)
				color_specific_this = None
				if type(color_specific)==dict and set(color_specific.keys())=={'eval'}:
					color_specific_this = eval(color_specific['eval'])
					sel['color_specific'] = True
				elif type(color_specific) in str_types:
					color_specific_this = str(color_specific)
					sel['color_specific'] = True
				if color_specific_this!=None: 
					view.set_color_cursor(color_specific_this)
				view.select(**sel)
			view.do('reset','xview')
			view.command('scale by %s'%cut_spec.get('zoom',1.2))
			if not cut_spec.get('debug',False): view.video()
			view.command('animate goto last')
			# render in higher detail
			view.command('display resize %d %d'%tuple(cut_spec.pop('resolution_snapshot',(4000,4000))))
			view['snapshot_filename'] = 'snap.%s.%s'%(cut_name,sn)
			view.do('snapshot')
			#---no repeats 
			#---! currently the video suffix is hard-coded
			vid_fn = os.path.join('vid.%s.%s.%s'%(cut_name,sn,'mp4'))
			is_complete = os.path.isfile(os.path.join(tempdir,vid_fn))
			if not is_complete and not cut_spec.get('debug',False): 
				view.show(quit=True)
				view.render(name='vid.%s.%s'%(cut_name,sn))
			elif cut_spec.get('debug',False):
				status('debug mode runs view.show to make a snapshot')
				view.show(quit=True)
			else: status('video is already complete at %s'%os.path.join(tempdir,vid_fn))
			if not cut_spec.get('debug',False):
				vid_fn_mp4 = os.path.join(tempdir,vid_fn)
				#---intervene to convert this to webm
				vid_fn = re.sub('\.mp4$','.webm',vid_fn_mp4)
				if not is_complete: bash('ffmpeg -y -i %s -b:v 0 -crf 20 %s'%(vid_fn_mp4,vid_fn))
				#---save the video file and cut details
				#---! should you resave more cut information?
				vids[sn][cut_name] = vid_fn
	#---files go back out to the interactive python
	with open(catalog_fn,'w') as fp: json.dump(vids,fp)
