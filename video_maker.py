#!/usr/bin/env python

__doc__ = """
Endpoint calling dynamics_videos.py with a scene file to make videos.
Note that we never delete temporary directories with images.
All videos are saved to vmdmake_videos subdirectory of the plot directory.
From inside an IPython notebook, the user will have to run `make_plots()` to render
and then `view_plots()` to view embedded videos linked into the `calcs` folder from the plot directory.
"""

import json,glob,os

# block: functions
@autoload(plotrun)
def load():
	"""Prepare the videomaker."""
	global dv
	from codes import dynamics_videos as dv
	if is_live: 
		status('INSTRUCTIONS: run `make_videos()` then `view_videos()` to display videos in the notebook. '
			'videos are rendered to %s once only'%os.path.join(work.plotdir,'vmdmake_videos')+
			'snapshots are written to temporary folders in the vmdmake_video directory. '
			'users must clean up the snapshots or delete the videos to remake them. ')
		return
	# export the plot environment to dv which resembles a plot script
	for key in _plot_environment_keys+['plotname','work']: dv.__dict__[key] = globals()[key]
	dv.prepare()

@autoplot(plotrun)
def view_videos():
	"""Display videos in an IPython Notebook."""
	if not is_live: return
	from IPython.display import HTML
	videodir = os.path.join(work.plotdir,'vmdmake_videos')
	if not os.path.isdir(videodir): raise Exception('cannot find %s'%videodir)
	vids = glob.glob(os.path.join(videodir,'*.webm'))
	from IPython.core.magics.display import Javascript
	Javascript('IPython.OutputArea.auto_scroll_threshold = 9999;')
	for vid in vids:
		vid_local = os.path.basename(vid)
		try: os.system('ln -s %s'%vid)
		except: pass
		display(HTML(('<h3>%s</h2><video width="50%%" '
			'controls><source src="%s" type="video/mp4"></video>')%(vid_local,vid_local)))

@autoplot(plotrun)
def make_videos():
	"""Make videos."""
	# read the scene instructions
	scene_fn = work.plotspec.specs.get('specs',{}).get('scene',None)
	# assume the scene file is in the calcs folder
	scene_fn = os.path.abspath(os.path.join('../calcs/' if is_live else 'calcs',scene_fn))
	if not os.path.isfile(scene_fn): raise Exception('cannot find %s'%scene_fn)
	scene = {}
	with open(scene_fn) as fp: exec(fp.read(),globals(),scene)
	if 'interpreter' not in scene: raise Exception('cannot find interpreter function in %s'%scene_fn)
	instructions = scene['interpreter']()
	# live notebook must call the CLI script (ask Ryan why)
	if is_live: 
		status('scene instructions from %s are as follows: %s'%(scene_fn,instructions))
		# live shell: cd .. && make plot video_maker make_videos
		return
	# write the instructions to JSON
	with open('video_requests.json','w') as fp: json.dump(instructions,fp)
	# render each video from the instructions
	dv.render_from_json()

# block: interactive video rendering
# live: plotrun.loader()
# live: make_videos()
# block: interactive videos
# live: view_videos()