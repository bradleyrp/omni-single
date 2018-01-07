#!/usr/bin/env python

__doc__ = """
Endpoint calling dynamics_videos.py with a scene file to make videos.
Note that we never delete temporary directories with images.
All videos are saved to vmdmake_videos subdirectory of the plot directory.
From inside an IPython notebook, the user will have to run `make_plots()` to render
and then `view_plots()` to view embedded videos linked into the `calcs` folder from the plot directory.
"""

import json,glob,os

# live: %%javascript
# live: IPython.OutputArea.auto_scroll_threshold = 9999;

# block: functions
@autoload(plotrun)
def load():
	"""Prepare the videomaker."""
	global dv
	from codes import dynamics_videos as dv
	if is_live: 
		status('run `make_videos()` then `view_videos()` to display videos in the notebook')
		return
	# export the plot environment to dv which resembles a plot script
	for key in _plot_environment_keys+['plotname','work']: dv.__dict__[key] = globals()[key]
	dv.prepare()

def view_videos():
	"""Display videos in an IPython Notebook."""
	if not is_live: return
	from IPython.display import HTML
	videodir = os.path.join(work.plotdir,'vmdmake_videos')
	if not os.path.isdir(videodir): raise Exception('cannot find %s'%videodir)
	vids = glob.glob(os.path.join(videodir,'*.webm'))
	for vid in vids:
		vid_local = os.path.basename(vid)
		try: os.system('ln -s %s'%vid)
		except: pass
		display(HTML(('<h3>%s</h2><video width="50%%" '
			'controls><source src="%s" type="video/mp4"></video>')%(vid_local,vid_local)))

@autoplot(plotrun)
def make_videos():
	"""Make videos."""
	# live notebook must call the CLI script
	if is_live: 
		# live shell: cd .. && make plot video_maker
		return
	# read the scene instructions
	scene_fn = work.plotspec.specs.get('specs',{}).get('scene',None)
	# assume the scene file is in the calcs folder
	scene_fn = os.path.join('../calcs/' if is_live else 'calcs',scene_fn)
	if not os.path.isfile(scene_fn): raise Exception('cannot find %s'%scene_fn)
	scene = {}
	with open(scene_fn) as fp: exec(fp.read(),globals(),scene)
	if 'interpreter' not in scene: raise Exception('cannot find interpreter function in %s'%scene_fn)
	instructions = scene['interpreter']()
	# write the instructions to JSON
	with open('video_requests.json','w') as fp: json.dump(instructions,fp)
	# render each video from the instructions
	dv.render_from_json()
