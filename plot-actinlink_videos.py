#!/usr/bin/env python

__doc__ = """
Endpoint calling dynamics_videos.py with a scene file to make videos.
Note that we never delete temporary directories with images.
All videos are saved to vmdmake_videos subdirectory of the plot directory.
"""

import json

@autoload(plotrun)
def load():
	"""Prepare the videomaker."""
	global dv
	from codes import dynamics_videos as dv
	# export the plot environment to dv which resembles a plot script
	for key in _plot_environment_keys+['plotname','work']: dv.__dict__[key] = globals()[key]
	dv.prepare()

@autoplot(plotrun)
def plot():
	# read the scene instructions
	scene_fn = work.plotspec.specs.get('specs',{}).get('scene',None)
	# assume the scene file is in the calcs folder
	scene_fn = os.path.join('calcs',scene_fn)
	if not os.path.isfile(scene_fn): raise Exception('cannot find %s'%scene_fn)
	scene = {}
	with open(scene_fn) as fp: exec(fp.read(),globals(),scene)
	if 'interpreter' not in scene: raise Exception('cannot find interpreter function in %s'%scene_fn)
	instructions = scene['interpreter']()
	# write the instructions to JSON
	with open('video_requests.json','w') as fp: json.dump(instructions,fp)
	# render each video from the instructions
	dv.render_from_json()
