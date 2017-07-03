#!/usr/bin/env python

import collections,shutil,tempfile
import scipy.ndimage
import scipy.spatial
import joblib

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from config import bash
from render.wavevids import print_birdseye_snapshot_render

sns = work.sns()
data,calc = plotload('undulations',work)
data_protein,calc = plotload('protein_abstractor',work)

#---prepare surfaces
surfs = [data[sn]['data']['mesh'].mean(axis=0)-data[sn]['data']['mesh'].mean(axis=0).mean() for sn in sns]

#---settings
print_review = False
ppn = 4
frameskipper = 1
cmap_name = 'seismic'
cmap_names = {'seismic':'seismic','redblue':'RdBu_r'}
cmap_mpl_name = cmap_names[cmap_name]
fs = {'axlabel':14,'title':20,'legend':14,'ticks':14,'tiny':10,'text_note':12,}
fs = dict([(i,j-4) for i,j in fs.items()])
nrows,ncols = 2,2
figsize = 8,8
panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[nrows,ncols],
	'hspace':0.4,'wspace':0.4}]},figsize=figsize)
handle = 'wavevid-'+'.'.join(sns)+'.%s'%cmap_name

figspot = work.paths['post_plot_spot']
outdir = figspot+'/vids/'
if not os.path.isdir(outdir): os.mkdir(outdir)
framecounts = [data[sn]['data']['nframes'] for sn in data]
nframes = min(framecounts)
if not all([framecounts[0]==framecounts[i] for i in range(len(framecounts))]): 
	status('[WARNING] uneven number of frames for the video: '+str(framecounts))
	status('[NOTE] defaulting to lowest number of frames: '+str(nframes))
frameset = np.arange(0,nframes,frameskipper)
surfs = np.array([s[:nframes] for s in surfs])
nprots_list = [work.meta[sn].get('nprots',1) for sn in sns]
protein_pts_all = np.array([data_protein[sn]['data']['points_all'] for sn in sns])
#---reformulate these for clarity in the joblib call below
protein_pts = [[protein_pts_all[ii][fr][...,:2] for ii,sn in enumerate(sns)] for fr in range(nframes)]
mvecs = np.array([data[sn]['data']['vecs'][:nframes] for sn in sns])
titles = [work.meta[sn]['label'] for sn in sns]

#---parallel render
tmpdir = os.path.join(outdir,handle)
if os.path.isdir(tmpdir): raise Exception('refusing to rerender to %s'%tmpdir)
else: os.mkdir(tmpdir)
extrema = surfs.min(),surfs.max()
#---symmetric colors
extrema = [(-1*i,i) for i in [max(np.abs(extrema))]][0]
if print_review:
	fi,fr = 10,frameset[10]
	print_birdseye_snapshot_render(
		surfs[:,fr],protein_pts[fr],mvecs[:,fr],nprots_list,
		handle=handle,outdir=tmpdir,fi=fi,fn=handle+'.fr.%04d'%fi,
		pbc_expand=1.0,smooth=1.,panelspecs=panelspecs,
		extrema=extrema,fs=fs,titles=titles,cmap_mpl_name=cmap_mpl_name)
	sys.exit(1)
joblib.Parallel(n_jobs=ppn,verbose=10)(
	joblib.delayed(print_birdseye_snapshot_render,joblib.pool.has_shareable_memory)(
		surfs[:,fr],protein_pts[fr],mvecs[:,fr],nprots_list,
		handle=handle,outdir=tmpdir,fi=fi,fn=handle+'.fr.%04d'%fi,
		pbc_expand=1.0,smooth=1.,panelspecs=panelspecs,
		extrema=extrema,fs=fs,titles=titles,cmap_mpl_name=cmap_mpl_name,
		) for fi,fr in enumerate(frameset))

#---print the film
ffmpegcmd =  ['ffmpeg','-i',tmpdir+'/'+handle+'.fr.%04d.png','-vcodec','mpeg2video','-q:v','0']+\
	(['-filter:v','setpts=2.0*PTS'] if False else [])+[' '+outdir+'/'+handle+'.mpeg']	
print(' '.join(ffmpegcmd))
#---! ffmpeg call is currently broken
bash(ffmpegcmd,cwd=work.paths['post_plot_spot'])
status('[STATUS] video rendered to '+outdir+'/'+handle+'.mpeg')
