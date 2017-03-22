#!/usr/bin/env python

"""
CURRENTLY PORTING THIS! WORK IN PROGRESS!

#---plot prep
from plotter import *
from base.store import plotload
import numpy as np
from base.store import plotload,load
sys.path.insert(0,'calcs')
from codes.vmdwrap import *
from base.gromacs import gmxpaths
import collections
import shutil
"""

#---settings
film_folder = 'films'
slice_name = 'medium2'
group = 'protein'
collection = 'all'

#---prepare directories and get simulation names
out_dn = work.plotdir+'/%s/'%film_folder
if not os.path.isdir(out_dn): os.mkdir(out_dn)
sns = work.vars['collections'][collection]

#---one video per simulation
for sn in sns:

	#---! add this to workspace
	work.cursor = (work.c,'xtc')
	gro,traj = [os.path.join(work.postdir,work.slice(sn)[slice_name][group][suf]) 
		for suf in ['gro',work.trajectory_format]]
	tpr = work.get_last_start_structure(sn,partname='tpr')
	#---! hacked
	outname = 'vid.%s'%re.findall('^.+\/(.+)\.xtc$',traj)[0]
	snapshot_outname = re.findall('^.+\/(.+)\.xtc$',traj)[0]
	#---render
	v = VMDWrap(output=out_dn,subdir=sn,
		gro=gro,xtc=traj,tpr=tpr,res=(800,800),last=None,backcolor='white')
	v.do(*'pbc load standard bonder'.split())
	v.select(protein='protein',style='NewCartoon',
		structure_color=True,smooth=True,xy=False,goodsell=True)
	v.do('reset','yview')
	v.do('align_backbone','ortho')
	v.video()
	v.show(render=outname,quit=True,clean=False,prompt=False,text=True)
	subprocess.call("ffmpeg -y -i %s.mp4 -pix_fmt rgb24 %s.gif"%(outname,outname),
		cwd=os.path.join(out_dn),shell=True)
	subprocess.call(
		'ffmpeg -y -ss 00:00:00 -i %s/%s.mp4 -frames:v 1 %s/fig.snapshot.%s.png'%
		(out_dn,outname,work.plotdir,snapshot_outname),
		shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=work.plotdir)
	shutil.rmtree(os.path.join(out_dn,sn))
