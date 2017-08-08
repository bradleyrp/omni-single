#!/usr/bin/env python

"""
Import regular membrane data for curvature coupling calculation.
Customized for the "dextran" project by "samrad".
"""

import os,re,sys
import numpy as np
from omnicalc import store,load
from base.tools import status

def repacker(source,source_out):
	"""
	Package regular membrane data to mimic the undulations calculation from omnicalc.
	"""
	data = {}
	#---read a dat file into the undulations data structure
	#---customized for data by samrad
	with open(source) as fp: text = fp.read()
	regex_frame = '(\d+)\n(.*?)(?=\n\d+\n|\Z)'
	frames = re.findall(regex_frame,text,flags=re.M+re.DOTALL)
	for frameno,frame in frames:
		data[int(frameno)] = np.array([i.split() for i in frame.splitlines()]).astype(float)

	mesh = []
	framenos = sorted(data.keys())
	xs_master,ys_master = [np.unique(data[0][:,j]) for j in [1,2]]
	#---reformulate the data into the mesh objects
	for frameno in framenos:
		frame = data[frameno]
		xs,ys = [np.unique(frame[:,j]) for j in [1,2]]
		if not np.all(xs==xs_master) or not np.all(ys==ys_master): 
			raise Exception('inconsistent grid points')
		regular = np.zeros((len(xs),len(ys)))
		inds = frame[:,1:2+1].astype(int)-1
		regular[(inds[:,0],inds[:,1])] = frame[:,3]
		if not np.all(regular!=0.0): raise Exception('some points are zero!')
		mesh.append(regular)

	result,attrs = {},{}
	#---the mesh object is per monolayer
	result['mesh'] = np.array([mesh])
	result['vecs'] = np.array([(xs_master.ptp()+1,ys_master.ptp()+1,0.0) for i in framenos])
	result['grid'] = result['vecs'][0].astype(int)
	result['nframes'] = np.array(len(framenos))
	result['framenos'] = np.array(framenos)
	attrs['grid_spacing'] = 1.0
	#---! check
	if False: plt.imshow(regular.T,origin='lower',interpolation='nearest');plt.show()
	#---save the new object
	store(obj=result,name=source_out,path='/home/localshare/factory/data/dextran/post',attrs=attrs)

#---list of source data
post_dn = '/home/localshare/factory/data/dextran/post'
sources = {
	'free_membrane_cl0':"/store-iota/MASTER/Samrad/Mem_Undulation/MEM_Position-0_Frame-0.dat",
	'bound_membrane_cl0':"/store-iota/MASTER/Samrad/Mem_Undulation/MEM_Position-0_Frame-1.dat",}

#---repackage all sources
for sn,source in sources.items():
	source_out = 'repacked_%s'%os.path.basename(source)
	drop_fn = os.path.join(post_dn,source_out)
	if not os.path.isfile(drop_fn): repacker(source,source_out)

def plotloader_for_dextran(calcname):
	"""
	Override for the `plotload` function in omnicalc to import data for curvature coupling from repackaged
	copies of regular membrane grids.
	"""
	calc,data = {'undulations':{}},{'undulations':{}}
	#---collect the repackaged data
	for sn,source in sources.items():
		'repacked_%s'%os.path.basename(source)
		data['undulations'][sn] = {'data':load(name=source_out,cwd=post_dn)}
	calc['calcs'] = {'undulations':{'specs':{'slice_name':'slice_name_placeholder'}}}
	return data,calc
