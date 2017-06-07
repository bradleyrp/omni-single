#!/usr/bin/python

import time
from numpy import *
import MDAnalysis
from base.tools import status

def protein_rmsd(grofile,trajfile,**kwargs):

	"""
	Compute the RMSD of a protein.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	
	#---prepare universe	
	slice_name = kwargs['slice_name']
	group = kwargs['group']
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	protein = uni.select_atoms('protein and name CA')

	#---reference frame
	uni.trajectory[0]
	r0 = protein.positions
	r0 -= mean(r0,axis=0)

	#---collect coordinates
	nframes = len(uni.trajectory)
	coords,times = [],[]
	for fr in range(0,nframes):
		uni.trajectory[fr]
		r1 = protein.positions
		coords.append(r1)
		times.append(uni.trajectory.time)

	#---simple RMSD code
	rmsds = []
	for fr in range(nframes):
		status('RMSD',i=fr,looplen=nframes)
		r1 = coords[fr]
		r1 -= mean(r1,axis=0)
		#---computation of RMSD validated against VMD but no reflection
		U,s,Vt = linalg.svd(dot(r0.T,r1))
		signer = identity(3)
		signer[2,2] = sign(linalg.det(dot(Vt.T,U)))
		RM = dot(dot(U,signer),Vt)
		rmsds.append(sqrt(mean(sum((r0.T-dot(RM,r1.T))**2,axis=0))))

	#---pack
	attrs,result = {},{}
	result['rmsds'] = array(rmsds)
	result['timeseries'] = array(times)
	return result,attrs	
