#!/usr/bin/python

import time
from numpy import *
from base.tools import status

def msd_fit_deprecated(time,displacement,leftcut=0.1,rightcut=0.9,dims=3):

	"""
	Fit mean-squared displacement curves to determine the diffusion coefficient.
	"""
	
	if dims==3: factor = 6.
	elif dims==2: factor = 4.
	else: raise Exception('[ERROR] dimensions must be 2 or 3 for diffusion')
	xvals,yvals = time,displacement
	datbuf = [int(len(yvals)*leftcut)+2,int(len(yvals)*rightcut)+1]
	[q,r], covariance = polyfit((xvals[datbuf[0]:datbuf[1]]),(yvals[datbuf[0]:datbuf[1]]),1,cov=True)
	return q/factor,r/factor
	
def msd_fit(time,displacement,cuts=None,dims=3,factor=None):

	"""
	Fit mean-squared displacement curves to determine the diffusion coefficient.
	Uses cutofss in the same units as time.
	"""
	
	if factor==None:
		if dims==3: factor = 6.
		elif dims==2: factor = 4.
	xvals,yvals = time,displacement
	#---! added log recently
	ins = where(all([xvals>=log(cuts[0]),xvals<=log(cuts[1])],axis=0))[0]
	try: q,r = polyfit(xvals[ins],yvals[ins],1)
	except: 
		print "ERROR"
		import pdb;pdb.set_trace()
	return q/factor,r/factor

def diffusion(coords,times,vecs,skip=1):

	"""
	Generic diffusion calculator for simulations with periodic boundary conditions.
	"""

	#---produce unbroken trajectory
	status('diffusion making unbroken trajectory',tag='compute')
	ncoords = shape(coords)[1]
	for d in range(3):
		vectile = tile(vecs[1:,d],(ncoords,1)).T
		coords[1:,:,d] -= (cumsum((coords[1:,:,d]-coords[:-1,:,d])>(vectile/2.),axis=0)-\
			cumsum((coords[:-1,:,d]-coords[1:,:,d])>(vectile/2.),axis=0))*vectile

	#---save delta-t and distances for all pairwise indices in the times list
	if type(times)==list: times = array(times)
	ntimes = len(times)
	dists,dts,msdraw = [],[],[]
	timechunks = [array([(i,i+d) for i in range(0,ntimes-d,skip)]) for d in range(1,ntimes)]
	start = time.time()
	for jj,j in enumerate(timechunks):
		status('diffusion distances',i=jj,looplen=len(timechunks),tag='compute',start=start)
		dists.append(linalg.norm(coords[j[:,1]]-coords[j[:,0]],axis=2))
		dts.append(times[j[:,1]]-times[j[:,0]])
		if 0: msdraw.append(dists[i]**2/tile(dts[i],(ncoords,1)).T)

	#---sort to account for jitter
	dtsu = unique(concatenate([unique(i) for i in dts]))
	msdsort = [[] for i in dtsu]
	start = time.time()
	for jj,j in enumerate(timechunks):
		status('diffusion, sorting',i=jj,looplen=len(timechunks),tag='compute',start=start)
		for kk,k in enumerate(dts[jj]):
			msdsort[where(dtsu==k)[0][0]].append(dists[jj][kk])
	
	#---average over observations for each timestep
	msdavg = zeros((len(dtsu),ncoords))
	start = time.time()
	for jj,j in enumerate(dtsu):
		status('diffusion, averaging',i=jj,looplen=len(dtsu),tag='compute',start=start)
		msdavg[jj] = mean(msdsort[jj],axis=0)

	attrs,result = {},{}
	result['msd'] = msdavg
	result['dt'] = dtsu
	del dists,dts,msdsort
	return result,attrs
	
def diffusion_zoned(coords,times,vecs,timechunks,ion_indices,skip=1):

	"""
	Generic diffusion calculator for simulations with periodic boundary conditions.
	"""

	#---produce unbroken trajectory
	status('diffusion, making unbroken trajectory',tag='compute')
	ncoords = shape(coords)[1]
	for d in range(3):
		vectile = tile(vecs[1:,d],(ncoords,1)).T
		coords[1:,:,d] -= (cumsum((coords[1:,:,d]-coords[:-1,:,d])>(vectile/2.),axis=0)-\
			cumsum((coords[:-1,:,d]-coords[1:,:,d])>(vectile/2.),axis=0))*vectile

	#---save delta-t and distances for all pairwise indices in the times list
	if type(times)==list: times = array(times)
	ntimes = len(times)
	dists,dts,msdraw = [],[],[]

	start = time.time()
	#---we assemble all motions in the zone together regardless of ion identity
	dists,dts = [],[]
	#---loop over ions
	for ii,ionnum in enumerate(ion_indices):
		status('diffusion, distances',i=ii,looplen=len(ion_indices),tag='compute',start=start)
		#---assemble the timestep pairs for this ion in the zone
		pairs = timechunks[ii]
		#---some ions never cross zones and hence are excluded from the calculation
		if len(pairs)!=0:
			#---slice the coordinates to return only this ion
			subcoords = coords[:,ionnum]
			these_dists = linalg.norm(subcoords[pairs[:,0]]-subcoords[pairs[:,1]],axis=1)
			these_dts = pairs[:,1]-pairs[:,0]
			dists.extend(these_dists)
			dts.extend(these_dts)
	dists,dts = array(dists),array(dts)
			
	#---sort to account for jitter
	dtsu = unique(dts)
	msdsort = [[] for i in dtsu]
	start = time.time()
	for jj,j in enumerate(dtsu):
		status('diffusion, sorting',i=jj,looplen=len(dtsu),tag='compute',start=start)
		msdsort[jj] = dists[where(dts==j)[0]]

	#---average over observations for each timestep
	msdavg = zeros((len(dtsu)))
	start = time.time()
	for jj,j in enumerate(dtsu):
		status('diffusion, averaging',i=jj,looplen=len(dtsu),tag='compute',start=start)
		msdavg[jj] = mean(msdsort[jj],axis=0)

	attrs,result = {},{}
	result['msd'] = msdavg
	result['dt'] = dtsu
	del dists,dts,msdsort
	return result,attrs

