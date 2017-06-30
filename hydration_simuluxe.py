#!/usr/bin/python

if 0: compsign,project_name = 'hydration','project-ptdins'
if 0: execfile('../../header.py')
from smx import *
from smxlocal import *
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
from numpy import *
import scipy
import scipy.spatial

def compute_hydration_shell(ipts,wpts,ld,vec=None,cutoff_hydration=2.3,cutoff=4.6,pbc=True,left_cutoff=None):

	"""
	"""

	if left_cutoff == None: inside = ld<cutoff
	else: inside = all((ld<cutoff,ld>=left_cutoff),axis=0)
	isbound = unique(transpose(where(inside))[:,0])
	ipts_subsel = ipts[isbound]
	if not pbc: 
		pd = scipy.spatial.distance.cdist(ipts_subsel,wpts)
		return sum(pd<cutoff_hydration,axis=1)
	else: 
		pd = array([scipy.spatial.distance.cdist(ipts_subsel[:,d:d+1],wpts[:,d:d+1]) for d in range(3)])
		for d in range(3):
			pd[d] -= (pd[d]>vec[d]/2.)*vec[d]
			pd[d] += (pd[d]<-1*vec[d]/2.)*vec[d]
		return sum(sqrt(sum(pd**2,axis=0))<cutoff_hydration,axis=1)

def compute_hydration(simname,metadat,grofile=None,trajfile=None,nrank=3,focus=None,panel=None,**kwargs):

	"""
	"""

	sn = simname
	
	headerdat = kwargs['headerdat']
	hydration_scan = headerdat['calculations']['hydration']['hydration_scan']
	#---hydration cutoff depends on the ion
	#---! note previous sn/simname error may have compromised those pickles
	cutoff_hydration = {'NA':3.05,'MG':2.3,'Cal':2.6,
		}[metadat[simname]['ion_name_positive']]

	mset = SimSetMembrane()
	print grofile,trajfile
	mset.load_trajectory((grofile,trajfile))

	solsel = mset.universe.selectAtoms('name OW')
	ionsel = mset.universe.selectAtoms('name '+metadat[simname]['ion_name_positive'])
	
	dat = loader(focus,headerdat,sn=sn,compsign='binding')[sn]
	ld = dat['lipid_distances']
	nframes = dat['nframes']

	st = time.time()
	icoords,wcoords = [],[]
	for fr in range(nframes):
		status('[LOAD] coordinates',i=fr,looplen=nframes,start=st)
		mset.gotoframe(fr)
		mset.vec(fr)
		icoords.append(ionsel.coordinates())
		wcoords.append(solsel.coordinates())
		
	fr = 100
	lcutoff = 2.2
	rcutoff = 3.0
	compute_hydration_shell(icoords[fr],wcoords[fr],ld[fr],
		vec=mset.vecs[fr],cutoff_hydration=cutoff_hydration,cutoff=rcutoff,left_cutoff=lcutoff)

	hydratedat = {}
	for ci in range(len(hydration_scan)-1):
		lcutoff,rcutoff = hydration_scan[ci:ci+2]
		status('[COMPUTE] hydration for cutoff='+str(rcutoff))
		outdat = Parallel(n_jobs=4,verbose=10)(
			delayed(compute_hydration_shell,has_shareable_memory)
			(icoords[fr],wcoords[fr],ld[fr],
			vec=mset.vecs[fr],cutoff_hydration=cutoff_hydration,cutoff=rcutoff,left_cutoff=lcutoff)
			for fr in range(nframes))
		for fr in range(nframes): hydratedat[str(lcutoff)+'-'+str(rcutoff)+':'+str(fr)] = array(outdat[fr])
		result = hydratedat
	attrs = {'nframes':nframes}
	return result,attrs
	
