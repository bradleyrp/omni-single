#!/usr/bin/python

from joblib import Parallel,delayed
# from joblib.pool import has_shareable_memory
from numpy import *
import scipy
import scipy.spatial
from scipy.optimize import curve_fit
from scipy.misc import factorial

def ontorus(pts1,pts2,vecs):

	"""
	Compute distances under periodic boundary conditions with the topology of a torus.
	"""

	cdx = scipy.spatial.distance.cdist(pts1[...,0:1],pts2[...,0:1])
	cdy = scipy.spatial.distance.cdist(pts1[...,1:2],pts2[...,1:2])
	cdz = scipy.spatial.distance.cdist(pts1[...,2:],pts2[...,2:])
	cdx = cdx-1*(cdx>vecs[0]/2.)*vecs[0]
	cdy = cdy-1*(cdy>vecs[1]/2.)*vecs[1]
	cdz = cdz-1*(cdz>vecs[2]/2.)*vecs[2]
	cd = sqrt(cdx**2+cdy**2+cdz**2)
	return cd
	
def poisson(k, lamb): 

	"""Poisson distribution."""

	return (lamb**k/factorial(k))*np.exp(-lamb)

def partnerfinder(pts1,pts2,vec,lipid_resids,nrank,includes=slice(None,None)):

	"""
	Given lipid and ion points this function returns a description of nearest binding partners and their
	identities. Note that it refers to lipid_pts_collected and ion_pts_collected, which are globals.
	"""
		
	nions = len(pts2)
	#---pairwise distances between lipid atoms and ions
	cd = ontorus(pts1,pts2,vec).T
	#---argsorted matrix of nearest atom-to-ion distances
	aa = argsort(cd,axis=1)
	#---lipid_resids maps atom indices onto the lipid residue id
	#---lipid_resids[aa].T is a map of the closest lipid residues
	#---the sort-unique command collects the first 3 unique lipid resids for each row 
	#---...where each row is a ranking of the lipid atoms nearest an ion
	#---...the result is the sorted-atom index for the closest 3 unique lipid resids
	partners_atoms_sorted = array([sort(unique(a,return_index=True)[1])[:nrank] for a in lipid_resids[aa]])
	#---recover the atom indices for the partners from which we determine resid and resname and distances
	partners_atoms = array([aa[i][partners_atoms_sorted[i]] for i in range(nions)])
	#---given the correct atom indices we can figure out the important identities and distances of partners
	lipid_distances = array([cd[i,partners_atoms[i]] for i in range(nions)])	
	return (lipid_distances,partners_atoms)
	
def zonelister(data,lims):

	"""
	Turns a trajectory into a discrete trajectory within particular limits. The dims input provides the 
	shape of the data.
	"""

	dims = shape(data)
	return argmax((data>
		transpose(tile(lims,(dims[0],dims[1],1)),(2,0,1)))*\
		transpose(tile(arange(len(lims)),(dims[0],dims[1],1)),(2,0,1)),
		axis=0)

