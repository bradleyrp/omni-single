#!/usr/bin/env python

from numpy import *

"""
Plotting functions for creating undulation spectra.
"""

def fftwrap(dat): 
	
	"""
	Wrapper function for 2D Fourier transform.
	"""

	return fft.fft2(array(dat))

def perfect_collapser(xs,ys,trim=False):

	"""
	Return two arrays for a "perfect" collapse.
	"""
	
	xsort = array(sort(list(set(xs))))
	if trim: xsort = xsort[1:]
	inds = argmax(array([xs==xsort[i] for i in range(len(xsort))]),axis=0)
	#---! future warning below
	if type(ys)==type(None): col = None
	else: col = array([mean(ys[where(inds==i)]) for i in range(len(xsort))])
	return xsort,col,inds

def blurry_binner(xs,ys,bin_width=0.05,trim=True):
	
	"""
	Group wavevectors by bins.
	"""
	
	blurred = (xs/bin_width).astype(int)
	xsort = array(sort(list(set(blurred))))
	if trim: xsort = xsort[1:]
	inds = argmax(array([(xs/bin_width).astype(int)==xsort[i] for i in range(len(xsort))]),axis=0)
	if type(ys)!=ndarray: coly = None
	else: coly = array([mean(ys[where(inds==i)]) for i in range(len(xsort))])
	colx = array([mean(xs[where(inds==i)]) for i in range(len(xsort))])
	return colx,coly,inds

def calculate_undulations_testing(surf,vecs,chop_last=False,lims=(0.0,1.0),perfect=False,raw=False):

	"""
	Compute undulation spectrum.
	"""

	nframes,m,n = shape(surf)
	frames = arange(nframes)
	Lx,Ly = mean(vecs,axis=0)[:2]
	lenscale = 1.
	qmagsshift = lenscale*array([[sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	surf = surf-mean(surf)
	hqs = array([fftwrap(surf[fr])/lenscale/double((m*n)) for fr in range(nframes)])
	y = reshape(mean(real(hqs)**2,axis=0),-1)
	x = reshape(qmagsshift,-1)
	#---!!! blurry binner appears broken
	if not perfect: x2,y2 = blurry_binner(x,x)[1],blurry_binner(x,y)[1]
	else: x2,y2 = perfect_collapser(x,x)[1],perfect_collapser(x,y)[1] 
	#import pdb;pdb.set_trace()
	x3,y3 = x2[1:(-1 if chop_last else None)],y2[1:(-1 if chop_last else None)]
	#---raw option to get back to basics
	if raw: x3,y3 = x,y
	goodslice = where(all((x3>lims[0],x3<lims[1]),axis=0))
	#---! error below x2 instead of x3 now fixed
	#import pdb;pdb.set_trace()
	#---dropped: kappa = mean(1/((x3[1:]*x3[1:]**4)[goodslice]*Lx*Ly/lenscale**2))
	kappa = mean(1/((y3*x3**4)[goodslice]*Lx*Ly/lenscale**2))
	#import pdb;pdb.set_trace()
	return {'y':y3,'x':x3,'kappa':kappa}

def calculate_undulations(surf,vecs,chop_last=False,lims=(0,1.0),perfect=False,raw=False):

	"""
	Compute undulation spectrum.
	"""

	nframes,m,n = shape(surf)
	frames = arange(nframes)
	Lx,Ly = mean(vecs,axis=0)[:2]
	lenscale = 1.
	qmagsshift = lenscale*array([[sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	surf = surf-mean(surf)
	hqs = array([fftwrap(surf[fr])/lenscale/double((m*n)) for fr in range(nframes)])
	y = reshape(mean(abs(hqs)**2,axis=0),-1)
	x = reshape(qmagsshift,-1)
	if not perfect: x2,y2 = blurry_binner(x,x)[1],blurry_binner(x,y)[1]
	else: x2,y2 = perfect_collapser(x,x)[1],perfect_collapser(x,y)[1] 
	x3,y3 = x2[1:(-1 if chop_last else None)],y2[1:(-1 if chop_last else None)]
	if lims!=None: goodslice = where(all((x3>lims[0],x3<lims[1]),axis=0))
	else: goodslice = arange(len(x3)).astype(int)
	kappa = mean(1/((y2[1:]*x2[1:]**4)[goodslice]*Lx*Ly/lenscale**2))
	return {'y':y3,'x':x3,'kappa':kappa}
