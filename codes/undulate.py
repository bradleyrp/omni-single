#!/usr/bin/env python

import re
import numpy as np
import scipy
import scipy.optimize
from base.tools import status

machine_eps = eps = np.finfo(float).eps
dotplace = lambda n : re.compile(r'(\d)0+$').sub(r'\1',"%3.5f"%float(n)).rjust(8)

"""
Plotting functions for creating undulation spectra.
"""

def fftwrap(dat): 
	"""
	Wrapper function for 2D Fourier transform.
	"""
	return np.fft.fft2(np.array(dat))

def perfect_collapser(xs,ys,trim=False,return_indices=False):
	"""
	Return two arrays for a "perfect" collapse.
	"""
	xsort = np.array(np.sort(list(set(xs))))
	if trim: xsort = xsort[1:]
	inds = np.argmax(np.array([xs==xsort[i] for i in range(len(xsort))]),axis=0)
	#---! future warning below
	if type(ys)==type(None): col = None
	else: col = np.array([np.mean(ys[np.where(inds==i)]) for i in range(len(xsort))])
	if return_indices: return xsort,col,inds
	else: return xsort,col

def blurry_binner_deprecated(xs,ys,bin_width=0.05,trim=True,return_mapping=False):
	"""
	Group wavevectors by bins.
	Retired this because it was confusing!
	"""
	blurred = (xs/bin_width).astype(int)
	xsort = np.array(np.sort(list(set(blurred))))
	if trim: xsort = xsort[1:]
	inds = np.argmax(np.array([blurred==xsort[i] for i in range(len(xsort))]),axis=0)
	if type(ys)!=np.ndarray: coly = None
	else: coly = np.array([np.mean(ys[np.where(inds==i)]) for i in range(len(xsort))])
	colx = np.array([np.mean(xs[np.where(inds==i)]) for i in range(len(xsort))])
	if not return_mapping: return colx,coly,inds
	else:
		#---! this takes a moment
		mapping = [np.where(blurred==i)[0] for i in np.sort(blurred)]
		return colx,coly,inds,mapping

def blurry_binner(xs,ys,bin_width=0.05,trim=True,return_mapping=False,mode='snapped'):
	"""
	Improved over original with sensible binning and interpolation.
	Removed the index return strategy in favor of handling things internally.
	"""
	#---get the range of the independent variable
	x0,x1 = xs.min(),xs.max()
	#---snap to the bin width
	x0b,x1b = bin_width*np.floor(x0/bin_width),bin_width*np.ceil(x1/bin_width)
	#---develop canonical bin markers
	bins = np.arange(x0b,x1b+bin_width,bin_width)
	#---bins are represented in the middle
	xmid = (bins[1:]+bins[:-1])/2.
	#---get the bin positions
	xbin = np.floor(xs/bin_width).astype(int) - np.floor(xs.min()/bin_width).astype(int)
	#---check that everything is in the right bin
	try: 
		if not (np.all(xs<=bins[1:][xbin]) and np.all(xs>=bins[:-1][xbin])):
			import ipdb;ipdb.set_trace()
			raise Exception('binning failure!')
	except:
		import ipdb;ipdb.set_trace()
	#---reindex each position into bins
	indices = np.unique(xbin)
	reindex = [np.where(xbin==i)[0] for i in indices]
	#---snapped mode just uses the middle of each bin and takes the average of the constituents
	if mode=='snapped':
		y_out = np.array([ys[r].mean() for r in reindex])
		x_out = xmid[indices]
	#---! recommend developing an interpolation method here to shift laterally if off-center constituents
	else: raise 
	#---for each point, return the indices of other points in the same bin
	#---...note that this is essential for developing weights in the blurry_explicit
	#---...where the weights end up: weights = 1./np.array([len(i) for i in q_mapping])[band]
	if return_mapping: 
		mapping = [np.where(xbin==i)[0] for i in np.sort(xbin)]
		return x_out,y_out,mapping
	else: return x_out,y_out

def undulation_fitter(q_raw,hqs,area,initial_conditions=(20.0,0.0),residual_form='linear'):
	"""Like bath fitter, but for undulations."""

	if residual_form == 'log':
		def residual(values): 
			return np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
	elif residual_form == 'linear': 
		def residual(values): 
			return np.sum((values-1.0)**2)/float(len(values))
	else: raise Exception('unclear residual form %s'%residual_form)

	def multipliers(x,y): 
		"""Multiplying complex matrices in the list of terms that contribute to the energy."""
		return x*np.conjugate(y)

	def callback(args):
		"""Watch the optimization."""
		global Nfeval
		name_groups = ['kappa','gamma']
		text = ' step = %d '%Nfeval+' '.join([name+' = '+dotplace(val)
			for name,val in zip(name_groups,args)+[('error',objective(args))]])
		status('searching! '+text,tag='optimize')
		Nfeval += 1

	def objective(args,mode='residual'):
		kappa,gamma = args
		termlist = [multipliers(x,y) for x,y in [(hqs,hqs)]]
		gamma = 0.0
		ratio = kappa/2.0*area*(termlist[0]*q_raw**4)+gamma/2.0*area*(termlist[0]*q_raw**2)
		if mode=='residual': return residual(ratio)
		elif mode=='ratio': return ratio
		else: raise Exception('invalid mode %s'%mode)

	global Nfeval
	Nfeval = 0
	test_ans = objective(initial_conditions)
	if not isinstance(test_ans,np.floating): 
		raise Exception('objective_residual function must return a scalar')
	fit = scipy.optimize.minimize(objective,
		#---note Newton-CG requires Jacobian, coblya is not strong enough, usw
		x0=tuple(initial_conditions),method='BFGS',callback=callback)
	return dict(fit,kappa=fit.x[0],gamma=fit.x[1])

def calculate_undulations(surf,vecs,fit_style=None,chop_last=False,lims=(0,1.0),
	perfect=False,raw=False,midplane_method=None,custom_heights=None,residual_form='log',fit_tension=False):
	"""
	Compute undulation spectrum.
	"""
	nframes,m,n = surf.shape
	frames = np.arange(nframes)
	Lx,Ly = np.mean(vecs,axis=0)[:2]
	lenscale = 1.
	qmagsshift = lenscale*np.array([[np.sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	area = (Lx*Ly/lenscale**2)

	#---remove the supposed "average" structure
	if midplane_method=='flat' or midplane_method==None:
		surf = surf-np.mean(surf)
	elif midplane_method=='average':
		surf = surf-np.mean(surf,axis=0)
	elif midplane_method=='average_normal' and type(custom_heights)==type(None):
		raise Exception('send custom_heights for average_normal')
	elif midplane_method=='average_normal': surf = custom_heights
	else: raise Exception('invalid midplane_method %s'%midplane_method)

	#---perform the FFT and line up the results
	hqs = np.array([fftwrap(surf[fr])/lenscale/np.double((m*n)) for fr in range(nframes)])
	y = np.reshape(np.mean(np.abs(hqs)**2,axis=0),-1)
	x = np.reshape(qmagsshift,-1)

	#---default method
	if fit_style==None: fit_style = 'band,perfect,simple'

	packed = {}
	#---choose a binning method, range method, and fitting method
	if fit_style in ['band,perfect,simple','band,perfect,basic',
		'band,perfect,fit','band,perfect,curvefit','band,blurry,curvefit','band,perfect,curvefit-crossover']:

		if lims==None: raise Exception('fit_style %s requires lims'%fit_style)
		#---collapse, perfectly
		x2,y2 = perfect_collapser(x,x)[1],perfect_collapser(x,y)[1] 
		goodslice = np.where(np.all((x2>lims[0],x2<lims[1]),axis=0))
		#---sample in the band and drop the zero mode
		#---! where is the zero mode dropped?
		x3,y3 = x2[goodslice],y2[goodslice]

		#---apply binning method
		if fit_style=='band,blurry,curvefit':
			x3,y3 = blurry_binner(x3,y3)


		#---perform the fit
		if fit_style=='band,perfect,simple': 
			#---the simple method is a crude way to do the fit, by fixing the exponent and then using the fact
			#---...that the y-axis centroid of the positions must determine the 
			#---we recommend against using this method for further calculaiton. useful as a comparison only.
			if fit_tension: raise Exception('cannot do simple with fit_tension')
			kappa,gamma = 2.0*np.mean(1/((y3*x3**4)*Lx*Ly/lenscale**2)),0.0
		elif fit_style=='band,perfect,basic': 
			if fit_tension: raise Exception('cannot do basic with fit_tension')
			#---in this method we use polyfit to do a linear fit in the log space
			#---note that exponent is free for this method and we have no specific residual form
			c0,c1 = np.polyfit(np.log10(x3),np.log10(y3),1)
			#---fit would be: 10**c1*x3**(c0)
			#---! just spitballing here
			kappa,gamma = 1./(10**c1*area)*2.0,0.0
			#---save for debugging
			packed['linear_fit_in_log'] = dict(c0=c0,c1=c1)
		elif fit_style in ['band,perfect,curvefit','band,blurry,curvefit']:
			#---in this method we use the scipy curve_fit function
			exponent = 4.0
			kwargs = {}
			kwargs.update(bounds=((5.0,-1.0),(10**2.0,1.0)))
			kwargs.update(maxfev=10**5)
			if residual_form=='linear':
				def hqhq(q_raw,kappa,sigma):	
					if not fit_tension: sigma = 0.0
					return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
				fit = scipy.optimize.curve_fit(hqhq,x3,y3,**kwargs)
				kappa,gamma = fit[0][0],fit[0][1]
				print(gamma)
			elif residual_form=='log':
				#---! deprecated but works
				if False:
					def hqhq(q_raw,kappa,sigma):	
						sigma = 0.0
						return np.log10(2.0/area/kappa)+-4.0*q_raw
					fit = scipy.optimize.curve_fit(hqhq,np.log10(x3),np.log10(y3),**kwargs)
				#---new method has the full expression fit in the log space
				else:
					def hqhq(q_raw,kappa,sigma):	
						if not fit_tension: sigma = 0.0
						return np.log10(1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2)))
					fit = scipy.optimize.curve_fit(hqhq,x3,np.log10(y3),**kwargs)
				kappa,gamma = fit[0][0],fit[0][1]
			else: raise Exception('invalid residual_form %s'%residual_form)

		elif fit_style=='band,perfect,curvefit-crossover':
			#---in this method we use the scipy curve_fit function
			exponent = 4.0
			kwargs = {}
			kwargs.update(bounds=((5.0,-1.0,lims[0]+0.001),(10**2.0,1.0,lims[1]-0.001)))
			kwargs.update(maxfev=10**5)
			if residual_form!='log': raise Exception('crossover only works with residual_form log')
			if fit_tension==False: raise Exception('crossover only works with fit_tension')
			def hqhq(q_raw,kappa,sigma,crossover):	
				if sigma==0.0: sigma = eps
				regime_tension = q_raw<=crossover
				regime_bending = ~regime_tension
				heights_tension = np.log10(1.0/(area/2.0*(sigma*q_raw[regime_tension]**2.0)))
				heights_bending = np.log10(1.0/(area/2.0*(kappa*q_raw[regime_bending]**4.0)))
				return np.concatenate((heights_tension,heights_bending))
			fit = scipy.optimize.curve_fit(hqhq,x3,np.log10(y3),(20.0,0.001,0.02),**kwargs)
			kappa,gamma,crossover = fit[0][0],fit[0][1],fit[0][2]
			packed.update(crossover=crossover)
		elif fit_style=='band,perfect,fit':
			#---the most advanced method uses scipy.optimmize in a separate function
			fit = undulation_fitter(x3,y3,area,residual_form=residual_form,
				initial_conditions=(20.0,0.0))
			kappa,gamma = fit['kappa'],fit['gamma']
		else: raise Exception('invalid fit_style %s'%fit_style)
		#---return the data
		packed.update(kappa=kappa,points=np.transpose((x3,y3)),sigma=gamma,
			q_raw=x,energy_raw=y,q_binned=x2,energy_binned=y2,area=area)

	else: raise Exception('invalid fit_style %s'%fit_style)
	return packed
