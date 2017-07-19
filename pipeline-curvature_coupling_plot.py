#!/usr/bin/env python

import scipy.ndimage
import scipy.interpolate
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects

#---settings
n_best = 20
cbar_nlevels = 5
cut_extent = 0,18
contour_lineskip = 4
contour_nlevels = 100
contour_interp_pts = 100
cut_curve = tweak['cut_curve']
error_ceiling = tweak.get('error_ceiling',None)
error_floor = None
routine = ['summary']

#---COLLECT
#-------------------------------------------------------------------------------------------------------------

def prepare_postdat(hypo_groups,**kwargs):
	"""
	Compile the data required to plot landscapes, spectra, etc.
	"""
	cut_ce = kwargs.get('cut_ce',None)
	contour_interp_pts = kwargs['contour_interp_pts']
	postdat,comparisons = {},[]
	#---group key indexes the curvature field comparison: protein/single, static/dynamic
	for group_key,tweak_cf_spec in enumerate(tweak_cf):
		#---in case you are running computations, we filter here
		hypo_groups = {k:[h for hh,h in enumerate(v) 
			if os.path.isfile(namer_cc(sessions['hypothesis'].query(Hypothesis).filter_by(
				**Hypothesis(**h).base()).one().id))] for k,v in hypo_groups.items()}
		#---organize the comparisons
		if tweak_cf_spec['mapping'] == 'single': 
			compspec = [('single',sn,'none') for sn in sns_all]
		elif tweak_cf_spec['mapping'] == 'protein': 
			control_sns = [sn for sn in sns_all if work.meta[sn].get('nprots',1)==0]
			if control_sns:
				sn_free, = [sn for sn in sns_all if work.meta[sn].get('nprots',1)==0]
				compspec = [('protein',sn_free,sn) for sn in sns_unfree]
			#---! former typo? compspec = [('protein',sn,'none') for sn in sns_unfree]
			compspec = [('protein',sn,'none') for sn in sns_all]
		else: raise
		comparisons.append(compspec)
		#---collate the results looping over comparisons
		for compnum,(mapping,sn,fallback) in enumerate(compspec):
			hypo_groups_sub = [h for hh,h in enumerate(hypo_groups[group_key]) 
				if h['sn']==sn and h['fallback']==fallback and h['mapping']==mapping]
			curvatures = sorted(set([i['curvature'] for i in hypo_groups_sub]))
			extents = sorted(set([i['sigma_a'] for i in hypo_groups_sub]))
			curvatures_all,extents_all = list(curvatures),list(extents)
			c0,c1 = min(curvatures),max(curvatures)
			e0,e1 = min(extents),max(extents)
			#---intervene here to apply a filter on curvature and extent
			#---the curvature/extent filter will apply for all downstream calculations
			if cut_ce: 
				c0,c1 = max(c0,cut_ce[0][0]),min(c1,cut_ce[0][1])
				e0,e1 = max(e0,cut_ce[1][0]),min(e1,cut_ce[1][1])
				curvatures = [c for c in curvatures if c>=c0 and c<=c1]
				extents = [e for e in extents if e>=e0 and e<=e1]
			#---prepare the raw error map
			land = 0*np.ones((len(curvatures_all),len(extents_all)))
			kappaland = -1*np.ones((len(curvatures_all),len(extents_all)))
			gammaland = -1*np.ones((len(curvatures_all),len(extents_all)))
			vibeland = -1*np.ones((len(curvatures_all),len(extents_all)))
			survey = {}
			for cc,c in enumerate(curvatures_all):
				for ee,e in enumerate(extents_all):
					subset = [h for hh,h in enumerate(hypo_groups_sub) 
						if h['sn']==sn and h['curvature']==c and h['sigma_a']==e]
					hypo = subset[0]
					#hacccccckkked hypo, = subset
					reduced = Hypothesis(**hypo).base()
					row = sessions['hypothesis'].query(Hypothesis).filter_by(**reduced).one()
					land[cc,ee] = row.error
					kappaland[cc,ee] = row.kappa
					gammaland[cc,ee] = row.gamma
					vibeland[cc,ee] = row.vibe
					survey[row.id] = {'kappa':row.kappa,'gamma':row.gamma,'vibe':row.vibe,'error':row.error}
			errormap_raw = land
			#---! extremly hacked. not sure why that pixel is not filled in?
			errormap_raw[errormap_raw==-1.0] = errormap_raw.max()
			postdat[(mapping,sn,fallback)] = {'errormap_raw':errormap_raw,
				'survey':survey,'land_kappa':kappaland,'land_gamma':gammaland,'land_vibe':vibeland}
			#---interpolate contour plots
			curvature_extent_error = np.array([(c,e,land[cc,ee]) 
				for cc,c in enumerate(curvatures) for ee,e in enumerate(extents)])
			finex = np.linspace(c0,c1,contour_interp_pts)
			finey = np.linspace(e0,e1,contour_interp_pts)
			grid_x,grid_y = np.meshgrid(finex,finey)
			errormap = scipy.interpolate.griddata(curvature_extent_error[:,:2],curvature_extent_error[:,2],
				(grid_x,grid_y),method='cubic')
			postdat[(mapping,sn,fallback)]['errormap'] = errormap
			postdat[(mapping,sn,fallback)]['curvatures'] = np.array(curvatures)
			postdat[(mapping,sn,fallback)]['extents'] = np.array(extents)
			#---pull some of the spectra
			errormap_raw_reduced = errormap_raw[np.array([c in curvatures for c in curvatures_all]),:][:,
				np.array([e in extents for e in extents_all])]
			choice_cuts = np.array([np.unravel_index(i,errormap_raw_reduced.shape) 
				for i in np.reshape(errormap_raw_reduced,-1).argsort()[:n_best]])
			st = time.time()
			for counter,key in enumerate(choice_cuts):
				status('pulling spectra',i=counter,looplen=len(choice_cuts),start=st,tag='reload')
				c,e = curvatures[key[0]],extents[key[1]]
				#---rhymes with the preparation of the errormap_raw above
				subset = [h for hh,h in enumerate(hypo_groups_sub) 
					if h['sn']==sn and h['curvature']==c and h['sigma_a']==e]
				try: hypo, = subset
				except:
					import ipdb;ipdb.set_trace()
				reduced = Hypothesis(**hypo).base()
				row = sessions['hypothesis'].query(Hypothesis).filter_by(**reduced).one()
				pk,error = row.id,row.error
				incoming = load(os.path.basename(namer_cc(pk)),cwd=os.path.dirname(namer_cc(pk)))
				postdat[(mapping,sn,fallback,'choice_cuts',counter)] = \
					(incoming['wavevectors'],incoming['best_energy'],error)
	postdat['comparisons'] = comparisons
	return postdat

def prepare_error_wells(errormap,weight_method,well_topcut=1.05,review=False,**kwargs):
	"""
	Use a well-finder with an error-weighter to find the "best" hypothesis on the contour plot
	along with its error bars and its nearest actually-observed best hypothesis.
	"""
	packed = {}
	dat = np.array(errormap)
	#---construct different error weighting functions
	elim0,elim1 = dat.min(),(dat.min()*well_topcut)
	weighters = {
		#---center of geometry for all points in the domain
		'cog':lambda x : 1.0 if x>=elim0 and x<= elim1 else 0.0,
		#---points are weighted by relative distance from the bottom
		'linear_interp':lambda x : (x-elim1)/(elim1-elim0),
		'min': lambda x : 1.0 if x==elim0 else 0.0,
		}
	vals = np.transpose(np.where(dat<=(dat.min()*well_topcut)))
	domains = np.zeros((contour_interp_pts,contour_interp_pts))
	domains[(vals.T[0],vals.T[1],)] = 1.0
	im,number_of_objects = scipy.ndimage.label(domains)

	#---convert contour back to regular units
	curvatures,extents = [postdat[pkey][i] for i in ['curvatures','extents']]
	c0,c1 = min(curvatures),max(curvatures)
	e0,e1 = min(extents),max(extents)
	fine_to_coarse_c = lambda x : float(x)/contour_interp_pts*(c1-c0)+c0
	fine_to_coarse_e = lambda x : float(x)/contour_interp_pts*(e1-e0)+e0
	blobs = scipy.ndimage.find_objects(im)
	vals,counts = np.unique(im,return_counts=True)
	#---the following code takes the second largest domain which should be the largest well
	#---the largest well might not contain the minimum, hence in the 'min' weighting method we pick different
	if weight_method == 'min' and not force_large_well: 
		best_well_domain = float(im[np.unravel_index(dat.argmin(),dat.shape)])
	else: best_well_domain = vals[counts.argsort()[-2]]
	dom = 1.0*(im==best_well_domain)
	possible_bottoms = {name:{} for name in weighters}
	#---we no longer loop over all weight methods
	if not force_large_well:
		weighter = weighters[weight_method]
	else:
		#---typically the weighters are based on the global minimum
		#---the force_large_well flag resets the apparent minimum to the minimum of the largest well
		#---...this is essential in the case of the 8xENTH system where the cubic interpolation on the 
		#---...contour plot causes a local (non-global minimum) in the raw data to be the apparent minimum
		#---...on the contour plot. the second well is valid but not globally minimum
		elim0 = dat[np.where(dom)].min()
		weighter = lambda x : 1.0 if x==elim0 else 0.0
	possible_bottoms[weight_method]['inds'] = np.transpose(np.where(dom))
	weights = np.array(map(weighter,dat[np.where(dom)]))
	if sum(weights)==0.0:
		print 'weights are zero'
		import pdb;pdb.set_trace()
	wwb = np.average(possible_bottoms[weight_method]['inds'],
		weights=weights,axis=0)[::-1]
	possible_bottoms[weight_method]['values'] = wwb
	cx = fine_to_coarse_c(possible_bottoms[weight_method]['values'][0])
	ex = fine_to_coarse_e(possible_bottoms[weight_method]['values'][1])
	best_coarse_inds,best_coarse = zip(*[(k,j[k]) for i,j in 
		zip([cx,ex],[curvatures,extents]) for k in np.argsort(np.abs(i-j))[:1]])
	#---get the nearest contour plot points
	best_fine_inds = np.round(possible_bottoms[weight_method]['values']).astype(int)
	dompts = np.transpose(np.where(dom))
	#---the following two lines are effectively lazy magic
	lims = [[(i.min(),i.max()) for i in 
		[dompts[np.where(dompts[:,j]==best_fine_inds[1-j])[0],1-j]]][0] for j in range(2)]
	error_lims = [(f(i),f(j)) for f,(i,j) in zip(*[[fine_to_coarse_c,fine_to_coarse_e],lims])]

	if review:
		ax = plt.subplot(131)
		ax.imshow(dat,interpolation='nearest',origin='lower',vmin=dat.min(),
			vmax=error_ceiling,extent=(0,contour_interp_pts,0,contour_interp_pts))
		ax.scatter(vals.T[1],vals.T[0],c='k',marker='.')
		ax.set_xlim(0,contour_interp_pts);ax.set_ylim(0,contour_interp_pts)
		ax = plt.subplot(132)
		ax.imshow(im,interpolation='nearest',origin='lower')
		ax = plt.subplot(133)
		ax.imshow(dom,interpolation='nearest',origin='lower')
		#---previously looped over weight_method
		wwb = possible_bottoms[weight_method]['values']
		ax.scatter([wwb[0]],[wwb[1]],marker='o',s=50,c='r')
		ax.set_xlim(0,contour_interp_pts);ax.set_ylim(0,contour_interp_pts)
		#---select a weighter
		plt.show()
	#---store the fine and coarse positions on the error landscapes for making error bars
	packed['error_lims_positions_fine'] = lims
	packed['error_lims_positions_coarse'] = [[jj(i) for i in lims[kk]] 
		for kk,jj in zip(range(2),[fine_to_coarse_c,fine_to_coarse_e])]
	packed['error_lims'] = error_lims
	packed['possible_bottoms'] = possible_bottoms
	packed['best_fine'] = cx,ex
	packed['best_fine_inds'] = best_fine_inds
	packed['best_coarse'] = best_coarse
	packed['best_coarse_inds'] = best_coarse_inds
	packed['fine_to_coarse'] = fine_to_coarse_c,fine_to_coarse_e
	packed['domain'] = dom
	return packed

def get_curvature_field(pkey):
	"""
	Look up a best-fitting curvature field. Loads the field into memory.
	"""
	best_c,best_e = memory[pkey+('packed',)]['best_coarse']
	lookup = {'mapping':pkey[0],'motion':tweak_cf_spec['motion'],'isotropy':1.0,
		'sigma_a':best_e,'sn':pkey[1],'fallback':pkey[2]}
	pk = sessions['field'].query(Field).filter_by(**Field(**(lookup)).present()).one().id
	key_cf = ('curvature',pk)
	if key_cf not in memory: 
		fn_cf = namer_cf(pk)
		memory[key_cf] = load(os.path.basename(fn_cf),rootdir_cf)
	return {'best_c':best_c,'best_e':best_e,'cfid':key_cf}

def fit_cooperativity(ax,plot_placer_sigmoid,fs=10):
	"""
	Perform the cooperativity fit and plot it.
	Requires an axis for plotting and a list of the pkeys in sns_all order.
	"""
	def volume_gauss2d(fac,sig,nprots):
		#---issue with the integral!
		if False: return np.pi*fac*sig**2*nprots/2.
		return 2*np.pi*fac*sig**2*nprots
	#---independent variable is the number of proteins
	xs = np.array([work.meta[sn].get('nprots',1) for sn in sns])
	#---dependent variable is the bending
	bends = ys = np.zeros(len(xs))
	for snum,sn in enumerate(sns_all):
		pkey = plot_placer_sigmoid[snum]
		nprots = xs[snum]
		cf_details = get_curvature_field(pkey)
		c,e = [cf_details['best_%s'%s] for s in 'ce']
		bends[snum] = volume_gauss2d(e,c,nprots)
	ys_normed = ys/ys.max()

	def binding(x,n,k):
		return n*k*numpy.power(x,n-1)/(1.0+n*k*numpy.power(x,n-1))
	def binding_obs(n,k):
		def kernel(n,k):
			return binding(xs,n,k)
		return kernel(n,k)
	def residual(args):
		n,k = args
		return np.mean(np.power(binding_obs(n,k)-ys_normed,2))
	def binding(x,n,k):
		return n*k*np.power(x,n-1)/(1.0+n*k*np.power(x,n-1))
	def binding_obs(n,k,):
		def kernel(n,k):
			return binding(xs,n,k)
		return kernel(n,k)
	def residual4(args):
		n,k,cinf,czero = np.abs(args)
		return np.mean(np.power(binding_obs(n,k)-(ys-czero)/(cinf-czero),2))

	#---optimizer constraints are fairly stable
	constraints = [{'type':'ineq','fun':lambda x: x[0]},
		{'type':'ineq','fun':lambda x: x[1]},
		{'type':'ineq','fun':lambda x: xs[-1]-x[1]**(-1./(x[0]-1))},
		{'type':'ineq','fun':lambda x: 0.005-x[3]},
		{'type':'ineq','fun':lambda x: x[3]}]
	#---different optimizations for different wellfinder choices
	if weight_method == 'cog' and not force_large_well:
		initial_guess = [3.0,0.7,0.03,0]
	elif weight_method == 'min' and not force_large_well:
		#---gives weird 2 cooperativity
		initial_guess = [4.0,10**-5,1.0,0]
		#---gives slightly better 6.5 cooperativity
		initial_guess = [4.0,0.7,1.0,0]
	elif weight_method == 'min' and force_large_well:
		initial_guess = [3.0,0.7,0.03,0]
	else: raise Exception('no optimization method here')
	#---! fixed volume integral
	initial_guess = [5.0,10**-31,1.0,0]
	fit4 = scipy.optimize.minimize(
		residual4,x0=np.array(initial_guess),
		options={'disp':True},constraints=constraints)

	n_fit4,k_fit4,cinf,czero = np.abs(fit4['x'])
	answer = {'n':n_fit4,'K':k_fit4,'Cinf':cinf,'C0':czero}
	print('[RESULT] %s'%fit4)

	#---sigmoid plotting
	xmin,xmax = (xs.min()-(xs.max()-xs.min())*0.1,xs.max()*1.2)
	finex = np.arange(xmin,xmax,(xmax-xmin)/100.)
	#---no need for fitted points to be emphasized since the line is there
	if False:
		ax.plot(xs,binding_obs(n_fit4,k_fit4),'o-',lw=0,mew=0,c='k',ms=10,zorder=3)
		ax.plot(xs,binding_obs(n_fit4,k_fit4),'o-',lw=0,mew=0,c='w',ms=15,zorder=2)
	ax.plot(finex,binding(finex,n_fit4,k_fit4),'-',lw=2,mew=0,c='k',ms=10,zorder=1)
	for altn in [2,3,4,5]:
		ax.plot(finex[np.where(finex>=0.0)],
			binding(finex,altn,k_fit4)[np.where(finex>=0.0)],'--',lw=1,mew=0,c='k',ms=10)
	#---plot the actual points
	ax.plot(xs,(ys-czero)/(cinf-czero),'o-',lw=0,mew=0,ms=5,c='w',zorder=3)
	ax.plot(xs,(ys-czero)/(cinf-czero),'o-',lw=0,mew=0,ms=10,c='k',zorder=2)
	ax.set_xlim(xmin,xmax)
	eqn = (r"\begin{eqnarray*} n_p &=& %.2f \\ K_n &=& %.1f \times 10^{3} \\ D_s(0) "+\
		r"&=& %.1f \\ D_s({\infty}) &=& %.1f \times 10^{2} \end{eqnarray*}")%(
		answer['n'],answer['K']*10**3,answer['C0'],answer['Cinf']*10**2)
	ax.text(0.95,0.05,eqn,transform=ax.transAxes,size=fs-1,verticalalignment='bottom',ha='right',
		bbox=dict(boxstyle='square,pad=0.4',facecolor='gainsboro',
		alpha=1.0,lw=1))
	ax.set_ylim(-0.1,1.1)
	ax.grid(True,linestyle='-',zorder=0,alpha=0.5)
	ax.tick_params(axis='y',which='both',left='off',right='off')
	ax.tick_params(axis='x',which='both',bottom='off',top='off')
	ax.set_axisbelow(True)
	if False: ax.set_title('bending',fontsize=fs+4)
	ax.set_xlabel('number of proteins',fontsize=fs+2)
	ax.set_xticklabels(['%d'%i for i in ax.get_xticks()],fontsize=fs+2)
	ax.set_yticklabels(ax.get_yticks(),fontsize=fs+2)
	ax.set_ylabel(
		r'$\mathrm{ \frac{ {D}_{s}({n}_{p}) - {D}_{s}(0)  }{ {D_s}({\infty}) - {D}_{s}(0) } }$',
		fontsize=fs+5)
	#---custom
	ax.set_xlim(xmin,xmax)

#---PLOT UTILITIES
#-------------------------------------------------------------------------------------------------------------

def tickoff(ax,grid=True,color='0.5',alpha=1.0):
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	ax.set_axisbelow(True)
	if grid:
		ax.yaxis.grid(b=True,which='both',color=color,linestyle='-',alpha=alpha)
		ax.xaxis.grid(b=True,which='both',color=color,linestyle='-',alpha=alpha)
	
def titlemaker(sn,fallback):
	return work.meta[sn]['label']+(
		'' if fallback=='none' else '\n(s.t. '+work.meta[fallback]['label']+')')

def clean_powers(vmin,vmax,strict=True,nticks=5,ticks=None):
	if not type(ticks)!=type(None): ticks = np.linspace(vmin,vmax,nticks)
	else: nticks = len(ticks)
	if strict: powlev = min([int(('%.2e'%t).split('e')[1]) for t in ticks])
	else: powlev = int(('%.2e'%ticks[-1]).split('e')[1])
	ticklabels = ['%.f'%(i/10**powlev) for i in ticks]
	powtext = r' $\times 10^{%d}$'%(-1*powlev)
	return ticks,ticklabels,powtext

def letters_above(ax,pts,letters,ymax,fs=10):
	letter_bbox_props = dict(boxstyle='square,pad=0.3',alpha=0)
	for letter,(x,y) in zip(letters,pts):
		tb = ax.text(x,y+(ymax/1.2)/20.,letter,ha='center',va='bottom',
			rotation=0,color='k',fontsize=fs)

#---PLOT FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

def whisker(ax,x,y,error=None,c='w',rc='k',lc='k',s=15,lw=2,rlw=3,
	ring=False,zorder=3,stance=0.2,no_center=False): 
	"""
	Plot a single whisker with lots of options.
	"""
	if not no_center:
		ax.scatter([x],[y],color=c,alpha=1.,marker='o',s=s,lw=0,
			zorder=zorder+ring+(not not error),clip_on=False)
	if ring and not no_center: 
		if False: ax.scatter([x],[y],color=rc,alpha=1.,marker='o',lw=0,s=s*4,zorder=zorder+(not not error))
		else: ax.scatter([x],[y],color=rc,alpha=1.,marker='o',
			lw=rlw,s=s*1,zorder=zorder+(not not error),clip_on=False)
	if error:
		ax.plot([x,x],error,color=lc,lw=lw,solid_capstyle='round',zorder=zorder,clip_on=False)
		for e in error: ax.plot([x+stance,x-stance],[e,e],color=lc,
			lw=lw,solid_capstyle='round',zorder=zorder,clip_on=False)

def accounting_stripev(cats,ylims=(0,1),width=1.0):
	"""
	Alternating vertical gray background stripes.
	"""
	for ii,i in enumerate([i%2 for i in cats]):
		ax.add_patch(mpl.patches.Rectangle(((ii-0.5)*width,ylims[0]),
			width,ylims[1],alpha=0.05,color=('k' if i==0 else 'w'),lw=0))

def plot_errormap(ax,pkey,fs=10,**kwargs):
	"""
	Plot the error landscape.
	"""
	mapping,sn,fallback = pkey
	error_max = kwargs['error_max']
	error_min = kwargs['error_min']
	cut_extent = kwargs['cut_extent']
	cut_curve = kwargs['cut_curve']
	show_best = kwargs.get('show_best',True)
	old_show_best = kwargs.get('old_show_best',False)
	best_coarse = kwargs.get('best_coarse',True)
	under_color = kwargs.get('under_color',None)
	levels = np.linspace(error_min,error_max,contour_nlevels)

	errormap = np.array(postdat[pkey]['errormap'])
	curvatures,extents = [postdat[pkey][key] for key in ['curvatures','extents']]
	c0,c1 = min(curvatures),max(curvatures)
	e0,e1 = min(extents),max(extents)
	finex,finey = np.linspace(c0,c1,contour_interp_pts),np.linspace(e0,e1,contour_interp_pts)
	X,Y = np.meshgrid(finex,finey)
	#---cheap way to show a domain by zeroing it out
	if 'domain' in kwargs: 
		errormap[np.where(kwargs['domain'])] = 0.0
	cs = ax.contourf(X,Y,errormap,levels=levels,vmax=error_max,vmin=error_min,
		extend='both',origin='lower',lw=2,zorder=3,cmap=mpl.cm.jet)
	cs.cmap.set_over('w')
	if under_color: cs.cmap.set_under(under_color)
	cs_lines = ax.contour(X,Y,errormap,vmax=error_max,
		vmin=error_min,levels=levels[::contour_lineskip],
		extend='both',origin='lower',linewidths=0.5,colors='k',zorder=4)
	this_cut_extent = 0,min([cut_extent[1],max(extents)])
	ax.set_aspect((cut_curve[1]-cut_curve[0])/(this_cut_extent[1]-this_cut_extent[0]))
	ax.set_title(titlemaker(sn,fallback),fontsize=fs+2)
	#---take a subset to identify the well that lies in frame
	finex,finey = np.linspace(c0,c1,contour_interp_pts),np.linspace(e0,e1,contour_interp_pts)
	xsub = np.where(np.all(np.array((finex>=cut_curve[0],finex<=cut_curve[1])),axis=0))[0]
	ysub = np.where(np.all(np.array((finey>=cut_extent[0],finey<=cut_extent[1])),axis=0))[0]
	Xsub,Ysub = np.meshgrid(finex[xsub],finey[ysub])
	errormap_sub = errormap[ysub][:,xsub].T
	#---the following "best" method is deprecated and instead we look it up
	if old_show_best:
		best = [(Ysub[i,j],Xsub[i,j]) for i,j in 
			[np.unravel_index(errormap_sub.T.argmin(),errormap_sub.T.shape)]][0][::-1]
	#---the best hypothesis is pulled from packed, determined by weight_method set elsewhere
	elif not old_show_best and best_coarse: best = memory[pkey+('packed',)]['best_coarse']
	elif not old_show_best and not best_coarse: best = memory[pkey+('packed',)]['best_fine']
	else: raise Exception('unclear dot-placing method in plot_errormap')
	if show_best: 
		ax.scatter([best[0]],[best[1]],zorder=6,s=70,lw=2,color='w',facecolor='k',clip_on=False)
	ax.set_xlabel(r'$\mathrm{C_{0,max}\,({nm}^{-1})}$',fontsize=fs)
	ax.set_ylabel(r'$\mathrm{\boldsymbol{\sigma_{a}}\,(nm)}$',fontsize=fs)
	#---remove 1 because it gets crowded at the bottom
	ax.set_xticks([e for e in curvatures if e!=0.018])
	#---remove some values because it's crowded
	if False: ax.set_yticks([e for e in extents if e!=1.0])
	ax.set_xlim(cut_curve)
	ax.set_ylim(0,min([cut_extent[1],max(extents)]))
	ax.set_xticklabels(ax.get_xticks(),rotation=-90,fontsize=fs-2)
	ax.set_yticklabels(ax.get_yticks(),fontsize=fs-2)
	return cs

def color_dictionary(comparison,descriptor):
	"""
	Master listing of colors for different comparisons.
	"""
	import brewer2mpl
	#---name the colors from brewer here, used as the default set for now
	clrs = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors
	clrs_names = {'red':0,'blue':1,'green':2,'purple':3,'orange':4,'yellow':5,'brown':6,'pink':7,'grey':8,}
	
	colordict = {'tension':{'high':'red','no':'blue',},
		'enth':{'membrane-v550':'grey','membrane-v653-free':'grey','membrane-v652-enthx1':'orange',
		'membrane-v650-enthx4-dev':'blue','membrane-v614-enthx4-12800':'blue',
		'membrane-v616-octamer-close':'green','membrane-v651-enthx8':'green',},
		'exo70':{'membrane-v701-exo70-anti-dilute':'red','membrane-v700-exo70-dilute':'blue',},}

	#---use more washed out colors
	clrs = brewer2mpl.get_map('Set2','qualitative',8).mpl_colors
	clrs_names = {'cyan':0,'blue':2,'magenta':3,'gray':7,}
	colordict['enth']['membrane-v653-free'] = 'gray'
	colordict['enth']['membrane-v652-enthx1'] = 'cyan'
	colordict['enth']['membrane-v650-enthx4-dev'] = 'blue'
	colordict['enth']['membrane-v651-enthx8'] = 'magenta'
	return clrs[clrs_names[colordict[comparison][descriptor]]]

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---prepare postdat specific for the landscapes (as opposed to postdat_master)
if 'postdat' not in globals():
	kwargs = {}
	kwargs['contour_interp_pts'] = contour_interp_pts
	kwargs['n_best'] = n_best
	postdat_master = postdat = prepare_postdat(hypo_groups=hypo_groups,**kwargs)
	comparisons = postdat['comparisons']

#---plot a review panel
if 'summary' in routine:

	#---plot settings
	fsbase = 14
	#---loop over top-level curvature field comparison: protein/single, static/dynamic
	for group_key,tweak_cf_spec in enumerate(tweak_cf):
		subcomp = comparisons[group_key]
		test_sig = '_'.join([tweak_cf_spec[k] for k in ['motion','mapping']])
		nrows,ncols = 7,len(subcomp)
		fig_subdim = 4,3
		figlayout = {'out':{'grid':[nrows,1],'hspace':0.6},
			'ins':[{'grid':[1,ncols],'hspace':0.35,'wspace':0.35} for i in range(nrows)]}
		for j in range(3): figlayout['ins'][-1-j-1].update(wspace=0.5)
		axes,fig = panelplot(figsize=(fig_subdim[0]*ncols,fig_subdim[1]*nrows),layout=figlayout)
		#---plot the top row of raw error maps
		error_min = min([postdat[(m,s,f)]['errormap_raw'].min() for snum,(m,s,f) in enumerate(subcomp)])
		error_max_abs = max([postdat[(m,s,f)]['errormap_raw'].max() for snum,(m,s,f) in enumerate(subcomp)])
		error_min_abs = min([postdat[(m,s,f)]['errormap_raw'].min() for snum,(m,s,f) in enumerate(subcomp)])
		error_max = error_max_abs
		levels_abs = np.linspace(error_min,error_max,contour_nlevels)
		for snum,(mapping,sn,fallback) in enumerate(subcomp):
			pkey = (mapping,sn,fallback)
			ax = axes[0][snum]
			errormap_raw = postdat[pkey]['errormap_raw']
			curvatures,extents = [postdat[pkey][key] for key in ['curvatures','extents']]
			im = ax.imshow(errormap_raw.T,interpolation='nearest',origin='lower',
				cmap=mpl.cm.jet,extent=(0,len(curvatures),0,len(extents)),
				vmin=error_min,vmax=error_max_abs)
			ax.set_xticks(np.arange(len(curvatures))+0.5)
			ax.set_xticklabels(curvatures,rotation=-90,fontsize=fsbase-4)
			ax.set_yticks(np.arange(len(extents))+0.5)
			ax.set_yticklabels(extents,fontsize=fsbase-4)
			ax.set_title(titlemaker(sn,fallback))
			tickoff(ax,grid=False)
			best_abs = np.unravel_index(postdat[pkey]['errormap_raw'].argmin(),
				postdat[pkey]['errormap_raw'].shape)
			ax.scatter([best_abs[0]+0.5],[best_abs[1]+0.5],zorder=6,s=70,lw=2,facecolor='k',color='w')
			ax.set_xlim(0,len(curvatures))
			ax.set_ylim(0,len(extents))
		#---colorbar for the image plots
		ax = axes[0][-1]
		axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.05,0.,1.,1.),
			bbox_transform=ax.transAxes,borderpad=0)
		#---plot the contour maps
		error_max = error_max_abs if not error_ceiling else error_ceiling
		error_min = error_min_abs if not error_floor else error_floor
		level_ticks = np.linspace(error_min,error_max,cbar_nlevels)
		powlev = int(('%.2e'%level_ticks[-1]).split('e')[1])
		cbar = plt.colorbar(im,cax=axins,orientation="vertical",ticks=level_ticks) 
		cbar.ax.set_yticklabels(['%.1f'%(i/10**powlev) for i in level_ticks])
		axins.set_ylabel(r'MSE '+r'$\times 10^{%d}$'%(-1*powlev),rotation=270,labelpad=20)
		if error_max < error_min: 
			print '[WARNING] error_ceiling<error_min'
			error_min = error_ceiling
			error_max = error_max_abs
		levels = np.linspace(error_min,error_max,contour_nlevels)
		best_results = {}
		for snum,(mapping,sn,fallback) in enumerate(subcomp):
			pkey = (mapping,sn,fallback)
			ax = axes[1][snum]
			errormap = postdat[pkey]['errormap']
			curvatures,extents = [postdat[pkey][key] for key in ['curvatures','extents']]
			c0,c1 = min(curvatures),max(curvatures)
			e0,e1 = min(extents),max(extents)
			finex,finey = np.linspace(c0,c1,contour_interp_pts),np.linspace(e0,e1,contour_interp_pts)
			X,Y = np.meshgrid(finex,finey)
			cs = ax.contourf(X,Y,errormap,levels=levels,vmax=error_max,vmin=error_min,
				extend='both',origin='lower',lw=2,zorder=3,cmap=mpl.cm.jet)
			cs.cmap.set_over('w')
			cs_lines = ax.contour(X,Y,errormap,vmax=error_max,
				vmin=error_min,levels=levels[::contour_lineskip],
				extend='both',origin='lower',linewidths=0.5,colors='k',zorder=4)
			this_cut_extent = 0,min([cut_extent[1],max(extents)])
			ax.set_aspect((cut_curve[1]-cut_curve[0])/(this_cut_extent[1]-this_cut_extent[0]))
			ax.set_xlim(cut_curve)
			ax.set_ylim(0,min([cut_extent[1],max(extents)]))
			ax.set_xticklabels(ax.get_xticks(),rotation=-90)
			ax.set_title(titlemaker(sn,fallback))
			#---take a subset to identify the well that lies in frame
			finex,finey = np.linspace(c0,c1,contour_interp_pts),np.linspace(e0,e1,contour_interp_pts)
			xsub = np.where(np.all(np.array((finex>=cut_curve[0],finex<=cut_curve[1])),axis=0))[0]
			ysub = np.where(np.all(np.array((finey>=cut_extent[0],finey<=cut_extent[1])),axis=0))[0]
			Xsub,Ysub = np.meshgrid(finex[xsub],finey[ysub])
			errormap_sub = errormap[ysub][:,xsub].T
			best = [(Ysub[i,j],Xsub[i,j]) for i,j in 
				[np.unravel_index(errormap_sub.T.argmin(),errormap_sub.T.shape)]][0][::-1]
			ax.scatter([best[0]],[best[1]],zorder=6,s=70,lw=2,color='w',facecolor='k')
			best_results[pkey] = best
		#---colorbar for the contour plots
		ax = axes[1][-1]
		axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.05,0.,1.,1.),
			bbox_transform=ax.transAxes,borderpad=0)
		level_ticks = np.linspace(error_min,error_max,cbar_nlevels)
		powlev = int(('%.2e'%level_ticks[-1]).split('e')[1])
		cbar = plt.colorbar(cs,cax=axins,orientation="vertical",ticks=level_ticks) 
		cbar.ax.set_yticklabels(['%.1f'%(i/10**powlev) for i in level_ticks])
		axins.set_ylabel(r'MSE '+r'$\times 10^{%d}$'%(-1*powlev),rotation=270,labelpad=20)
		#---plot best spectra
		for snum,(mapping,sn,fallback) in enumerate(subcomp):
			pkey = (mapping,sn,fallback)
			ax = axes[2][snum]
			for cnum in range(n_best)[::-1]:
				wv,energy,error = postdat[(mapping,sn,fallback,'choice_cuts',cnum)]
				energy = np.abs(energy)
				wv_red = couplecalc.perfect_collapser(wv,wv)
				energy_red = couplecalc.perfect_collapser(wv,energy)
				kwargs = dict(lw=2.0,alpha=0.02,zorder=2+n_best-cnum,mew=0,ms=5,
					color=mpl.cm.jet((error-error_min)/(error_max-error_min)))
				if cnum==0: kwargs.update(path_effects=
					[matplotlib.patheffects.withStroke(linewidth=4,foreground='w')])
				ax.plot(wv,energy,'o',**kwargs)
				ln, = ax.plot(wv_red,energy_red,'-',lw=2.0,alpha=1.0,zorder=n_best+1+2+cnum,
					color=mpl.cm.jet((error-error_min)/(error_max-error_min)))
				ln.set_solid_capstyle('round')
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xlim(0.05,2*tweak['hicut'] if False else 10**1)
			ax.axvline(tweak['hicut'],c='k',lw=2)
			tickoff(ax,alpha=0.5)
			ax.set_ylabel('$\mathrm{%s}$'%'energy\,(k_B T)')
			ax.set_title(titlemaker(sn,fallback))
		#---scatter parameters vs error
		labels = {'kappa':r'\boldsymbol{\kappa} \,(k_B T)',
			'gamma':r'\boldsymbol{\gamma} \,(k_B T {nm}^{-2})',
			'vibe':'C'}
		maxima = {'kappa':100,'gamma':10,'vibe':10}
		param_lims = [(10**15,0) for j in range(3)]
		material_props = ['kappa','gamma','vibe']
		if 'kappa_gamma_prefactor' in tweak: prefactor = tweak['kappa_gamma_prefactor']
		else: prefactor = dict([(k,1.0) for k in material_props])
		for pnum,param in enumerate(material_props):
			for snum,(mapping,sn,fallback) in enumerate(subcomp):
				pkey = (mapping,sn,fallback)
				ax = axes[3+pnum][snum]
				survey = postdat[pkey]['survey']
				keys,params,errors = zip(*[(key,val[param],val['error']) for key,val in survey.items()])
				moderate = np.where(np.absolute(params)<maxima[param])[0]
				params,errors = [np.array(i)[moderate] for i in [params,errors]]
				params = np.array(params)*prefactor[param]
				subsel = np.where(np.all([np.array(errors)<=error_max,np.array(errors)>=error_min],axis=0))
				errors,params = np.array(errors)[subsel],np.array(params)[subsel]
				if len(params)==0:
					print('filter is too strong on params')
					import pdb;pdb.set_trace()
				param_lims[pnum] = [min([param_lims[pnum][0],min(params)]),
					max([param_lims[pnum][1],max(params)])]
				ax.set_ylabel('$\mathrm{%s}$'%labels[param])
				ax.set_title(titlemaker(sn,fallback))
				ax.scatter(errors,params,lw=0,zorder=3,
					color=[mpl.cm.jet((e-error_min)/(error_max-error_min)) for e in errors])
				ax.scatter(errors[errors.argmin()],params[errors.argmin()],zorder=4,lw=2,s=100,color='k',
					facecolor=[mpl.cm.jet((errors.min()-error_min)/(error_max-error_min))])
				ax.set_xticks(level_ticks)
				ax.set_xticklabels(['%.1f'%(i*10**(-1*powlev)) for i in level_ticks])
				ax.set_xlabel(r'MSE '+r'$\times 10^{%d}$'%(-1*powlev))
				ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(nbins=4,prune='both'))
				tickoff(ax)
		widen = lambda x,w=0.1 : [(x[0]-abs(j),x[1]+abs(j)) for j in [(x[1]-x[0])*w]][0]
		param_lims = [widen(p) for p in param_lims]
		wide_error = widen((error_min,error_max))
		for ii,i in enumerate(subcomp): 
			[axes[3+j][ii].set_ylim(param_lims[j]) for j in range(3)]
			[axes[3+j][ii].set_xlim(wide_error) for j in range(3)]
		#---plot the curvature fields
		for snum,(mapping,sn,fallback) in enumerate(subcomp):
			pkey = (mapping,sn,fallback)
			best = best_results[pkey]
			best_curvature = curvatures[np.abs(best[0]-curvatures).argmin()]
			best_extent = extents[np.abs(best[1]-extents).argmin()]
			lookup = {'mapping':tweak_cf_spec['mapping'],'motion':tweak_cf_spec['motion'],'isotropy':1.0,
				'sigma_a':best_extent,'sn':sn,'fallback':fallback,'curvature':1.0}
			pk = sessions['field'].query(Field).filter_by(**Field(**(lookup)).present()).one().id
			key_cf = ('curvature',pk)
			if key_cf not in memory: 
				fn_cf = namer_cf(pk)
				memory[key_cf] = load(os.path.basename(fn_cf),rootdir_cf)
			ax = axes[6][snum]
			max_curvature = max(zip(*best_results.values())[0])
			ax.imshow(best_curvature*memory[key_cf]['fields'].mean(axis=0).T,
				interpolation='nearest',origin='lower',
				vmin=-max_curvature,vmax=max_curvature,cmap=mpl.cm.RdBu_r,
				extent=next([0,i[0],0,i[1]] for i in [memory[(sn,'vecs')].mean(axis=0)]))
		#---print
		picturesave('fig.couple.contours_%s-%s_%s'%(roundsig,hcv,test_sig),
			work.plotdir,backup=False,version=True,meta={})
		plt.close()
