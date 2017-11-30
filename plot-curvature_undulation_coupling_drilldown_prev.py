#!/usr/bin/env python

"""
DRILLDOWN to investigate curvature coupling results
"""

from codes.undulate import calculate_undulations
from codes.undulate_plot import undulation_panel,add_undulation_labels,add_axgrid,add_std_legend
from codes.curvature_coupling.curvature_coupling_plots import individual_reviews_plotter

@autoload(plotrun)
def load():
	"""Load the data."""
	#---condition for loading the data
	if 'data' not in globals(): 
		#---load once
		global data_down,data,calc,calc_up
		data_down,calc = plotload('curvature_undulation_coupling_drilldown',work)
		data,calc_up = plotload('curvature_undulation_coupling',work)

@autoplot(plotrun)
def drilldown():
	"""
	"""
	pass

#---! iterative development
if 'data' in globals() and 'sns' not in globals():
	sns = work.sns()
	sn = sns[0]
	dat = data_down[sn]['data']

	#---! pick up a new dat file?
	base_fn = 'drilldown.curvature_coupling_s04'
	design_loop_tag = 'v14'

	#---compute then plot
	if not os.path.isfile(os.path.join(work.paths['post_data_spot'],'%s.dat'%base_fn)):
		#---unpacking calc specs from an unusual location
		calc_specs = calc['calcs']['specs']['specs']
		#---! softcode the loop position in the drilldown specs. !!! should this match the dat file
		design = calc_specs['design']['loop'][design_loop_tag]
		fitting = calc_specs.get('fitting',{})
		remember = {
			(sn,'drop_gaussians_points'):dat['drop_gaussians_points'],
			(sn,'fields_unity'):dat['fields_unity']}
		from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
		ic = InvestigateCurvature(sn=sn,work=work,
			design=design,fitting=fitting,remember=remember,
			protein_abstractor=data['protein_abstractor'][sn]['data'],
			undulations=data['undulations'][sn]['data'])
		import copy
		from omnicalc import store
		attrs,result = copy.deepcopy(ic.finding['attrs']),copy.deepcopy(ic.finding['result'])
		try: attrs['bundle'][sn].pop('final_simplex')
		except: pass
		store(obj=result,name='%s.dat'%base_fn,path=work.paths['post_data_spot'],
			attrs=attrs,verbose=True)
	import h5py
	fn  = os.path.join(work.paths['post_data_spot'],'%s.dat'%base_fn)
	fl = h5py.File(fn,'r')
	dat_incoming = dict([(key,np.array(val)) for key,val in fl.items()])
	dat.update(**dat_incoming)
	fl.close()
	#---note no subsampling in the views
	#---! height profile only has a dot for eight proteins right now
	plotspec = {'coupling_review_new':{'viewnames':[
		'average_height','spectrum']},}
	postdat = {sn:data['protein_abstractor'][sn]['data']}
	postdat[sn].update(points_protein=data_down[sn]['data']['drop_gaussians_points'])
	undulations_name,protein_abstractor_name = 'undulations','protein_abstractor'
	datas = dict(placeholder_tag={sn:dat})
	datas['placeholder_tag'][sn]['points_protein'] = data['protein_abstractor'][sn]['data']['points_protein']
	seep = dict(undulations_name=undulations_name,protein_abstractor_name=protein_abstractor_name,
		data=data,datas=datas,postdat=postdat,
		calcs=dict(placeholder_tag={sn:dict(calcs=calc['calcs'])}))
	for out_fn,details in plotspec.items(): individual_reviews_plotter(out_fn=out_fn,seep=seep,**details)
	# rm -f /home/rpb/sigma/analyze-project-curvature/post/drilldown.curvature_coupling_s02.dat && make plot curvature_undulation_coupling_drilldown

	#---from the previous drilldown
	if False:
		#---assemble a mapping between each element of the design loop and the set of optimizations we want
		#---! currently deterministic but we should make this extensible later
		#---! unfortunately we are recapping a lot of the compute funcitonality from omnicalc here (wasteful)
		table = []
		for tag in ['v13']:#datas:
			for sn in work.sns():
				for hnum,hypo in enumerate(hypos):
					status('computing %s,%s with the following hypothesis'%(tag,sn),tag='compute')
					asciitree(hypo)
					if False: ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M')
					#---add design back to the hypothesis
					hypo_more = dict(design=calcs[tag][sn]['calcs']['specs']['design'])
					hypo_more['design'].update(**hypo['design'])
					hypo_more['fitting'] = hypo['fitting']
					table.append(dict(hypo=hypo_more,sn=sn,tag=tag,num=hnum))
		#---run all hypotheses
		for item in table:
			tag = item['tag']
			sn = item['sn']
			num = item['num']
			hypo = item['hypo']
			fn = 'curvature_coupling_drilldown.%s.%s.%d'%(tag,work.namer.short_namer(sn,spot='sims'),num)
			if not os.path.isfile(os.path.join(work.paths['post_data_spot'],'%s.dat'%fn)):
				status('computing and saving %s'%fn,tag='compute')
				try: 
					dat = datas[tag][sn]
					remember = {
						(sn,'drop_gaussians_points'):dat['drop_gaussians_points'],
						(sn,'sampling'):dat['sampling'],
						(sn,'fields_unity'):dat['fields_unity']}
					ic = InvestigateCurvature(sn=sn,work=work,
						remember=remember,
						protein_abstractor=data['protein_abstractor'][sn]['data'],
						undulations=data['undulations'][sn]['data'],
						**hypo)
					#---package and save by timestamp for now
					attrs,result = copy.deepcopy(ic.finding['attrs']),copy.deepcopy(ic.finding['result'])
					#---final_simplex is not saveable
					try: attrs['bundle'][sn].pop('final_simplex')
					except: pass
					store(obj=result,name='%s.dat'%fn,path=work.paths['post_data_spot'],
						attrs=attrs,verbose=True)
				except: status('failed %s'%fn,tag='fail')

		#---collect everything
		if 'collected' not in globals():
			collected = []
			for item in table:
				tag = item['tag']
				sn = item['sn']
				num = item['num']
				hypo = item['hypo']
				fn = 'curvature_coupling_drilldown.%s.%s.%d'%(tag,work.namer.short_namer(sn,spot='sims'),num)
				fl = h5py.File(os.path.join(work.paths['post_data_spot'],'%s.dat'%fn),'r')
				dat_incoming = dict([(key,np.array(val)) for key,val in fl.items()])
				fl.close()
				collected.append(dat_incoming)

		figsize = (10,10)
		collection_inds = [ii for ii,i in enumerate(table) if i['sn']=='membrane-v1005']

		axes,fig = square_tiles(1,figsize)
		colors = [mpl.cm.__dict__['jet'](i) for i in np.linspace(0,1.,len(collection_inds))]
		
		for inum,cnum in enumerate(collection_inds):
			item = collected[cnum]
			high_cutoff = table[cnum]['hypo']['fitting']['high_cutoff']
			color = colors[inum]
			ax = axes[0]
			#for jj in range(len(item['qs'])):
			ax.scatter(item['qs'],item['ratios'],marker='.',s=2,color=color)
			ax.axvline(high_cutoff,color=color)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.axhline(1.0,c='k')
		picturesave('fig.WHAT',work.plotdir,backup=False,version=True,meta={})
