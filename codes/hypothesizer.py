#!/usr/bin/env python

def hypothesizer(*sweep,**kwargs):
	"""
	Code for sweeping an arbitrarily deep dictionary over many dimensions in combinations.
	"""
	from copy import deepcopy
	default = kwargs.get('default',None)
	def delve(o,*k): return delve(o[k[0]],*k[1:]) if len(k)>1 else o[k[0]]

	#---create the default hypothesis
	if default==None: default = {}
	for pathway in sweep:
		if pathway['route'][0] not in default: default[pathway['route'][0]] = {}
		for i in range(1,len(pathway['route'])-1):
			level = delve(default,*pathway['route'][:i])
			if pathway['route'][i] not in level: level[pathway['route'][i]] = {}
		if len(pathway['route'])>1: 
			delve(default,*pathway['route'][:-1])[pathway['route'][-1]] = pathway['values'][0]
	for i in default: default[i] = None if default[i] == {} else default[i]
	#---values to sweep over
	t = [i['values'] for i in sweep]
	allcombos = list([[i] for i in t[0]])
	for s in t[1:]:
		for bi in range(len(allcombos)):
			b = allcombos.pop(0)
			for r in list(s): allcombos.append(b + [r])
	#---accumulate hypotheses by traverse/replace
	hypotheses = []
	for combo in allcombos:
		newhypo = deepcopy(default)
		for routenum in range(len(sweep)):
			tmp = newhypo[sweep[routenum]['route'][0]]
			if type(newhypo[sweep[routenum]['route'][0]]) != dict:
				newhypo[sweep[routenum]['route'][0]] = combo[routenum]
			else:
				for i in sweep[routenum]['route'][1:-1]: tmp = tmp[i]
				tmp[sweep[routenum]['route'][-1]] = combo[routenum]
		hypotheses.append(newhypo)	
	return hypotheses
