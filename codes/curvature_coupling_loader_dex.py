#!/usr/bin/env python

"""
???
"""

def prepare_mesoscale_postdat(data):
	"""
	"""
	#---hard-coding the keys
	dat_prot = data['import_readymade_meso_v1_nanogel']
	dat_memb = data['import_readymade_meso_v1_membrane']
	sns = dat_prot.keys()
	postdat = dict([(sn,{}) for sn in sns])
	for sn in sns:
		postdat[sn]['vecs'] = dat_memb[sn]['data']['vecs'].mean(axis=0)
		postdat[sn]['points_protein'] = dat_prot[sn]['data']['points_all']
		postdat[sn]['points_protein_mean'] = postdat[sn]['points_protein'].mean(axis=1).mean(axis=0)
	return postdat
