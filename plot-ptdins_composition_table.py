#!/usr/bin/env python

"""
"""

# autoplot settings
plotrun.routine = None

@autoload(plotrun)
def loader():
	"""Load data."""
	data = plotload(plotname)
	sns = work.sns()

if __name__=='__main__':

	if 'toc' not in globals():

		gro_files = dict([(sn,
			os.path.join(work.postdir,data[1]['extras'][sn]['slice_path']+'.gro'))
			for sn in sns])

		import makeface
		try: mod = makeface.import_remote('amx/amx')
		except: raise Exception('please clone a copy of automacs next to omni in `amx`')
		GMXStructure = mod['GMXStructure']

		toc = {}
		for sn,fn in gro_files.items():
			print('reading %s'%sn)
			struct = GMXStructure(fn)
			#! minor mods to amx for this
			toc[sn] = struct.detect_composition()

	do_latex_table = True

	# get all residue names
	resnames_this = sorted(list(set([i[0] for j in toc.values() for i in j])))
	# override this ordering
	resnames = 'POPC DOPE DOPC DOPS CHL1 PI2P P35P PIPP PIPU SAPI K NA MG Cal CL SOL'.split()
	if set(resnames_this)!=set(resnames): raise Exception
	['CHL1', 'CL', 'Cal', 'DOPC', 'DOPE', 'DOPS', 'K', 'MG', 'NA', 'P35P', 'PI2P', 'PIPP', 'PIPU', 'POPC', 'SAPI', 'SOL']
	counts = np.zeros((len(resnames),len(sns))).astype(int)
	for snum,sn in enumerate(sns):
		for i,j in dict(toc[sn]).items():
			counts[resnames.index(i)][snum] = int(j)
	if do_latex_table:
		print(r'\begin{tabular{ l | c | r } \hline')
	else:
		#! you're on your own for the titles
		print(' '*20+''.join([('%s'%r).rjust(10) for r in resnames]))
	for rnum,row in enumerate(counts.T):
		if not do_latex_table:
			dat = ''.join([('%d'%i).rjust(10) for i in row])
			print(('%s'%sns[rnum]).rjust(20)+dat)
		else:
			dat = ' & '.join([str(i) for i in row])+' \\\\'
			print(dat)
	if do_latex_table:
		print(r'\hline \end{tabular}')
