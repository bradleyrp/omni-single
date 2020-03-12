#!/usr/bin/env python

@autoload(plotrun)
def load():
	"""Load everything for the plot only once."""
	data = plotload(plotname)
	sn = 'membrane-v566'
	dat = data.this[sn]

if __name__=='__main__':
	if 0:
		mn = 0
		nframes = dat['nframes']
		resnames = dat['resnames']
		imono = dat['monolayer_indices']
		simplices = [dat['%d.%d.simplices'%(mn,fr)] for fr in range(nframes)]
		ghosts = [dat['%d.%d.ghost_ids'%(mn,fr)] for fr in range(nframes)]
		subjects = np.where(resnames[np.where(imono==mn)[0]]=='PI2P')[0]
		counts_triplets,triplets_this = [],[]
		for fr in range(nframes):
			print(fr)
			s = simplices[fr]
			g = ghosts[fr]
			# subjects ghost for this frame
			sg = np.where(np.in1d(g,subjects))[0]
			# framewise count of the lipids in each simplex which are PIP2
			counts = np.sum(np.array([np.in1d(i,sg) for i in s]),axis=1)
			# count the triplets
			counts_triplets.append(np.sum(np.sum(np.array([np.in1d(i,sg) for i in s]),axis=1)==3))
			# save the triplets at each frame
			triplets_this.append(s[np.sum(np.array([np.in1d(i,sg) for i in s]),axis=1)==3])
	# get unique triplets
	tu = []
	for fr,this in enumerate(triplets_this):
		for i in this:
			if tuple(i) not in tu: tu.append(tuple(i))
	whiches = np.array([[i in t for t in triplets_this] for i in tu])
	if 1:
		ax = plt.subplot(111)
		#for i in range(len(whiches))
		#	ax.plot(w)
		for i in range(len(whiches)):
			ax.plot(
				i+whiches[i],c=['r','g','b','k'][i],alpha=0.5
			)
		plt.savefig(work.plotdir+'TMP.jpg')
		plt.close()