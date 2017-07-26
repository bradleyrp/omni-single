#!/usr/bin/env python

"""
Retrieve an ITP file for a simulation.

Note that many calculations may require connectivity information from a TPR or an ITP. The ITP files are 
not currently tracked by the treeparser. You can find them systematically by retrieving GRO files. To do this,
use your YAML file to point to this function, which can replace a direct looking in the YAML.

Add the following to the `vars` dictionary in your YAML file to use this feature. 
Change the paths where appropriate.

  protein_itp_loader: 
    module: "codes.protein_itp_loader"
    function: "protein_itp_loader"
"""

import os,glob

def protein_itp_loader(sn,**kwargs):
	work = kwargs.get('work')
	last_start = work.raw.get_last_structure(sn)
	fns = glob.glob(os.path.join(os.path.dirname(last_start),'*.itp'))
	return sorted(fns,key=lambda x:os.path.getmtime(x))[-1]
