#!/usr/bin/python

from tools import *

commons = """
pbc: package require pbc
bonder:|
	source CGBONDSPATH
	cg_bonds -cutoff BOND_CUTOFF_ANGSTROMS -gmx GMXDUMP -tpr TPR
load: mol new GRO
load dynamic:|
	mol new GRO
	animate delete beg 0 end 0 skip 0 0
	mol addfile XTC FRAMES step STEP waitfor all 
ortho:
	display projection Orthographic
standard:|
	display projection Orthographic
	color Display Background BACKCOLOR
	axes location off
	display resize VIEWX VIEWY
	mol delrep 0 top
dev:|
	mol selection "SELECT"
	mol addrep top
play: animate forward
reset:|
	display resetview
	after 3000
xview:|
	mouse stoprotation
	rotate x to -90
	rotate y by -90
yview:|
	rotate z to 180
	rotate x by -90
zview:|
	mouse stoprotation
	rotate z to 180
backbonename: name CA
align_backbone:|
<<<<<<< HEAD
	set reference [atomselect top "BACKBONE_SELECT" frame 0]
	set compare [atomselect top "BACKBONE_SELECT"]
=======
	set reference [atomselect top "BACKBONENAME" frame 0]
	set compare [atomselect top "BACKBONENAME"]
>>>>>>> c3db4a14cbba03e94a15950b51f09c06acff8b9a
	set num_steps [molinfo top get numframes]
	for {set frame 0} {$frame < $num_steps} {incr frame} {
		$compare frame $frame
		set trans_mat [measure fit $compare $reference]
		$compare move $trans_mat
	}
snapshot:|
	set filename SNAPSHOT_FILENAME
	render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" \\
	-aasamples 12 %s -format TARGA -o %s.tga -res XRES YRES
	exec convert $filename.tga $filename.png
	exec rm $filename.tga $filename
draw_arrow:|
	proc vmd_draw_arrow {mol start end} {
		set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
		graphics $mol cylinder $start $middle radius 0.15
		graphics $mol cone $middle $end radius 0.25
	}
"""

drawing = """
beads: Beads 1.000000 12.000000
vdw: VDW 0.600000 12.000000
licorice: Licorice 0.300000 12.000000 12.000000
licorice_coarse: Licorice 1.300000 12.000000 12.000000
cartoon: NewCartoon 0.300000 10.000000 4.100000 0
cpk_water: CPK 1.000000 0.300000 10.000000 10.000000
lines: Lines 4.000000
points: Points 20.000000
"""

extras = """
update: mol selupdate REP top 1
smooth: mol smoothrep top REP 5
rainbow: mol modcolor REP top ResID
resname_color: mol modcolor REP top ResName
structure_color: mol modcolor REP top Structure
beta_color: mol modcolor REP top Beta
transparent: mol modmaterial REP top Transparent
edgyglass: mol modmaterial REP top EdgyGlass
glass1: mol modmaterial REP top Glass1
goodsell: mol modmaterial REP top Goodsell
hide: mol showrep top REP 0
color_specific: mol modcolor REP top ColorID $color_cursor
glass: mol modmaterial REP top Diffuse
xy:|
	mol showperiodic top REP x
	mol numperiodic top REP 1
	mol showperiodic top REP xy
	mol numperiodic top REP 1
	mol showperiodic top REP xyX
	mol numperiodic top REP 1
	mol showperiodic top REP xyXY
	mol numperiodic top REP 1
pbcz:|
	mol showperiodic top REP z
	mol numperiodic top REP 1
	mol showperiodic top REP zZ
	mol numperiodic top REP 1
"""

defaults = """
res: (1000,1000)
viewbox: (800,800)
backcolor: white
bond_cutoff_angstroms: 20
step: 1
script_fn: script-vmd.tcl
vmd_prepend: "VMDNOCUDA=1"
backbone_select: "name CA"
"""

recipes = """
atomistic bilayer licorice:|
	{'lipids':'not water and not ions and not name Cal and not protein and not hydrogen',
	'smooth':True,'style':'licorice','goodsell':True}
atomistic bilayer lines:|
	{'lipids':'not water and not ions and not protein and not hydrogen',
	'smooth':True,'style':'lines'}
atomistic protein:|
	{'protein_subject':'protein','smooth':True,'style':'cartoon','structure_color':True,'goodsell':True}
atomistic all:|
	{'everything':'all','smooth':True,'style':'cartoon','structure_color':True,'goodsell':True}
waterdots:|
	{'water':'water and name OW','smooth':True,'style':'cpk_water','goodsell':True}
coarse bilayer lines:|
	{'lipids':'not resname W and not resname ION and '+
		'not name BB and not name SC1 and not name SC2 and not name SC3 and not name SC4',
	'smooth':True,'style':'lines'}
coarse bilayer licorice:|
	{'lipids':'not resname W and not resname ION and '+
		'not name BB and not name SC1 and not name SC2 and not name SC3 and not name SC4',
	'smooth':True,'style':'licorice_coarse','goodsell':True}
coarse protein:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB','SC1','SC2','SC3','SC4']]),
	'smooth':True,'style':'licorice_coarse','color_specific':True}
coarse protein backbone:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB']]),
	'smooth':True,'style':'licorice_coarse','color_specific':True,'goodsell':True}
"""

recipes_collect = """
video aamd atomistic bilayer protein:['atomistic_bilayer licorice','atomistic_protein']
live aamd atomistic bilayer protein:['atomistic_bilayer lines','atomistic_protein']
video cgmd bilayer protein:['coarse bilayer licorice','coarse protein']
video cgmd bilayer protein backbone:['coarse bilayer licorice','coarse protein backbone']
video aamd atomistic all waterdots:['atomistic_bilayer licorice','waterdots']
video aamd generic:['atomistic_bilayer licorice']
"""

import os,sys,time,re,subprocess,tempfile,glob

def remove_box_jumps(step,last_step,last_part):

	"""
	Convert the trajectory to reassemble broken molecules.
	"""

	#---run get_last_frame beforehand
	assert os.path.isfile(step+'system-input.tpr')
	xtc =last_step+'md.part%04d.xtc'%last_part
	bash(gmxpaths['trjconv']+' -f %s -s %s -o %s -pbc mol'%(
		os.path.join('../',xtc),
		'system-input.tpr',
		'md.part%04d.pbcmol.xtc'%last_part),
		log='trjconv-pbcmol',cwd=step,
		inpipe='0\n')
	bash(gmxpaths['trjconv']+
		' -f %s -s %s -o %s -pbc mol'%(
		'system-input.gro','system-input.tpr','system-input.pbcmol.gro'),
		log='trjconv-pbcmol-gro',cwd=step,inpipe='0\n')
	gro = 'system-input.gro'
	tpr = 'system-input.tpr'
	xtc = 'md.part%04d.pbcmol.xtc'%last_part
	return {'xtc':xtc,'gro':gro,'tpr':tpr}

class VMDWrap:

	#---hard-coded vmd colors
	vmd_colors = dict([(i,ii) for ii,i in enumerate(
		('blue red gray orange yellow tan silver green white pink cyan purple lime '+
		'mauve ochre iceblue black').split())])

	def __init__(self,name='vmd',**kwargs):

		"""
		Interface with VMD to make videos and execute common rendering functions.
		Need to ensure viewbox corresponds to resolution.
		"""

		global commons,defaults,extras,drawing
		self.commons = yamlparse(commons,style='tabbed')
		self.defaults = yamlparse(defaults,style='tabbed')
		self.extras = yamlparse(extras,style='tabbed')
		self.drawing = yamlparse(drawing,style='tabbed')
		self.recipes = yamlparse(recipes,style='tabbed')
		self.recipes_collect = yamlparse(recipes_collect,style='tabbed')

		#---add any tcl scripts in this folder to the commons
		for fn in glob.glob(os.path.join(os.path.dirname(__file__),'*.tcl')):
			tcl_lib = re.search('^(.+)\.tcl$',os.path.basename(fn)).group(1)
			assert tcl_lib not in self.commons,'already found "%s" in the commons'%tcl_lib
			self.commons[tcl_lib] = open(fn).read()

		#---load defaults
		self.__dict__.update(**self.defaults)
		#---all kwargs are saved to the dictionary as substitutions
		self.__dict__.update(**kwargs)
		self.start_time = time.time()
		self.viewx,self.viewy = self.viewbox
		self.xres,self.yres = kwargs.get('res',self.viewbox)

		#---handle paths where the 'site' is the top level folder
		self.cwd = kwargs.get('site',os.getcwd())

		#---counters
		self.script = []
		self.selections = {}
		self.selection_order = []
		self.set_color_cursor('black')

	def __setitem__(self,key,value):

		"""
		Use the class like a dictionary for updating substitutions in the commons.
		"""

		self.__dict__[key] = value


	def do(self,*names):

		"""
		Use a set of prescribed code blocks defined in commons.
		"""

		for name in names:
			chunk = self.commons[name].strip('\n')
			for key,val in self.__dict__.items(): 
				chunk = re.sub('(\s|\")(%s)(\s|$|\")'%key.upper(),r"\g<1>%s\g<3>"%str(val),chunk)
			self.script.append(chunk)

	def select(self,**kwargs):

		"""
		Create a new selection with string substitutions.
		"""

		style = kwargs.pop('style') if 'style' in kwargs else None
		extras = []
		for key,val in kwargs.items():
			if type(kwargs[key]) == bool:
				kwargs.pop(key)
				if val: extras.append(key)
		dynamics = kwargs.pop('dynamics') if 'dynamics' in kwargs else None
		for key in kwargs:
			selstring = str(kwargs[key])
			for sel in self.selections:
				selstring = re.sub(sel.upper(),'(%s)'%self.selections[sel],selstring)
			self.script += ["mol selection \"(%s)\""%selstring,"mol addrep top"]		
			self.script += ['set %s [atomselect top "%s"]'%(key,selstring)]
			self.selection_order.append(key)
			self.selections[key] = kwargs[key]
		if style != None:
			for key in kwargs:
				style = self.drawing[style] if style in self.drawing else style
				self.script += ['mol modstyle %d top %s'%(self.selection_order.index(key),style)]
		for extra in extras:
			if extra in self.extras and extra:
				for key in kwargs:
					instruct = str(self.extras[extra])
					instruct = re.sub('REP','%d'%self.selection_order.index(key),instruct)
					self.script += [instruct]
			else: raise Exception('missing extra setting: %s'%extra)			

	def set_color_cursor(self,color):

		"""
		The $color_cursor variable in tcl is set to black (16) by the standard package but can be changed
		to use coloring by a particular ColorID, on the fly.
		"""

		self.script.append('set color_cursor %d'%self.vmd_colors[color])

	def command(self,text): 

		"""
		Add a raw command to the script if necessary.
		"""

		self.script.append(text)

	def write(self):

		"""
		Write the script.
		"""

		with open(self.cwd+'/'+self.script_fn,'w') as fp: 
			for line in self.script: fp.write(line+'\n')

	def recipe(self,recipe,**kwargs):

		"""
		Add from a set of predefined selections.
		"""

		recipe = dict(self.recipes[recipe])
		recipe.update(**kwargs)
		self.select(**recipe)

	def show(self,**kwargs):

		"""
		Run VMD.
		The pipe option will feed commands directly to VMD for systems with a wonky setup.
		!option to remove everything when done (clean)?
		!dispdev text option (text)?
		"""

		quit = kwargs.get('quit',True)
		pipe = kwargs.get('pipe',False)
		#---text mode if quit on complete
		text = kwargs.get('text',quit)


		#---quit after the script is done
		if quit: self.script += ['quit']
		self.write()
		#---! arrest execution if you on light running dark: import ipdb;ipdb.set_trace()
		#---open VMD for the user		
		if quit == False:
			#---if we are not quitting then we run vmd directly
			#---note this might fail if the weird segfault that motivated "not run_script" happens
			os.system('cd %s && %s vmd %s -e %s'%(self.cwd,self.vmd_prepend,
				'-dispdev text' if text else '',self.script_fn))
		#---quit after running
		else: 
			#---feed the commands directly to VMD to preclude the weird segfault errors
			if not pipe:
				text = '\n'.join(self.script)
				proc = subprocess.Popen('vmd -dispdev text',stdin=subprocess.PIPE,
					shell=True,executable='/bin/bash',stdout=sys.stdout,stderr=sys.stdout,
					cwd=self.cwd)
				proc.communicate(input=text)
			else:
				#---call the VMD script (note that this fails sometimes on systems with a poor VMD setup)
				subprocess.check_call('%s vmd %s -e %s'%(self.vmd_prepend,'-dispdev text' 
					if text else '',self.script_fn),shell=True,cwd=self.cwd)
		return True

	def video(self,traces=None,snapshot=False,pause=False,cwd='',dn=None):

		"""
		Prepare to render a video by writing an add-on script.
		"""

		#---ensure that the viewbox and the resolution have the same proportions
		if float(self.viewbox[0])/self.viewbox[1]!=float(self.res[0])/self.res[1]:
			self.viewbox = tuple([int(round(float(max(self.viewbox))/max(self.res)*i)) for i in self.res])
			self.script.append('display resize %d %d'%self.viewbox)
		#---add a "trace" or tail to some elements
		if traces != None:
			trace_commands = []
			if type(traces) != list: raise Exception('traces must be a list or None')
			for key in traces:
				trace_commands.append('mol drawframes top %d $lag:$lead'%self.selection_order.index(key))
		else: trace_commands = []
		#---generic header for making a video
		video_lines = [
			'set num [molinfo top get numframes]',
			'for {set i 0} {$i < $num} {incr i} {',
			'animate goto $i',
			'set lag [expr $i-15]',
			'set lead [expr $i]']
		video_lines.extend(trace_commands)
		#---all videos have snapshots rendered in a temporary folder
		if not dn: self.video_dn = tempfile.mkdtemp(dir=self.cwd)
		else:
			self.video_dn = dn
			dropspot = os.path.join(self.cwd,self.video_dn)
			if not os.path.exists(dropspot): os.mkdir(dropspot)
			assert glob.glob(os.path.join(dropspot,'*')) == []
		dropspot = os.path.join(os.path.basename(self.video_dn),'')
		#---snapshot is the cheap/fast method otherwise we render
		self.snapshot = snapshot
		if not self.snapshot:
			video_lines.extend([
				'set filename '+dropspot+'snap.[format "%04d" $i]',
				'render TachyonInternal $filename.tga' if False else
				'render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" '+\
				'-aasamples 12 %s -format TARGA -o %s.tga -res '+'%d %d'%self.res,
				'exec convert $filename.tga $filename.png',
				'exec rm $filename.tga $filename',
				'}'])
		else:
			video_lines.extend([
				'set filename '+dropspot+'snap.[format "%04d" $i]',
				'render snapshot $filename.ppm',
				'}'])
		#---write script-video.tcl which will be sourced by script-vmd.tcl
		with open(self.cwd+'/script-video.tcl','w') as fp:
			for line in video_lines: fp.write(line+'\n')
		if not pause: self.script.append('source %sscript-video.tcl'%os.path.join(cwd,''))
		return self.video_dn

	def render(self,name='video',rate=1.0,codec=None,size=50.0,duration=0.0):

		"""
		Render a video from snapshots.
		"""

		def is_installed(command):
			return not re.search('command not found','\n'.join(
			subprocess.Popen(command,shell=True,executable='/bin/bash',
			stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()))

		encoders = ['ffmpeg','avconv']
		use_encoder = next((i for i in encoders if is_installed(i)),None)
		if not use_encoder: raise Exception('cannot find any encoders in %r'%encoders)

		#---when using ffmpeg confirm that the codec is available 
		if not codec and use_encoder=='ffmpeg':
			#---priority codecs
			ffmpeg_codec_priorities = [('mpeg4','mp4'),('mpeg2video','mpeg')]
			proc = subprocess.Popen('ffmpeg -codecs',shell=True,
				stdout=subprocess.PIPE,stdin=subprocess.PIPE)
			std,err = proc.communicate()
			ffmpeg_codecs = [j.split()[1] for j in std.splitlines()[next(ii for ii,i in 
				enumerate(std.splitlines()) if re.match('^\s*\-+\s*$',i))+1:]]
			codec = next((i for i in zip(*ffmpeg_codec_priorities)[0] if i in ffmpeg_codecs),None)
			if not codec: 
				raise Exception('could not find good codecs in ffmpeg from %r'%
					zip(*ffmpeg_codec_priorities)[0])
			video_suffix = dict(ffmpeg_codec_priorities)[codec]

		assert not (size and use_encoder!='ffmpeg')
		two_pass = '2pass' if size else '1pass'

		#---index commands by encoder,self.snapshot
		dropspot = os.path.join(os.path.basename(self.video_dn),'')
		#---bitrate for 2pass is currently somewhat conservative and undershoots the desired size
		nframes = float(len(glob.glob(os.path.join(self.cwd,self.video_dn,'*.png'))))
		frate = 24.0
		if not size: bitrate = None 
		elif duration==0.0: bitrate = float(size)*8192/(nframes/frate)
		#---duration overrides rates
		else: 
			bitrate = float(size)*8192/(float(duration))
			rate = float(duration)/(nframes/frate)
		commands = {
			('ffmpeg','snapshot','1pass'):r"ffmpeg -i "+dropspot+"snap.%04d.ppm "+"-vcodec %s -q 0"%codec+
				"-filter:v 'setpts=%.1f*PTS' "%rate+name+"."+video_suffix,
			('ffmpeg','tachyon','1pass'):r"ffmpeg -i "+dropspot+"/snap.%04d.png "+"-vcodec %s -q 0 "%codec+
				"-filter:v 'setpts=%.1f*PTS' "%rate+name+"."+video_suffix,
			#---codec options for avconv are currently hardcoded
			('avconv','snapshot','1pass'):r"avconv -r %.1f -i "%(rate*4)+dropspot+
				"/snap.%04d.png "+self.cwd+'/'+name+".mp4",
			('avconv','tachyon','1pass'):r"avconv -r %.1f -i "%(rate*4)+dropspot+
				"/snap.%04d.png "+self.cwd+'/'+name+".mp4",
			#---! currently the 2pass method is hard-coded
			('ffmpeg','tachyon','2pass'):r"ffmpeg -y -i "+dropspot+"snap.%04d.png "+
				"-framerate %d "%frate+"-filter:v 'setpts=%.1f*PTS' "%round(rate,3)+
				"-c:v libx264 -b:v %dk -pass 1 -f mp4 /dev/null && "%int(bitrate)+
				"ffmpeg -y -i "+dropspot+"snap.%04d.png -framerate "+"%d "%frate+
				"-c:v libx264 -preset medium -b:v %dk -pass 2 "%int(bitrate)+
				"-filter:v 'setpts=%.1f*PTS' "%round(rate,3)+name+".mp4",
			}
		command_key = (use_encoder,'snapshot' if self.snapshot else 'tachyon',two_pass)
		assert rate==1.0 or command_key==('ffmpeg','tachyon','2pass'),'rates must be 1.0 for 2pass'
		if command_key not in commands: 
			raise Exception('cannot parse the following command options for the encoder: "%r"'%command_key)
		subprocess.check_call(commands[command_key],cwd=self.cwd,shell=True)
		print('[NOTE] ffmpeg command was: "%s"'%commands[command_key])
		if command_key==('ffmpeg','tachyon','2pass'):
			print('[NOTE] cleaning up ffmpeg2pass files')
			for fn in glob.glob(self.cwd+'/ffmpeg2pass*'): os.remove(fn)
			
	def martini_secondary_structure(self,itp,style='all',name='protein_subject',back_color=None,mers=1):

		"""
		Color a MARTINI protein by secondary structure.
		"""

		import re
		#---set the jet colorscale from colorscales.tcl
		self.do('colorscales')
		self.script.append('colorscale_jet')
		#---backbone forces secondary structure coloring
		if style=='backbone': back_color = None
		#---select secondary structure colors
		coloring = {'C':'gray','H':'red','3':'pink','1':'pink','S':'blue','2':'pink','T':'yellow'}
		#---! HACKED
		coloring = {'C':0.5,'H':1.0,'3':1.0,'1':1.0,'S':0.0,'2':1.0,'T':0.5}
		#---hard-coded bead numbers generated from lib_helices.CoarseGrainedModel
		bead_counts = {'ILE':2,'GLN':2,'GLY':1,'PHE':4,'GLU':2,'CYS':2,'HIS':4,'SER':2,
			'LYS':3,'PRO':2,'HISH':4,'ASN':2,'VAL':2,'THR':2,'ASP':2,'ASP0':2,'LYS0':3,
			'ARG0':3,'GLU0':2,'ALA':1,'MET':2,'LEU':2,'ARG':3,'TYR':4,'TRP':5}
		aa_codes3 = ['TRP','TYR','PHE','HIS','ARG','LYS','CYS','ASP','GLU',
			'ILE','LEU','MET','ASN','PRO','HYP','GLN','SER','THR','VAL','ALA','GLY']
		aa_codes1 = 'WYFHRKCDEILMNPOQSTVAG'
		assert os.path.isfile(itp),'need an ITP but I cannot find "%s"'%itp
		with open(itp) as fp: text = fp.read()
		regex_ss = '^;\s*(S|s)econdary structure\n;\s*([A-Z0-9]+)\s'
		ss = re.search('^;\s*(?:S|s)econdary (?:S|s)tructure.*?\n;\s*(.*?)$',text,re.M+re.DOTALL).group(1)
		seq = re.search('^;\s*(?:S|s)equence.*?\n;\s*(.*?)$',text,re.M+re.DOTALL).group(1)
		mapping = dict([(i,float(ii)) for ii,i in enumerate(['C','H','3','1','S','2','T'])])
		obs_bead_counts = [bead_counts[dict(zip(aa_codes1,aa_codes3))[i]] for i in seq]
		#betas = [(self.vmd_colors[back_color] if i==0 and back_color else self.vmd_colors[coloring[s]]) 
		#	for ii,s in enumerate(ss) for i in range({'all':obs_bead_counts[ii],'backbone':1}[style])]
		betas = [coloring[s] for ii,s in enumerate(ss) 
			for i in range({'all':obs_bead_counts[ii],'backbone':1}[style])]
		#---if the ITP is a single monomer and there is more than one, we use the "mers" flag
		if mers>1: betas = betas*mers
		betas_string = "{%s}"%' '.join(['%.2f'%i for i in betas])
		self.script.append('set beta_%s %s'%(name,betas_string))
		self.script.append('$%s set beta $beta_%s'%(name,name))
		self.script.append('mol modcolor %d top Beta'%self.selection_order.index('protein_subject'))

def view_routine(**kwargs):

	"""
	Function which interprets incoming video requests for VMDWrap and then executes them.
	Note that this function prepares a local viewer script for each video so users can tune/rerender.
	"""

	#---distinguish incoming paths based on amx/omni calls
	caller = kwargs.pop('caller','amx')
	assert caller in ['omni','amx'],'caller must be omni or amx'
	#---if automacs calls, the wordspace is fully populated
	wordspace = kwargs.get('wordspace',None)
	#---if omnicalc calls, we populate the wordspace
	if not wordspace:
		if caller=='amx': raise Exception('must pass wordspace if running from amx')
		wordspace = WordSpace()
		wordspace['trajectory_details'] = kwargs.pop('trajectory_details',{})
	#---fold settings into the wordspace
	wordspace.update(**yamlparse(kwargs['settings']))
	snap_dn_abs = kwargs.pop('snap_dn_abs','snapshots')
	video_name = kwargs.pop('video_name',wordspace.get('video_name','video'))
	#---amx uses a local wordspace for iterative rerendering
	if caller=='amx': wordspace_out = wordspace.step+'/'+'wordspace.json'
	elif caller=='omni': 
		film_root = os.path.abspath(os.path.join(snap_dn_abs,'..'))
		wordspace_out = os.path.join(film_root,'%s.wordspace.json'%video_name)
		#---if the wordspace file exists we might be doing a rerender
		if os.path.isfile(wordspace_out):
			print('[NOTE] found wordspace for this video: "%s"'%wordspace_out)
			incoming_wordspace = json.load(open(wordspace_out))
			wordspace.update(**incoming_wordspace)
			#---current settings override the previous wordspace
			#---previous wordspace really just tracks the development state for rerendering
			wordspace.update(**yamlparse(kwargs['settings']))

	#---recapitulate script-viewer here with amx/omni switching
	try:
		if wordspace.view_mode=='video': wordspace.show_trajectory = True
		#---stage 1: create trajectory with unbroken molecules
		if not wordspace['under_development']:
			if caller=='amx':
				wordspace['last_step'],wordspace['last_part'] = detect_last()
				start(wordspace.step)
				#---always get the last frame for viewing
				get_last_frame(tpr=True,cpt=False,ndx=True,top=False,itp=False)
				wordspace.trajectory_details = remove_box_jumps(step=wordspace.step,
					last_step=wordspace.last_step,last_part=wordspace.last_part)
			else:
				if 'trajectory_details' not in wordspace:
					raise Exception('view_routine called with caller="amx" '+
						'but we need trajectory_details')
			wordspace.under_development = 1
			write_wordspace(wordspace,outfile=wordspace_out)
		#---stage 2: prepare the viewer scripts
		if caller=='amx': site = wordspace.step
		elif caller=='omni': site = film_root
		view = VMDWrap(site=site,
			gro='system-input.gro' if caller=='amx' else wordspace.trajectory_details['gro'],
			res=wordspace.resolution,viewbox=wordspace.viewbox)
		if wordspace['cursor_color']: view.set_color_cursor(wordspace.cursor_color)
		#---custom variables
		if wordspace['custom_variables']:
			assert type(wordspace.custom_variables)==dict
			view.__dict__.update(**wordspace.custom_variables)
		#---if dynamics we prepare a trajectory with whole molecules
		if wordspace.show_trajectory:
			#---send trajectory files names
			view.__dict__.update(**wordspace.trajectory_details)
			#---use MDAnalysis to get the right frame counts
			import MDAnalysis
			uni = MDAnalysis.Universe(
				(wordspace.step+'/' if caller=='amx' else '')+view.gro,
				(wordspace.step+'/' if caller=='amx' else '')+view.xtc)
			stepskip = max(int(float(len(uni.trajectory))/wordspace.max_frames+1),1)
			view.__dict__.update(step=stepskip,frames='')
		view.do('load_dynamic' if wordspace.show_trajectory else 'load','standard')
		if wordspace['backbone_align_name']:
			view.__dict__['backbonename'] = wordspace.backbone_align_name
			view.do('align_backbone')
		if wordspace['cg_bonds'] or wordspace.bonder:
			cg_bonds_tcl = os.path.abspath(os.path.expanduser(str(wordspace['cg_bonds'])))
			if not os.path.isfile(cg_bonds_tcl): 
				raise Exception('\n'+''.join(['[ERROR] %s\n'%i for i in 
					['cannot find cg_bonds.tcl','get this file from MARTINI',
					'and set the settings variable called "cg_bonds" to that path']]))
			view.__dict__['cgbondspath'] = cg_bonds_tcl
			#---! need a protocol for GMXDUMP possibly tempfile or alias
			view.__dict__['gmxdump'] = 'gmxdump'
		if wordspace.bonder: view.do('bonder')
		view.do(wordspace.which_view)
		if wordspace['scale']: view.command(wordspace.scale)
		for select in wordspace.selections: view.select(**select)
		if wordspace.recipe_collection:
			for recipe in view.recipes_collect[nospaces(wordspace.recipe_collection)]: 
				view.recipe(nospaces(recipe))
		#---apply special MARTINI protein colors
		if wordspace['martini_structure_color']: 
			if wordspace['recipe_collection']=='video cgmd bilayer protein backbone':
				martini_style = 'backbone'
			else: martini_style = 'all'
			view.martini_secondary_structure(itp=wordspace.itp,
				style=martini_style,back_color=wordspace.backbone_color)
			view.martini_secondary_structure(itp=wordspace.itp,style=martini_style,
				back_color=wordspace.backbone_color,mers=wordspace.get('mers',1))
		#---custom scripts directly from the wordspace
		if wordspace['custom_script']:
			script = wordspace.custom_script
			assert type(script)==list,'custom script must be a list'
			for line in script: eval(line)
		#---stage 3: make snapshots if this is a video
		if wordspace.view_mode == 'video' and not wordspace['video_snaps']: 
			wordspace.video_snaps = view.video(dn=snap_dn_abs,snapshot=wordspace['use_snapshot'])
			wordspace.under_development = 2
			write_wordspace(wordspace,outfile=wordspace_out)
		elif wordspace.view_mode == 'video': 
			view.video_dn,view.snapshot = wordspace.video_snaps,wordspace['use_snapshot']
		#---stage 3: run show and create snapshots if making a video
		if wordspace.under_development < 3:
			text_mode = wordspace.view_mode in ['video','snapshot']
			if wordspace.view_mode == 'snapshot': view.do('snapshot')
			view.show(quit=text_mode,text=text_mode)
			wordspace.under_development = 3
			write_wordspace(wordspace,outfile=wordspace_out)
		#---stage 4: render the video
		if wordspace.view_mode == 'video': 
			view.render(name=video_name,
				size=wordspace.video_size,duration=wordspace.get('duration',0))
		wordspace.under_development = 4 if wordspace.view_mode == 'video' else 2
	except KeyboardInterrupt as e: exception_handler(e,wordspace,all=True)
	except Exception as e: exception_handler(e,wordspace,all=True)
