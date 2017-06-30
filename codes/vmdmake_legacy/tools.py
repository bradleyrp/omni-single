#!/usr/bin/python

import os,sys,re,json

def yamlparse(text,style=None):

	"""
	A function which reads the settings files in yaml format.
	DEVELOPMENT NOTE: this will ignore lines if your forget a colon. Needs better error checking.
	"""
	
	unpacked = {}
	#---evaluate code blocks first
	regex_block_standard = r"^\s*([^\n]*?)\s*(?:\s*:\s*\|)\s*([^\n]*?)\n(\s+)(.*?)\n(?!\3)"
	regex_block_tabbed = r"^\s*([^\n]*?)\s*(?:\s*:\s*\|)\s*\n(.*?)\n(?!\t)"
	if style == 'tabbed': regex_block = regex_block_tabbed
	else: regex_block = regex_block_standard
	regex_line = r"^\s*(.*?)\s*(?:\s*:\s*)\s*(.+)$"
	while True:
		blockoff = re.search(regex_block,text,re.M+re.DOTALL)
		if not blockoff: break
		if style == 'tabbed': key,block = blockoff.groups()[:2]
		else: 
			#---collect the key, indentation for replacement, and value
			key,indent,block = blockoff.group(1),blockoff.group(3),''.join(blockoff.groups()[1:])
		#---alternate style does multiline blocks with a single tab character
		if style == 'tabbed': compact = re.sub("(\n\t)",r'\n',block.lstrip('\t'),re.M)
		#---remove indentations and newlines (default)
		else: compact = re.sub('\n','',re.sub(indent,'',block))
		unpacked[re.sub(' ','_',key)] = compact
		#---remove the block
		text = re.sub(re.escape(text[slice(*blockoff.span())]),'',text)
	while True:
		line = re.search(regex_line,text,re.M)
		if not line: break
		key,val = line.groups()
		assert key not in unpacked
		unpacked[re.sub(' ','_',key)] = val
		text = re.sub(re.escape(text[slice(*line.span())]),'',text)
	#---evaluate rules to process the results
	for key,val in unpacked.items():
		#---store according to evaluation rules
		try: val = eval(val)
		except: pass
		if type(val)==list: unpacked[key] = val
		elif type(val)==str:
			if re.match('^(T|t)rue$',val): unpacked[key] = True
			elif re.match('^(F|f)alse$',val): unpacked[key] = False
			#---! may be redundant with the eval command above
			elif re.match('^[0-9]+$',val): unpacked[key] = int(val)
			elif re.match('^[0-9]*\.[0-9]*$',val): unpacked[key] = float(val)
			else: unpacked[key] = val
		else: unpacked[key] = val
	return unpacked

def nospaces(text):

	"""
	Working with keys processed by yamlparse enforces a no-spaces rule.
	"""

	return re.sub(' ','_',text)

"""
NOTE
we recapitulate the WordSpace, bash, and ... here
this allows the vmdmake.view_routine to mimic some of the automacs behavior 
	that is useful for managing the video directories
"""

class WordSpace(dict):
	def __getattribute__(self,key):
		if key in self: return self[key]
		else: return dict.__getattribute__(self,key)
	def __setattr__(self,key,value):
		if key=='under_development' and type(value)==int:
			#from amx.base.metatools import write_wordspace
			#write_wordspace(wordspace,wordspace.wordspace_location)
			dict.__setitem__(self,key,value)
		elif key not in dict.__dict__: dict.__setitem__(self,key,value)
		else: raise Exception('[ERROR] cannot set the %s key in the wordspace'%key)
	def __getitem__(self,key):
		if key=='last' and key not in self:
			raise Exception("".join([
				"[ERROR] wordspace['last'] is not defined ...",
				"[ERROR] it is likely that you started AMX from a downstream step"]))
		elif key not in self:
			if key == 'watch_file':
				print('[ERROR] watch_file not defined yet so returning "WATCHFILE"')
				return "WATCHFILE"
			elif key == 'under_development': return False
		return dict.get(self,key)

gmxpaths = {}

def bash(command,log=None,cwd=None,inpipe=None):

	"""
	Run a bash command
	"""
	
	cwd = wordspace['step'] if cwd == None else cwd
	if log == None: 
		if inpipe: raise Exception('under development')
		kwargs = dict(cwd=cwd,shell=True,executable='/bin/bash',
			stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		subprocess.call(command,**kwargs)
	else:
		output = open(cwd+'log-'+log,'w')
		kwargs = dict(cwd=cwd,shell=True,executable='/bin/bash',
			stdout=output,stderr=output)
		if inpipe: kwargs['stdin'] = subprocess.PIPE
		proc = subprocess.Popen(command,**kwargs)
		if not inpipe: proc.communicate()
		else: proc.communicate(input=inpipe)

#---CORRECT LOCATION?

def write_wordspace(wordspace,outfile=None):

	"""
	In addition to saving checkpoints permanently in the log, we also drop the wordspace into a json file
	for rapid development.
	"""

	outfile = 'wordspace.json' if not outfile else outfile
	with open(outfile,'w') as fp: json.dump(wordspace,fp)

def exception_handler(e,wordspace,all=False,outfile=None):

	"""
	Report an error concisely to the terminal to avoid overwhelming the user.
	"""

	exc_type, exc_obj, exc_tb = sys.exc_info()
	fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
	#---modified from automacs to work without report
	def report(x,tag=None): print('%s%s'%(tag.upper() if tag else '',x))
	report('%s in %s at line %d'%(str(exc_type),fname,exc_tb.tb_lineno),tag='error')
	report('%s'%e,tag='error')
	if all:
		import traceback
		report(re.sub('\n','\n[TRACEBACK] ',traceback.format_exc()),tag='traceback')
	write_wordspace(wordspace,outfile=outfile)
	sys.exit(1)
