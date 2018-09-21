#!/usr/bin/env python

"""
AUTOPLOT DEMONSTRATION

Setting `autoplot: True` in the metadata `director` or inside an individual plot will allow you to use a few
neat features which we call "autoplot" and make use of decorated functions to manage plotting.

# How to run this demo.

1. Add the following entry to a specs file. Replaces collections and calculation with existing ones. 
   Note that this demo is not seamless with a unit test (sorry!).

plots:
  demo_autoplot:
    script: demonstrate_autoplot.py
    autoplot: True
    collections: focus
    calculation: undulations

2. Run `make go demo_autoplot`. The omnicalc output tells you what has been executed.
3. Change the function.


"""

# The routine below tells us which plots to run. Note that None means everything, 
# ... [] means nothing, otherwise this list names the decorated plot functions to run, by default.
plotrun.routine = None 

import re,textwrap

def say(msg): status(tag='tutorial',string=textwrap.fill(re.sub('(\n|\t|\s)+',' ',msg.strip()),80))

@autoload(plotrun)
def load_ignored():
	"""Dummy autoload to illustrate the perils of decorating twice."""
	print("I will never run because I am not the last-decorated autoload functio. Sad!")


@autoload(plotrun)
def load():
	"""
	Load and post-process data.
	"""
	import re,textwrap
	message = """
	Designate a "load" function of any name with the decorator `@autoload(plotrun)`.
	This function will run once when the script is started, but it will be ignored when you use the 
	`replot()` function to run the script again. The purpose of `replot()` is to redo any plots you are
	making *without* loading or performing any post-processing. This saves time. This function also
	exports all locals to globals so they do not need to be passed. Rerun this function during development
	with its name `load()`. It will only be updated after `replot()`.
	"""
	numbers = np.random.random(10)
	status(tag='load',string='The `load` function has generated some numbers:\n\n%s\n'%numbers)
	say(message)


@autoplot(plotrun)
def plot():
	"""Plot the data. Make changes and `replot()` without running the loader again."""
	message = """I am a plot function called `plot` which has been decorated with `@autoplot(plotrun)`. 
	I found numbers randomly generated in the load function and automatically exposed to globals. If you 
	run `replot()` then I will run again, even if you changed the code. Try uncommenting the line at the 
	end and then run `replot()`."""
	status(tag='load',string='The numbers are:\n\n%s\n'%numbers)
	say(message)
	#! print('edit the script and run `replot()` without exiting python or reloding. neat!')
	# uncomment the above and run `replot()` to see how it works

@autoplot(plotrun)
def plot2():
	"""Another plot function to demonstrate the routine."""
	message = """I am another plot function creatively called `plot2`. I will always run on `replot()` 
	as long as `plotrun.routine = None` is found in the script. If you set `plotrun.routine = []` instead, 
	then nothing will run. If you set `plotrun.routine = ["plot2"]` then only I will run, and the other plot 
	function `plot` will not. Note that any function which is decorated with `@autoplot(plotrun)` will run 
	automatically unless you change plotrun.routine. Any of these functions can be individually executed 
	from the command line without using the interactive mode by using `make go plot_name plot_function`. 
	"""
	say(message)

if __name__=='__main__':
	message = """If you want to iteratively code up a plot in globals, you can use the usual method, by 
	placing it in `if __name__=='__main__':` so that it will only run after we run the load function and 
	any other plot functions. When it's working correctly, it's best to make it into a plot function by
	decorating it with `@autoplot(plotrun)`. Any changes to this main loop will also take effect when 
	running `replot()`."""
	say('Now we are in the main loop.')
	say(message)
