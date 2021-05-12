# How to reproduce these calculations

This README was added on 2021.05.10 in order to provide general instructions for reproducing these calculations. The status of the "factory" codes (specifically the calculation portion which we called "omnicalc") is as follows:

- the calculations for specific projects (atomistic bilayers, coarse-grained bilayers, and some mesoscale calculations) are av

## General instructions

### 1. Get the data

We start with existing simulation data, and produce two additional kinds of postprocessing data summarized below:

1. The original simulation data (typically produced with "automacs") have a folder with individual simulation steps named e.g. `s02-production` (regex: `s([0-9]{2})-(.*?)$`) inside of a parent folder. The segments of each simulation are named with the typical GROMACS naming scheme, for example `md.part0001.xtc` or `md.part0004.tpr`.
2. There are two kinds of post-processed data. The first are subsampled trajectories which take the form `v1021.2000000-12000000-1.proteins.pbcmol.gro` and are paired with an `xtc` and `ndx` file. These will be deposited in the `post` folder later and can be viewed directly in VMD. Sometimes we call these simulation "slices". The name includes the selection (`proteins`) and a GROMACS-specific wrapping parameter (almost always `pbcmol`, which prevents molecules from wrapping across periodic boundaries).
3. The second post-processing data consist of pairs of files with either a `dat` or `spec` suffix. The `spec` files include metadata that describes the calculation. If you perform a parameter sweep, the files are numbered and the specifics input parameters for each calculation in the sweep (for example a cutoff value, a VMD-style selection, etc). The corresponding `dat` files are hdf5 files produced by our analysis scripts. An example file is `v1021.2000000-12000000-1.lipid_abstractor.n0.dat`. This example file contains the results of the `lipid_abstractor` calculation (which comes from a script at `calcs/lipid_abstractor.py` which we will discuss later). The `n0` means this is "row 0" in our metaphorical database. The metadata are found in `v1021.2000000-12000000-1.lipid_abstractor.n0.spec`.

The following is an example of the `spec` file:

```
{"meta": {"spec_version": 3, "sn": "v1021"}, "slice": {"skip": 1, "end": 12000000, "start": 2000000, "pbc": "mol", "group": "lipids", "sn": "v1021", "slice_name": "current"}, "calc": {"name": "undulations", "specs": {"upstream": {"lipid_abstractor": {"selector": "lipid_com", "separator": {"cluster": true}}}, "grid_spacing": 0.5}}}
```

This is written in the JSON format, and provides information on the simulation "slice", the grid spacing for the interpolation we use to compute undulations, and a few other related parameters (for example, we select the center of each lipid by mass `lipid_com`).

I typically name the simulations with numbers. In the example that motivates these instructions, we will discuss a BAR-domain simulation in `v1021` and reproduce undulation spectra from the `v650` simulation labeled `4xENTH`, which was one of subjects of [Bradley-2016](https://www.pnas.org/content/113/35/E5117.abstract).

In the following steps, we will start with item #1 abve (the original simulation data), and reproduce the undulation calculation.

### 2. Install an anaconda environment

There are many ways to install an environment. In the past Ryan Bradley and Joe Jordan have used the "factory" to manage software environments, serve a web portal for starting simulations, and lastly, to glue together the simulation code. At some point we forked the code into two pices. In this step, we will install an Anaconda environment.

You are welcome to install the environment manually with the following dependency file (in YAML). You can do this with a pre-installed Anaconda distribution with `conda env update --file <env.yaml> -p ./path/to/env`.

```
name: bpcfac
dependencies:
  - python==3.8.5
  - pip
  - pip:
    - django<2
    - django-extensions
    - pyyaml
    - nbformat
    - image # not sure which is PIL
    - Pillow # this might also be PIL
    - numpy
    - ipdb
    - mdanalysis
    - h5py
    - scikit-learn
    - seaborn
```

Alternately, you can use the new version of the factory to automatically install Anaconda and make the environments with a couple commands. Note that I have labeled this the `factory` because it uses an updated code, but it will be distinct from the older version below, which we call the `legacy-factory`. The environment name here is `bpcfac` but you are welcome to use a different name. In the following instructions we use a cat/EOF trick to write a text file from the shell so the instructions are (more) concise.

```
git clone https://github.com/bradleyrp/factory -b streamline ./factory
cd factory
cat > specs/env_bpcfac.yaml << EOF
name: bpcfac
dependencies:
  - python==3.8.5
  - pip
  - pip:
    - django<2
    - django-extensions
    - pyyaml
    - nbformat
    - image # not sure which is PIL
    - Pillow # this might also be PIL
    - numpy
    - ipdb
    - mdanalysis
    - h5py
    - scikit-learn
    - seaborn
EOF
make conda specs/env_bpcfac.yaml
# to use this environment, you can run this command from anywhere
source path/to/factory env.sh bpcfac
```

You will need to activate this environment for all subsequent steps with `conda activate path/to/env`.

### 3. Clone the (legacy) factory

In order to use some of these codes, you should make a copy of the `legacy-factory` code. In the past, we used this for a few different purposes (namely a small simulation portal). Now it still functions as a method for gluing together our data and simulation scripts.

```
# select a location for storing both data and code (this depends on your local machine)
cd /data/rbradley
git clone http://github.com/biophyscode/factory ./legacy-factory
cd legacy-factory
# this is the factory root path
export FACTORY_ROOT=$PWD
```

Next we create a "connection" file to set some paths and point to the calculation repositories. You can use the cat/EOF trick below to write the `connections/banana.yaml` file or just make the file with vim. Before you write this file, change the `route_to_data` items and `spot_directory` so that they provide the parent directory and folder name for your source data. In the example below, my source data are located at `/data/rbradley/data/2021.02.15-banana-set`. We are calling this project "banana" but the factory was designed to host many projects, and the name is arbitrary.

```
mkdir -p connections
cat > connections/banana.yaml <<EOF
banana:
  enable: true 
  site: site/PROJECT_NAME  
  calc: calc/PROJECT_NAME
  calc_meta_filters: []
  repo: http://github.com/bradleyrp/omni-single
  database: data/PROJECT_NAME/db.factory.sqlite3
  post_spot: data/PROJECT_NAME/post
  plot_spot: data/PROJECT_NAME/plot
  simulation_spot: data/PROJECT_NAME/sims
  spots:
    # colloquial name for the default "spot" for new simulations given as simulation_spot above
    sims:
      # name downstream postprocessing data according to the spot name (above) and simulation folder (top)
      # the default namer uses only the name (you must generate unique names if importing from many spots)
      namer: "lambda name,spot=None: name"
      # parent location of the spot_directory (may be changed if you mount the data elsewhere)
      route_to_data: /data/rbradley/data/
      # path of the parent directory for the simulation data
      spot_directory: 2021.02.15-banana-set
      # rules for parsing the data in the spot directories
      regexes:
        # each simulation folder in the spot directory must match the top regex
        top: '(.+)'
        # each simulation folder must have trajectories in subfolders that match the step regex (can be null)
        # note: you must enforce directory structure here with not-slash
        step: '([stuv])([0-9]+)-([^\/]+)'
        # each part regex is parsed by omnicalc
        part: 
          xtc: 'md\.part([0-9]{4})\.xtc'
          trr: 'md\.part([0-9]{4})\.trr'
          edr: 'md\.part([0-9]{4})\.edr'
          tpr: 'md\.part([0-9]{4})\.tpr'
          # specify a naming convention for structures to complement the trajectories
          structure: '(system|system-input|structure)\.(gro|pdb)'
EOF
```

Next we have to run some commands:

```
make set omnicalc="http://github.com/biophyscode/omnicalc"
make set automacs="http://github.com/biophyscode/automacs"
mkdir -p logs
make connect
```

Now the factory is ready. You will use the following path to run all of the scripts:

```
cd $FACTORY_ROOT/calcs/banana
export $CALC_ROOT=$PWD
```

We have a couple more files to set. First we set a config file with a GROMACS binary suffix. If you are using an HPC system, you probably need to specify `gmx_mpi` to run GROMACS. We also have a flag for wrapping GROAMCS with `mpiexec`, which is required on the cluster I tested this with.

```
cat > gromacs_config.py <<EOF
#!/usr/bin/python
machine_configuration = {
	"LOCAL": dict(gpu_flag="auto",
		use_mpiexec=True,
		gmx_series=5,
		suffix="_mpi"),}
EOF
```

### 5. Refresh the codes

In case we need to update this code, you can pull new changes with these commands:

```
git pull
git -C calcs pull
```

### 6. Configure the calculation

This code uses a YAML file to set up the calculation in order to allow for more flexible loops and parameter sweeps. For this example, I have prepared a file named `calcs/specs/curvature_repro_20210501.yaml`. I recommend inspecting it. To use this file, run this command:

```
make set meta_filter curvature_repro_20210501.yaml
```

I strongly recommend maintaining your own specs files as you add data to your dataset and expand and refine your analyses. These individual files, combined with the git commit hashes for calculation codes, can tell you exactly how to reproduce a calculation later on. 

### 7. Compute undulation spectra

At this point, we now have several components in place:

- We have installed an environment with the right supporting libraries.
- The `legacy-factory` program is holding our calculation codes in a subfolder (`calcs/banana`).
- The calculation codes are located in `calcs/banana/calcs` and come from a shared repository.
- These codes include a "specs" (specifications) file at `calcs/specs/curvature_repro_20210501.yaml`.

To reproduce the calculations, we will run the following commands, ideally with a few processors in a batch job.

```
make compute
make plot undulations
```

Further instructions are included in the specs file for this specific example.
