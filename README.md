# Current Features

## Interactivity

- Start interactive QE session in Jupyter or connect to remote Jupyter server

## QE 

- Structural relaxation
- Band gap calculation 

# To-Do

## QE

- Phonon band structure / DOS calculation workflow
- Electron band structure / DOS calculation workflow
- Seebeck coefficient calculation
- 

## JAX PySCF
- Re-write critical functions of PySCF in JAX for automatic differentiation and TPU acceleration support

## Hybrid KMC
- Automatically compute activation energies, rate constants using ab-initio calculations

## Structure
- Class for building atomic environments
    - Automatically pull structures from Materials Project if they exist
- Config class for QE input parameters
- Config class for QE output parameters

## Productivity Tools

### LLM
- OpenAI Embedding for QE documentations
- ChatGPT interface for generating QE input files for specific tasks (e.g. band gap calculations)
    - Linter / sandbox for testing QE input file validity
- Automatic plotting of QE output files

### Compute Cluster
- Queueing system for sending jobs to compute clusters
- Mobile notifications and job submission

## Documentation 
- Document how to install, set up QE and point to the pw.x executable

## Supercloud 
- Build QE from source (Github)
- Include the following in .jupyter/llsc_notebook_bashrc:
    export PYTHONNOUSERSITE=True
    export PYTHONPATH=/home/gridsan/$USER/.conda/envs/$ENV/bin/python
    export PATH="${PATH}:/home/gridsan/$USER/.conda/envs/$ENV/bin"
- Launch Jupyter notebook with as many cores as desired
- Create .env file with PWSCF_COMMAND and PSEUDOPOTENTIALS pointing to the correct files/directories
- Run this at the top of the notebook:
    os.environ['SLURM_TASKS_PER_NODE'] = os.environ['SLURM_CPUS_PER_TASK']
    os.environ['SLURM_CPUS_PER_TASK'] = "1"
