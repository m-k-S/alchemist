from ase import Atoms
import ase.units as units
from ase.io import read
from ase.build import bulk
from ase.calculators.espresso import Espresso

import numpy as np
import matplotlib.pyplot as plt

import os
import dotenv

dotenv.load_dotenv()
PWSCF_COMMAND = os.environ.get("PWSCF_COMMAND")
PPDIR = os.environ.get("PSEUDOPOTENTIALS")

outdir = "./results"
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Define system and electron parameters
system_params = {
    'smearing': 'mp',
    'occupations': 'smearing',
    'starting_magnetization(1)': 0.7,
    'nat': 2,
    'degauss': 0.02,
    'nspin': 2,
    'ntyp': 1,
    'ibrav': 0,
    'ecutwfc': 30 * units.eV,  # PW cutoff
    'ecutrho': 240 * units.eV,  # Charge cutoff
}

electron_params = {
    'diagonalization': 'david',
    'mixing_beta': 0.5,
    'conv_thr': 1e-07,
}

# Merge parameters into the input_data dictionary
input_data = {
    'system': system_params.copy(),
    'electrons': electron_params.copy()
}

# Lattice constants in Angstroms
a_bcc = 2.87 * units.Angstrom
a_hcp = 2.45 * units.Angstrom 
c_hcp = 3.93 * units.Angstrom

k_points = 3

# Building the crystal structures
iron_bcc = bulk('Fe', 'bcc', a=a_bcc)
iron_hcp = bulk('Fe', 'hcp', a=a_hcp, c=c_hcp)

# Setup DFT calculation parameters
pseudopotentials = {'Fe': 'Fe.pbe-nd-rrkjus.UPF'}  # Example, change as needed

calc_params = {
    'calculation': 'vc-relax',
    'tprnfor': True,
    'tstress': True,
    'input_data': {
        'system': {
            'ecutwfc': 30 * units.eV,  # PW cutoff
            'ecutrho': 240 * units.eV,  # Charge cutoff
        },
        'electrons': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-07,
        }
    },
    'pseudopotentials': pseudopotentials,
    'kpts': (k_points, k_points, k_points),  # k-points grid, adjust as needed
    'parallel': 'all',
    'directory': outdir,  # Custom directory for calculation files
    'label': "Fe",  # Prefix for the filenames
    'logfile': "Fe.log",  # Logfile name
    'command': "mpirun -np 42 " + PWSCF_COMMAND + " -in Fe.pwi > Fe.pwo",
    'pseudo_dir': PPDIR,
}

energies_bcc = []
energies_hcp = []
lattice_params_bcc = []
lattice_params_hcp = []
for lattice in ['bcc', 'hcp']:
    for k in [3, 4, 5, 6]:
        calc_params['kpts'] = (k, k, k)

        # Attach the calculator to the structures
        calculator = Espresso(**calc_params)

        if lattice == 'bcc':
            iron_bcc.set_calculator(calculator)
            energy_bcc = iron_bcc.get_potential_energy()
            energies_bcc.append(energy_bcc)

            final_structure = read('results/Fe.pwo') 
            lattice_params_bcc.append(final_structure.cell.lengths())


        elif lattice == 'hcp':
            iron_hcp.set_calculator(calculator)
            energy_hcp = iron_hcp.get_potential_energy()
            energies_hcp.append(energy_hcp)

            final_structure = read('results/Fe.pwo') 
            lattice_params_hcp.append(final_structure.cell.lengths())
