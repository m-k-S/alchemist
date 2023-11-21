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
    'command': PWSCF_COMMAND + " -in Fe.pwi > Fe.pwo",
    'pseudo_dir': PPDIR,
}