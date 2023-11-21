import ase.units as units

# Lattice constants in Angstroms
a_bcc = 2.87 * units.Angstrom
a_hcp = 2.45 * units.Angstrom 
c_hcp = 3.93 * units.Angstrom

k_points = 3

# Building the crystal structures
iron_bcc = bulk('Fe', 'bcc', a=a_bcc)
iron_hcp = bulk('Fe', 'hcp', a=a_hcp, c=c_hcp)