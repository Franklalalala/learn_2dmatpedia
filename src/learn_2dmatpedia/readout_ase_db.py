from ase.visualize import view
from ase.db import connect
from ase.io import read, write
import os


with connect('materials.db') as db:
    # os.makedirs('cif')
    # os.chdir('cif')
    for idx, a_row in enumerate(db.select()):
        an_atoms = a_row.toatoms()
        a=1
        write(str(idx) + '.cif', an_atoms, format='cif')


        # view(an_atoms)
