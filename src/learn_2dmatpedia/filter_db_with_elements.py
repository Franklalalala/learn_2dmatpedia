from ase.visualize import view
from ase.db import connect
from ase.io import read, write
import os


def is_material_eligible(elements_list):
    """
    Check if a material meets the screening criteria:
    1. Contains only elements from the first 4 rows of the periodic table
    2. No noble gas elements
    3. No transitional elements from Sc to Ni

    Parameters:
    -----------
    elements_list : list
        List of elements in the material

    Returns:
    --------
    bool
        True if the material meets all criteria, False otherwise
    """
    # First 4 rows of the periodic table (excluding H)
    first_4_rows = set([
        'H', 'He',  # Row 1
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',  # Row 2
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',  # Row 3
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  # Row 4 (part 1)
        'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'  # Row 4 (part 2)
    ])

    # Noble gases
    noble_gases = set(['He', 'Ne', 'Ar', 'Kr'])

    # Transition metals from Sc to Ni
    transition_metals = set(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni'])

    # Check if all elements are in the first 4 rows
    if not all(element in first_4_rows for element in elements_list):
        return False

    # Check if there are no noble gases
    if any(element in noble_gases for element in elements_list):
        return False

    # Check if there are no transition metals from Sc to Ni
    if any(element in transition_metals for element in elements_list):
        return False

    return True


if os.path.exists('filtered_2DMat.db'):
    os.remove('filtered_2DMat.db')

with connect('all_2DMat.db') as db, connect('filtered_2DMat.db') as tgt_db, connect('filtered_2DMat_sample.db') as tgt_db_sample:
    for a_row in db.select():
        a_row.material_id = a_row.data.material_id
        an_atoms = a_row.toatoms()
        eligible_elements_flag = is_material_eligible(list(an_atoms.symbols))
        if eligible_elements_flag:
            tgt_db.write(a_row)
            if tgt_db.count() < 5:
                tgt_db_sample.write(a_row)
    print(tgt_db.count())
    print(tgt_db_sample.count())

