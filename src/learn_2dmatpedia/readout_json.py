import json
import numpy as np
from ase import Atoms
from ase.db import connect
import os

# Load the data from db.json
materials = []
with open('db.json', 'r') as file:
    for line in file:
        materials.append(json.loads(line))
print(f"Found {len(materials)} materials")

# Display the first 10 entries (simplified view)
print("\nFirst 10 materials:")
for i in range(min(10, len(materials))):
    print(f"{i + 1}. {materials[i].get('material_id', 'N/A')} - {materials[i].get('formula_pretty', 'N/A')}")

# Convert to ASE database
ase_db_filename = 'all_2DMat.db'
# Remove existing database if it exists
if os.path.exists(ase_db_filename):
    os.remove(ase_db_filename)

# Create a new database
db = connect(ase_db_filename)

# Convert first 10 entries to ASE database
for i in range(min(1000000000, len(materials))):
    material = materials[i]

    # Extract structure information
    structure = material.get('structure', {})
    lattice = structure.get('lattice', {})
    sites = structure.get('sites', [])

    # Extract cell matrix from lattice
    cell = lattice.get('matrix', [[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    # Extract atomic positions and symbols
    positions = []
    symbols = []

    for site in sites:
        xyz = site.get('xyz', [0, 0, 0])
        positions.append(xyz)

        species = site.get('species', [])
        if species and len(species) > 0:
            element = species[0].get('element', '')
            symbols.append(element)

    # Create ASE Atoms object
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=[True, True, True])

    # Collect metadata to store with the atoms
    key_value = {
        'material_id': material.get('material_id', ''),
        'formula_pretty': material.get('formula_pretty', ''),
        'formula_reduced_abc': material.get('formula_reduced_abc', ''),
        'formula_anonymous': material.get('formula_anonymous', ''),
        'nelements': material.get('nelements', 0),
        'elements': material.get('elements', []),
        'chemsys': material.get('chemsys', ''),
        'sg_number': material.get('sg_number', 0),
        'sg_symbol': material.get('sg_symbol', ''),
        'source_id': material.get('source_id', ''),
        'creation_task_label': material.get('creation_task_label', ''),
        'discovery_process': material.get('discovery_process', ''),
    }

    # Add energy properties if available
    if 'thermo' in material:
        thermo = material['thermo']
        key_value['energy'] = thermo.get('energy', None)
        key_value['energy_per_atom'] = thermo.get('energy_per_atom', None)
        key_value['energy_vdw'] = thermo.get('energy_vdw', None)
        key_value['energy_vdw_per_atom'] = thermo.get('energy_vdw_per_atom', None)

    # Add bandgap information if available
    if 'bandgap' in material:
        key_value['bandgap'] = material['bandgap']

    # Add magnetism information if available
    if 'total_magnetization' in material:
        key_value['total_magnetization'] = material['total_magnetization']

    # Add decomposition and exfoliation energy if available
    if 'decomposition_energy' in material:
        key_value['decomposition_energy'] = material['decomposition_energy']

    if 'exfoliation_energy_per_atom' in material:
        key_value['exfoliation_energy_per_atom'] = material['exfoliation_energy_per_atom']

    # Store atoms in database with metadata
    db.write(atoms, data=key_value)

print(f"\nConverted first 10 entries to ASE database: {ase_db_filename}")
print("To access this database, use:")
print("from ase.db import connect")
print(f"db = connect('{ase_db_filename}')")
print("for row in db.select():")
print("    print(row.material_id, row.formula_pretty)")
print("    atoms = row.toatoms()")