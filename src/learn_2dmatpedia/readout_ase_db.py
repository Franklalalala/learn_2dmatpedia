from ase.visualize import view
from ase.db import connect
from ase.io import read, write
import os


source_id_set = set()
material_id_set = set()
with connect('filtered_2DMat_sample.db') as db:
    for idx, a_row in enumerate(db.select()):
        an_atoms = a_row.toatoms()
        # a_row.material_id = a_row.data.material_id
        if 'material_id' in a_row.data.keys():
            material_id_set.add(a_row.data['material_id'])
        if 'source_id' in a_row.data.keys():
            source_id_set.add(a_row.data['source_id'])
    print(db.count())

print(len(source_id_set))
print(len(material_id_set))
print(idx)