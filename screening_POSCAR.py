import os
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from io import StringIO
from collections import defaultdict
import math


def process_poscar(file_path):
    datasets = []
    with open(file_path, 'r') as f:
        while True:
            header = [f.readline().strip() for _ in range(6)]
            if not header[0]:
                break
            elements_count_line = f.readline().strip().split()
            elements_count = list(map(int, elements_count_line))
            total_atoms = sum(elements_count)
            direct_line = f.readline().strip()
            atoms_info = [f.readline().strip() for _ in range(total_atoms)]
            elements_count_str = '   '.join(map(str, elements_count))
            dataset = header + [elements_count_str] + [direct_line] + atoms_info
            datasets.append(dataset)
    return datasets


def dataset_to_structure(dataset):
    poscar_str = "\n".join(dataset)
    poscar_io = StringIO(poscar_str)
    print(poscar_io.getvalue())
    structure = Structure.from_str(poscar_io.getvalue(), fmt="poscar")
    return structure


def get_crystal_system(structure):
    lattice = structure.lattice
    a, b, c = lattice.a, lattice.b, lattice.c
    alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma
    if alpha != 90 and beta != 90 and gamma != 90:
        return "triclinic"
    elif (alpha != 90 or beta != 90 or gamma != 90):
        return "monoclinic"
    elif alpha == 90 and beta == 90 and gamma == 90:
        if a == b == c:
            return "cubic"
        elif a == b or b == c or a == c:
            return "tetragonal"
        else:
            return "orthorhombic"
    elif alpha == 90 and beta == 90 and gamma == 120:
        return "hexagonal"
    elif alpha == beta == gamma != 90:
        return "trigonal"
    else:
        return "unknown"


def compare_structures(struct1, struct2):
    species = sorted(set(site.species_string for site in struct1))
    for specie in species:
        struct1_sites = sorted((site for site in struct1 if site.species_string == specie),
                               key=lambda site: site.frac_coords[0])
        struct2_sites = sorted((site for site in struct2 if site.species_string == specie),
                               key=lambda site: site.frac_coords[0])
        for site1, site2 in zip(struct1_sites, struct2_sites):
            if not (abs(site1.frac_coords[0] - site2.frac_coords[0]) < 1e-3 and abs(
                    site1.frac_coords[1] - site2.frac_coords[1]) < 1e-3):
                return False
        z_coords1 = sorted([(site.species_string, site.frac_coords[2]) for site in struct1_sites], key=lambda x: x[1])
        z_coords2 = sorted([(site.species_string, site.frac_coords[2]) for site in struct2_sites], key=lambda x: x[1])
        z_order1 = [atom[0] for atom in z_coords1]
        z_order2 = [atom[0] for atom in z_coords2]
        if z_order1 != z_order2:
            return False
    return True


def filter_similar_structures(crystal_system_structures):
    unique_structures = []
    for i, (index, dataset, struct) in enumerate(crystal_system_structures):
        is_unique = True
        for _, _, unique_struct in unique_structures:
            if compare_structures(struct, unique_struct):
                is_unique = False
                break
        if is_unique:
            unique_structures.append((index, dataset, struct))
    return unique_structures


def save_structures_to_file(filtered_crystal_systems):
    output_dir = 'results'
    output_file = os.path.join(output_dir, 'de_ini_1')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_file, 'w') as f:
        for crystal_system, poscars in filtered_crystal_systems.items():
            for poscar_info in poscars:
                idx, dataset, structure = poscar_info
                poscar_str = "\n".join(dataset)
                f.write(poscar_str + "\n")

def is_nan(value):
    try:
        return math.isnan(float(value))
    except:
        return False
def contains_nan(data):
    for item in data:
        parts = item.split()
        for part in parts:
            if is_nan(part):
                return True
    return False

file_path = 'results.dat'
datasets = process_poscar(file_path)

crystal_systems = defaultdict(list)
for i, dataset in enumerate(datasets, start=1):
    def deal_lattice(dataset):
        dataset[2] = dataset[2].replace('-',' -')
        dataset[3] = dataset[3].replace('-',' -')
        dataset[4] = dataset[4].replace('-',' -')
        return dataset
    dataset = deal_lattice(dataset)
    if not contains_nan(dataset):
        structure = dataset_to_structure(dataset)
    else:
        continue
    crystal_system = get_crystal_system(structure)
    crystal_systems[crystal_system].append((i, dataset, structure))

filtered_crystal_systems = {cs: filter_similar_structures(structs) for cs, structs in crystal_systems.items()}

save_structures_to_file(filtered_crystal_systems)

