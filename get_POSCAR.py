from collections import Counter
import numpy as np
import read_poscar
import re
import os
"""
this file mainly  about generate POSCAR file from de_ini_1 and check the centrosymmtry of 
the POSCARS, owing to the two-dimensional boundary problem
"""
def get_q2d_status(file_path='../test.sh'):
    with open(file_path, 'r',encoding='utf-8',errors='ignore') as file:
        for line in file:
            if line.strip().startswith('Q2D='):
                value = line.strip().split('=')[1].strip()
                if value.upper().startswith('T'):
                    return True
                elif value.upper().startswith('F'):
                    return False
    return False  

def generate_equal_atoms(atom_info,is_2D):
    if is_2D:
        if atom_info[0] == 0 and atom_info[1] != 0:
            equal_atom_one = atom_info + np.array([1,0,0])
            sub_atom_info = (atom_info+equal_atom_one)/2
        elif atom_info[1] == 0 and atom_info[0] != 0:
            equal_atom_one = atom_info + np.array([0,1,0])
            sub_atom_info = sub_atom_info = (atom_info+equal_atom_one)/2
        elif atom_info[0] == 0 and atom_info[1] == 0:
            equal_atom_one = atom_info + np.array([0, 1, 0])
            equal_atom_two= atom_info + np.array([1, 0, 0])
            equal_atom_three = atom_info + np.array([1, 1, 0])
            sub_atom_info = (atom_info+equal_atom_one+equal_atom_three+equal_atom_two)/4
        else:
            sub_atom_info = atom_info
    else:
        equal_array = np.array([0, 0, 0])
        if atom_info[0] == 0:
            equal_array += np.array([1, 0, 0])
        if atom_info[1] == 0:
            equal_array += np.array([0, 1, 0])
        if atom_info[2] == 0:
            equal_array += np.array([0, 0, 1])

        if np.any(equal_array != 0):
            sub_atom_info = atom_info + equal_array / 2
        else:
            sub_atom_info = atom_info
    return sub_atom_info

def check_mapping_atoms(atoms_info_classified,center_point,atoms_num):
    def remove_even_duplicates(nums):
        counts = Counter(nums)
        no_even_duplicates = [num for num, count in counts.items() if count % 2 != 0]
        return no_even_duplicates
    for atoms_info in atoms_info_classified:
        #cal the L2 norm table
        l2_norm_table = []
        for atom_info in atoms_info:
            l2_norm_table.append(round(np.linalg.norm(atom_info - center_point),2))
        norm = remove_even_duplicates(l2_norm_table)
        try:
            norm.remove(0.00)
            if len(norm) != 0:
                return False
        except(ValueError):
            if len(norm) != 0:
                return False
    return True

if __name__ == '__main__':
    file_path = '../test.sh'
    population_value = None
    is_2D = get_q2d_status()
    #find atoms accounts
    regex_pattern = r'NumberOfElements=([\d\s]+)'
    sum_of_numbers = 0
    with open(file_path, 'r',encoding='utf-8',errors='ignore') as file:
        file_contents = file.read()
        matches = re.finditer(regex_pattern, file_contents)
        for match in matches:
            numbers = match.group(1).split()
            sum_of_numbers += sum(int(number) for number in numbers)
    with open('de_ini_1','r') as file:
        content = file.readlines()
        population_value = int(len(content)/(8+sum_of_numbers))

    with open('de_ini_1', 'r') as f1:
        content = f1.readlines()
        for i in range(population_value + 1)[1:]:
            POSCAR = content[(8 + sum_of_numbers) * (i - 1):(8 + sum_of_numbers) * i]
            f2 = open("POSCAR" + str(i), 'w')
            f2.writelines(POSCAR)
            f2.close()

    #check symmtry
    for i in range(population_value):
        POSCAR = read_poscar.read_poscar('POSCAR'+str(i+1))
        atoms_info_list = POSCAR['atoms_axis']
        atoms_info_classified = []
        atoms_info_unclassified = []
        for atoms in atoms_info_list.values():
            elements_atoms_info = []
            for atom in atoms:
                elements_atoms_info.append(generate_equal_atoms(atom,is_2D))
                atoms_info_unclassified.append(generate_equal_atoms(atom,is_2D))
            atoms_info_classified.append(elements_atoms_info)
        center_point =np.mean(atoms_info_unclassified,axis=0)
        judge = check_mapping_atoms(atoms_info_classified,center_point,len(atoms_info_unclassified))
        if not judge:
            os.system("rm POSCAR"+str(i+1))

