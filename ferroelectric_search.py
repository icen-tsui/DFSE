import os
import random
import re
import shutil
import sys
import numpy as np
import yaml
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import math
import read_poscar
from get_POSCAR import generate_equal_atoms


def perturbation_axis_mapping(perturbation_value=0.05,FQFE:bool=False,crystallization:int=0,Q2D:bool=False):
    if FQFE:
        x,y,z = 0,0,0
        perturbation_value_list = [[0.30, 0.30], [0.25, 0], [0.50, 0], [0.13, 0.13]]
        if crystallization == 7:
            perturbation_value_list.append([0.25,0.13])
        perturbation_value = random.choice(perturbation_value_list)
        perturbation_value = [value * random.choice([1, -1]) for value in perturbation_value]
        if Q2D:
            z = 0
            if random.choice([True,False]):
                x,y = perturbation_value
            else:
                y,x = perturbation_value
        else:
            zero_index = random.choice([0,1,2])
            if zero_index == 0:
                y,z = perturbation_value
            elif zero_index == 1:
                x,z=perturbation_value
            else:
                x,y =perturbation_value
    else:
        x, y, z = 0, 0, 0
        axis = random.choice(['x','y','z'])

        if axis == 'x':
            x = random.choice([1, -1])*perturbation_value
        elif axis == 'y':
            y = random.choice([1, -1])*perturbation_value
        elif axis == 'z':
            z = random.choice([1, -1])*perturbation_value
    return np.array([x, y, z])
def add_perturbation(origin_POSCAR:str,Q2D:bool,FQFE:bool,loop_index:int,coefficient:float=0.05):
    global folder_name
    POSCAR = read_poscar.read_poscar(origin_POSCAR)

    ## choose random atom
    if Q2D==True:
        pertubation_value = perturbation_axis_mapping(perturbation_value=coefficient,FQFE=FQFE,crystallization=None,Q2D=True)
    else:
        structure = Structure.from_file(origin_POSCAR)
        def check_cubic(structure):
            lattice = structure.lattice
            a, b, c = lattice.a, lattice.b, lattice.c
            alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma
            if (alpha == 90 and beta == 90 and gamma == 90) and (a==b==c):
                return 7
            else:
                return 0
        crystallization = check_cubic(structure)
        pertubation_value = perturbation_axis_mapping(perturbation_value=coefficient,FQFE=FQFE,crystallization=crystallization,Q2D=False)

    dir_corrdination_axiss = []
    atoms_kind = []
    atoms_num = []
    for atoms_axis in POSCAR['atoms_axis']:
        for atom_axis in POSCAR['atoms_axis'][atoms_axis]:
            try:
                if  atoms_axis[0] != atoms_kind[-1]:
                    atoms_kind.append(atoms_axis[0])
                    atoms_num.append(atoms_axis[1])
            except(IndexError):
                atoms_kind.append(atoms_axis[0])
                atoms_num.append(atoms_axis[1])
            dir_coordination_axis = atom_axis
            dir_corrdination_axiss.append(dir_coordination_axis)
    #random atom random add perturbation value
    atoms_ID = None
    mapping_info = None
    old_atoms_id = None
    with open(origin_POSCAR.split('/')[0]+r'/sym.log') as sym_info:
        sym_info = yaml.safe_load(sym_info)
        try:
            mapping_info = sym_info['atom_mapping']
        except(TypeError):
            mapping_info= None
    if mapping_info == None:
        atoms_id = random.randint(0,sum(atoms_num)-1)
    else:
        atoms_ID = list(set(mapping_info.values()))
        atoms_id = random.choice(atoms_ID)
    choice = random.randint(0,1)
    if choice == 0 or FQFE == True:
        dir_corrdination_axiss[atoms_id] = dir_corrdination_axiss[atoms_id] + pertubation_value
    if choice == 1 and FQFE == False:
        try:
            old_atoms_id = atoms_id
            dir_corrdination_axiss[atoms_id] = dir_corrdination_axiss[atoms_id] + (pertubation_value/2)
            temp = list(range(sum(atoms_num)))
            temp.remove(atoms_id)
            atoms_id = random.choice(temp)
            dir_corrdination_axiss[atoms_id] = dir_corrdination_axiss[atoms_id] + (pertubation_value/2)
        except(TypeError):
            dir_corrdination_axiss[atoms_id] = dir_corrdination_axiss[atoms_id] + pertubation_value

    #make new folder to storage POSCAR and make vasp_calculation
    match = re.match(r'.*POSCAR(\d+).*',origin_POSCAR)
    if match:
        inner_number = match.group(1)
        folder_name=f'POSCAR-{inner_number}-{loop_index}'
        with open(f'./perturbation.log','a+') as f:
            f.write(folder_name+'\n')
            if choice == 1 and FQFE == False:
                f.write('atoms_index1:'+str(atoms_id)+'\n')
                f.write(str(list(pertubation_value/2))+'\n') 
                f.write('atoms_index2:'+str(old_atoms_id)+'\n')
                f.write(str(list(pertubation_value/2))+'\n') 
            f.write('atom_index:'+str(atoms_id)+'\n')
            f.write(str(list(pertubation_value))+'\n') 
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
        else:
            shutil.rmtree(folder_name)
            os.makedirs(folder_name)

    # write inforamtion into new POSCAR , create new files to storage
    with open(f'./{folder_name}/POSCAR','w') as f1:
        f1.write("POSCAR"+'\n')
        f1.write(str(POSCAR['scale_factor'])+'\n')
        for lattice in POSCAR['lattice_matrix']:
            formatted_lattice = '        '.join([f"{x:.11f}" for x in lattice])
            f1.write(formatted_lattice+'\n')
        f1.write(' '.join([str(x) for x in atoms_kind])+'\n')
        f1.write(' '.join([str(x) for x in atoms_num])+'\n')
        f1.write("Direct"+'\n')
        for atom_axis in dir_corrdination_axiss:
            formatted_axis = '        '.join([f"{x:.11f}" for x in atom_axis])
            f1.write(formatted_axis+'\n')
        f1.close()
    return None

if sys.argv[1]  == 'generate_refer_phase':
    file_path = '../test.sh'
    #find atoms accounts
    regex_pattern = r'NumberOfElements=([\d\s]+)'
    sum_of_numbers = 0
    with open(file_path, 'r',errors='ignore') as file:
        file_contents = file.read()
        matches = re.finditer(regex_pattern, file_contents)
        for match in matches:
            numbers = match.group(1).split()
            sum_of_numbers += sum(int(number) for number in numbers)
    with open('de_ini_1',errors='ignore') as file:
        content = file.readlines()
        population_value = int(len(content)/(8+sum_of_numbers))

    with open('de_ini_1', 'r',errors='ignore') as f1 :
        content = f1.readlines()
        for i in range(population_value+1)[1:]:
            POSCAR = content[(8+sum_of_numbers)*(i-1):(8+sum_of_numbers)*i]
            f2 = open("POSCAR"+str(i),'w')
            f2.writelines(POSCAR)
            f2.close()
            os.system("mkdir refer-POSCAR"+str(i))
            os.system("cp POSCAR"+str(i)+" refer-POSCAR"+str(i)+"/POSCAR")
            os.system("rm POSCAR"+str(i))
        f1.close()

elif sys.argv[1]  == 'generate_ferroelectric_phase':
    ferro_phase_num = sys.argv[2]
    file_path = '../test.sh'
    # find atoms accounts
    regex_pattern = r'NumberOfElements=([\d\s]+)'
    sum_of_numbers = 0
    with open(file_path, 'r',errors='ignore') as file:
        file_contents = file.read()
        matches = re.finditer(regex_pattern, file_contents)
        for match in matches:
            numbers = match.group(1).split()
            sum_of_numbers += sum(int(number) for number in numbers)


    # search Q2D switch
    Q2D = False
    regex_pattern_Q2D = r'^\s*Q2D\s*=\s*\.\s*true\s*'
    pattern_Q2D = re.compile(regex_pattern_Q2D, re.IGNORECASE)

    # search FQFE switch
    FQFE = False
    regex_pattern_FQFE = r'^\s*FQFE\s*=\s*\.\s*true\s*'
    pattern_FQFE = re.compile(regex_pattern_FQFE, re.IGNORECASE)

    # search coefficient value
    coefficient = 0.05
    regex_pattern_coefficient = r'^\s*COEFFICIENT\s*=\s*([\d.]+)'
    pattern_coefficient = re.compile(regex_pattern_coefficient, re.IGNORECASE)
    with open(file_path, 'r',errors='ignore') as file:
        file_contents = file.readlines()
        for line in file_contents:
            if not line.strip().startswith('#'):
                # search Q2D switch
                if pattern_Q2D.search(line):
                    Q2D = True

                # switch FQFE switch
                if pattern_FQFE.search(line):
                    FQFE = True

                # search coefficient value
                match = pattern_coefficient.search(line)
                if match:
                    coefficient = float(match.group(1))

    with open('de_ini_1',errors='ignore') as file:
        content = file.readlines()
        population_value = int(len(content)/(8+sum_of_numbers))
    for index in range(int(population_value) + 1)[1:]:
        for loop_index in range(int(ferro_phase_num)+1)[1:]:
            try:
                add_perturbation("refer-POSCAR" + str(index) + '/CONTCAR',Q2D=Q2D,FQFE=FQFE,loop_index=loop_index)
            except(IndexError):
                print('the file might be empty, will jump this run, this run is refer-POSCAR'+str(index))

elif sys.argv[1]  == './EIGENVAL':
        def analyze_eigenval(eigenval_path):
            occupations_info = []
            not_metal = False
            with open(eigenval_path, 'r',errors='ignore') as file:
                lines = file.readlines()[7:]
                for line in lines:
                    info = line.strip().split()
                    if len(info) == 3:
                        occupation = np.float64(info[2])
                        occupations_info.append(round(occupation, 3))
                    elif len(info) == 5:
                        occupation_spin_up = np.float64(info[3])
                        occupation_spin_down = np.float64(info[4])
                        occupations_info.append(occupation_spin_up)
                        occupations_info.append(occupation_spin_down)
                        continue
            for occupation in occupations_info:
                if occupation == 1.00 or occupation == 0.00:
                    not_metal = True
                else:
                    not_metal = False
                    break

            if not_metal:
                print("not_metal")
            else:
                print("metal")

        analyze_eigenval('./EIGENVAL')

elif sys.argv[1]  == 'process_data':
    # 2D judge
    Q2D = False
    regex_pattern_Q2D = r'^\s*Q2D\s*=\s*\.\s*true\s*'
    pattern_Q2D = re.compile(regex_pattern_Q2D, re.IGNORECASE)
    file_path = '../test.sh'
    with open(file_path, 'r') as file:
        Q2D = any(pattern_Q2D.match(line) for line in file)


    def generate_equal_atoms_3D(atom_info):
        sub_atom_info=[0.00,0.00,0.00]
        for i in range(len(atom_info)):
            if atom_info[i] == 0.00:
                sub_atom_info[i]=0.50
            else:
                sub_atom_info[i]=atom_info[i]
        return sub_atom_info

    def cal_metal_center_info(line):
        structure = Structure.from_file(r'./' + line + '/POSCAR')
        # check mass center
        POSCAR = read_poscar.read_poscar(r'./' + line + '/POSCAR')
        atoms_info_list = POSCAR['atoms_axis']
        center_point_per_elements = []
        for atoms in atoms_info_list.values():
            atoms_equaled = []
            for atom in atoms:
                if Q2D:
                    atom = generate_equal_atoms(atom)
                else:
                    atom = generate_equal_atoms_3D(atom)
                for i in range(len(atom)):
                    if atom[i] < 0:
                        atom[i] += 1
                atoms_equaled.append(atom)
            center_point_per_elements.append(np.mean(atoms_equaled, axis=0))
        center = np.mean(center_point_per_elements, axis=0)
        analyzer = SpacegroupAnalyzer(structure)
        space_group = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()
        crystal_system = analyzer.get_crystal_system()
        return {'crystal_type':'metal','mass_center': [round(center[0],2),round(center[1],2),round(center[2],2)], 'space_group': space_group,
                'space_group_number': space_group_number, 'crystal_system': crystal_system}

    def process_polar(polar_cal_value, line):
        def cal_pro_value_double(value,lattice):
            divided_value = value / lattice
            fraction_part,_ = math.modf(divided_value)
            if fraction_part < 0.00:
                upper = fraction_part + 1.00
                lower = fraction_part
            elif fraction_part > 0.00 :
                upper = fraction_part
                lower = fraction_part - 1.00
            else:
                upper , lower = 0,0
            return upper, lower
        p_elec = np.array([float(value) for value in polar_cal_value[0].split()])
        p_ion = np.array([float(value) for value in polar_cal_value[1].split()])
        polar_value = p_elec + p_ion
        # obtain the lattic parameters for polar quantum
        structure = Structure.from_file(r'./' + line + '/POSCAR')
        a = round(structure.lattice.a,2)
        b = round(structure.lattice.b,2)
        c = round(structure.lattice.c,2)
        polar_value_a = round(polar_value[0],2)
        polar_value_b = round(polar_value[1],2)
        polar_value_c = round(polar_value[2],2)
        polar_proportion_double_a = cal_pro_value_double(polar_value_a,a)
        polar_proportion_double_b = cal_pro_value_double(polar_value_b,b)
        polar_proportion_double_c = cal_pro_value_double(polar_value_c,c)

        # check mass center
        POSCAR = read_poscar.read_poscar(r'./' + line + '/POSCAR')
        atoms_info_list = POSCAR['atoms_axis']
        center_point_per_elements = []
        for atoms in atoms_info_list.values():
            atoms_equaled = []
            for atom in atoms:
                atom = generate_equal_atoms(atom)
                for i in range(len(atom)):
                    if atom[i] < 0:
                        atom[i] += 1
                atoms_equaled.append(atom)
            center_point_per_elements.append(np.mean(atoms_equaled, axis=0))
        center = np.mean(center_point_per_elements,axis=0)

        #judge polar direction
        if center[0] > 0.50:
            polar_a = round(polar_proportion_double_a[0] * a , 2)
        elif center[0] < 0.50:
            polar_a = round(polar_proportion_double_a[1] * a , 2)
        else:
            polar_a = 0.00

        if center[1] > 0.50:
            polar_b = round(polar_proportion_double_b[0] * b, 2)
        elif center[1] < 0.50:
            polar_b = round(polar_proportion_double_b[1] * b, 2)
        else:
            polar_b = 0.00

        if Q2D:
            max_value = float('-inf')
            min_value = float('inf')
            for values in atoms_info_list.values():
                for sub_value in values:
                    z_value = sub_value[2]
                    if z_value > max_value:
                        max_value = z_value
                    if z_value < min_value:
                        min_value = z_value
            mid = round(max_value-min_value,2)
            if center[2] > mid:
                polar_c = round(polar_proportion_double_c[0] * c, 2)
            elif center[2] < mid:
                polar_c = round(polar_proportion_double_c[1] * c, 2)
            else:
                polar_c = 0.00
        else:
            if center[2] > 0.50:
                polar_c = round(polar_proportion_double_c[0] * c, 2)
            elif center[2] < 0.50:
                polar_c = round(polar_proportion_double_c[1] * c, 2)
            else:
                polar_c = 0.00

        #get crystal information
        analyzer = SpacegroupAnalyzer(structure)
        space_group = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()
        crystal_system = analyzer.get_crystal_system()
        return {'polar_value':[polar_a,polar_b,polar_c],'space_group':space_group,'space_group_number':space_group_number,
                'crystal_system':crystal_system,'lattice_parameters':[round(a,2),round(b,2),round(c,2)]}

    def read_data(file_path, file_type):

        data = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if file_type == 'refer' and line.startswith('refer-POSCAR'):
                index = line.split('refer-POSCAR')[-1].strip()
                next_line = lines[i + 1].strip()
                if next_line.startswith("Sorry"):
                    data[index] = [next_line]
                    i += 2
                else:
                    # Assumes next two lines are the data lines
                    data[index] = [lines[i + 1].strip(), lines[i + 2].strip()]
                    i += 3
            elif file_type == 'target' and line.startswith('POSCAR-'):
                index = line.split('POSCAR-')[-1].strip()
                try:
                    next_line = lines[i + 1].strip()
                except(IndexError):
                    break
                if next_line.startswith("Sorry"):
                    data[index] = [next_line]
                    data[index] = cal_metal_center_info(line)
                    i += 2
                elif next_line.startswith('POSCAR'):
                    i += 1
                else:
                    # Assumes next two lines are the data lines
                    data[index] = [lines[i + 1].strip(), lines[i + 2].strip()]
                    data[index] = process_polar(data[index], line)
                    i += 3                    
            else:
                i += 1
        return data


    def merge_data(refer_data, target_data, energy_data):
        merged_data = {}

        def get_base_key(key):
            return key.split('-')[0] if '-' in key else key

        #deal energy data
        filtered_energy_data = {k: round(float(v),3) for k, v in energy_data.items() if v != ''}
        sorted_energy_data = sorted(filtered_energy_data.items(), key=lambda item: item[1])
        ranked_energy_data = {k: (v, f"rank_{i + 1}") for i, (k, v) in enumerate(sorted_energy_data)}

        # Merge energy data with target_data
        for key,value in ranked_energy_data.items():
            if key in target_data:
                try:
                    target_data[key]['energy'] = value
                except TypeError:
                    target_data[key].append(value)
            else:
                continue
        del ranked_energy_data,sorted_energy_data,filtered_energy_data
        # Create a dictionary to group keys by their base number
        base_refer_keys = {}
        base_target_keys = {}
        for key in refer_data.keys():
            base_key = get_base_key(key)
            if base_key not in base_refer_keys:
                base_refer_keys[base_key] = []
            base_refer_keys[base_key].append(key)

        for key in target_data.keys():
            base_key = get_base_key(key)
            if base_key not in base_target_keys:
                base_target_keys[base_key] = []
            base_target_keys[base_key].append(key)

        # Merge the dictionaries based on the grouped keys
        all_base_keys = set(base_refer_keys.keys()).union(base_target_keys.keys())

        for base_key in all_base_keys:
            merged_data[base_key] = {
                'Refer': {k: refer_data[k] for k in base_refer_keys.get(base_key, [])},
                'Target': [{k: target_data[k]} for k in base_target_keys.get(base_key, [])]
            }
        return merged_data,target_data


    def write_data(merged_data,target_merged_data, output_path,res_target_energy_path):
        with open(output_path, 'w') as file:
            for key in sorted(merged_data.keys(), key=int):
                file.write(f"POSCAR {key}:\n")
                file.write("Refer Data:\n")
                for subkey, line in merged_data[key]['Refer'].items():
                    file.write(f"{subkey}: {line}\n")
                file.write("Target Data:\n")
                for item in merged_data[key]['Target']:
                    for subkey, line in item.items():
                        if len(line) == 1:
                            continue
                        else:
                            file.write(f"{subkey}: {line}\n")
                file.write("\n")

        with open(res_target_energy_path, 'w',encoding= 'utf-8') as file:
            sorted_res = sorted(target_merged_data.items(), key=lambda item: item[1]['energy'][0])
            sorted_res_re = [{'id':i[0],'info':i[1]} for i in sorted_res]
            for data in sorted_res_re:
                file.write(f"POSCAR {data['id']} energy per atoms {data['info']['energy'][0]} {data['info']['energy'][1]}:\n")
                if 'polar_value' in data['info'].keys():
                    file.write(f"Polar_value {data['info']['polar_value'][0]} {data['info']['polar_value'][1]} {data['info']['polar_value'][2]} unit:e*Å \n")
                else:
                    file.write(f"Crystal type: {data['info']['crystal_type']} \n")
                    file.write(f"Crystal center: {data['info']['mass_center']} \n")
                file.write("crystal info:\n")
                file.write(f"crystal_system: {data['info']['crystal_system']} \n")
                file.write(f"space_group: {data['info']['space_group']} \n")
                file.write(f"spacegroup number: {data['info']['space_group_number']} \n")
                try:
                    file.write(f"lattice_parameters: {data['info']['lattice_parameters'][0]} {data['info']['lattice_parameters'][1]} "
                               f"{data['info']['lattice_parameters'][2]}\n")
                except KeyError:
                    pass
                file.write(f"--------------------------------------------------------------\n")

    def extract_valid_target_data(input_file, output_file):
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            lines = infile.readlines()
            index = None
            capturing = False
            buffer = []

            for line in lines:
                line = line.strip()
                if line.startswith('POSCAR'):
                    index = line.split('POSCAR')[1].strip(':')
                    # If there's content to save and the last section was about target data with values
                    if buffer and 'Sorry, seems' not in buffer[-1]:
                        outfile.write('\n'.join(buffer) + '\n\n')
                    buffer = [f"POSCAR{index}:"]  # Start a new buffer for a new section
                elif line == "Refer Data:" or line == "Target Data:":
                    if buffer and 'Target Data:' in buffer:
                        # Check if previous Target Data section should be saved
                        if 'Sorry, seems' not in buffer[-1]:
                            outfile.write('\n'.join(buffer) + '\n\n')
                        buffer = buffer[:1]  # Keep the POSCAR header
                    buffer.append(line)
                    capturing = True
                elif capturing:
                    if line:
                        buffer.append(line)
                    else:
                        capturing = False  # Stop capturing if empty line

            # Write the last buffered content if valid
            if buffer and 'Sorry, seems' not in buffer[-1]:
                outfile.write('\n'.join(buffer) + '\n')
    def extract_poscar_files(indices_file, source_dir, target_dir):
        # Ensure the target directory exists
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        # Read indices from the res.info file, ensuring that each line is valid
        refer_indices = []
        target_indices = []

        with open(indices_file, 'r') as file:
            lines = file.readlines()

            current_poscar = None
            in_target_data = False

            for line in lines:
                stripped_line = line.strip()

                if stripped_line.startswith("POSCAR"):
                    current_poscar = stripped_line.split()[1].replace(":", "")
                    refer_indices.append(current_poscar)
                    in_target_data = False

                elif stripped_line.startswith("Target Data:"):
                    in_target_data = True

                elif in_target_data and ':' in stripped_line:
                    target_key = stripped_line.split(':')[0]
                    if target_key != "Energy":
                        target_indices.append(target_key)

        # Iterate over each index, copy and rename POSCAR files
        #for index in refer_indices:
        #    refer_poscar_path = os.path.join(source_dir, f"refer-POSCAR{index}", "POSCAR")
        #    if os.path.exists(refer_poscar_path):
        #        shutil.copy(refer_poscar_path, os.path.join(target_dir, f"POSCAR-refer-{index}"))

        for index in target_indices:
            target_poscar_path = os.path.join(source_dir, f"POSCAR-{index}", "POSCAR")
            if os.path.exists(target_poscar_path):
                shutil.copy(target_poscar_path, os.path.join(target_dir, f"POSCAR-target-{index}"))


    def read_energy_data(file_path):
        energy_data = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('POSCAR-'):
                index = line.split('POSCAR-')[-1].strip()
                try:
                    energy_value = lines[i + 1].strip()
                except(IndexError):
                    break
                energy_data[index] = energy_value
                if energy_value == '':
                    energy_value = float('inf')
                i += 2
            else:
                i += 1
        return energy_data

    refer_path = r'./refer_polar_info.dat'
    target_path = r'./target_polar_info.dat'
    output_path = r'./search_result.info'
    energy_path = r'./target_polar_energy.dat'
    res_path = r'./res.info'
    res_target_energy_path = r'./results.info'
    res_dir = r'./ferro_search_results'

    refer_data = read_data(refer_path, "refer")
    target_data = read_data(target_path, "target")
    energy_data = read_energy_data(energy_path)

    merged_data,target_merged_data = merge_data(refer_data, target_data,energy_data)
    write_data(merged_data,target_merged_data, output_path,res_target_energy_path)

    extract_valid_target_data(output_path, res_path)
    extract_poscar_files(res_path, './', res_dir)
    #shutil.move(res_path, res_dir)
    shutil.move(res_target_energy_path, res_dir)

    # os.remove(refer_path)
    # os.remove(target_path)

else:
    print("something seems to be wrong,check your scripts")
