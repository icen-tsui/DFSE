import numpy as np
class POSCAR:
    def __init__(self,path):
       self.POSCAR_info = read_poscar(path) 

def read_poscar(path):
    with open(path) as f:
            poscar = f.readlines()
            ##get primitive cell lattice information
            scale_factor = float(deal_space(poscar[1])[0])
            lattice_matrix = get_lattice_matrix(poscar[2:5])
            element_consist = deal_space(poscar[5])
            nums_of_element = deal_space(poscar[6])
            element_info = dict(zip(element_consist,[int(x) for x in nums_of_element]))
            if 'Direct' in poscar[7] or  'Cartesian' in poscar[7]:
                Coordinate_system_type = deal_space(poscar[7])
                atom_axis_info = poscar[7:8+sum(map(int,element_info.values()))]
            elif 'Direct' in poscar[8] or  'Cartesian' in poscar[8] and (poscar[7][0]=='s'or'S'):
                Coordinate_system_type = deal_space(poscar[8])[0]
                atom_axis_info = poscar[8:9+sum(map(int,element_info.values()))]
            axis_type=atom_axis_info[0]
            atoms_axis=get_atoms_axis_list(atom_axis_info[1:],element_consist,[int(num) for num in nums_of_element])
            f.close()
            data = {
                "scale_factor": scale_factor,
                "lattice_matrix": lattice_matrix,
                "Coordinate_system_type":Coordinate_system_type,
                "element_info": element_info,
                "atoms_axis": atoms_axis,
            }

            return data

def get_lattice_matrix(matric_info:list):
    #get a,b,c lattice_matrix
    res=[]
    for line in matric_info:
        line = line.lstrip(' ').rstrip('\n').split()
        line = [float(x) for x in line]
        res.append(line)
    res = np.array(res).astype(float)
    return res


def deal_space(content:str):
    #deal_space and delete 'T' 'F'
    content=content.strip(' ').rstrip('\n').split(' ')
    content=[ x for x in content if x!=' 'and x!='T' and x!='F' and x!='' ]
    return content


def get_atoms_axis_list(atoms_info:list,atoms_kind:list,atoms_nums:list):
    #take accurate atom_axis information 
    #eg.{'P4':array[[a_lattice,b_lattice,c_lattice],**]}
    start_point = 0
    end_point = 0
    res={}
    atoms_info=[deal_space(atom_info) for atom_info in atoms_info]
    for i in range(len(atoms_kind)):
        end_point += atoms_nums[i]
        key = tuple((atoms_kind[i],atoms_nums[i]))
        res[key] = np.array(atoms_info[start_point:end_point]).astype(float)
        start_point = end_point
    return res        





    