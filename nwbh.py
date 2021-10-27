'''
This program does 4 things:
1: Analyzes the previous basin hops, if any
2: Perturbs the coordinates of the previous hop. If this is impossible, it 
   generates new coordinates
3: Initializes a new calculation
4. Detects if the search is finished, possibly performing an action if it is
'''
import random, math, csv, sys, copy, subprocess, re
import nwchem_interface
import numpy as np

def outname(name, num_calc, num_hop):
    return name + '_' + str(num_calc) + '_' + str(num_hop) + '.out'
def inpname(name, num_calc, num_hop):
    return name + '_' + str(num_calc) + '_' + str(num_hop) + '.inp'


def get_rand_bond_length(atom1, atom2, radius_dict):
    '''
    Generates a uniformly random bond length between the sum of the minimum 
    atomic radii and the maximum atomic radii, both of which are given by the 
    radius dictionary
    
    atom1, atom2: strings
        The atomic symbols for the atoms between which the bond length is found
    radius_dict: dict
        The radius dictionary. Key is atomic symbol, element is tuple of 
        (<minimum radius>, <maximum radius>)
    '''
    return random.uniform(radius_dict[atom1][0] + radius_dict[atom2][0], radius_dict[atom1][1]+radius_dict[atom2][1])

def append_atom(coordinates, atom_name, radius_dict):
    '''
    Appends a new atom using a semiempirical semistochastic algorithm. 
    By repeating this, we can grow a semiplausible randomly generated cluster 
    to initialize our global minimum search. This is preferable to a truly 
    unbiased stochastic generation, since it avoids fragentation and lessens 
    the chance of a calculation failing.
    
    coordinates: list of lists
        Each element has the form 
        [<atomic symbol (string)>, <x (float)>, <y (float)>, <z (float)>].
    atom_name: string
        The name of the atom to be appended.
    radius_dict: dict
        The radius dictionary. Key is atomic symbol, element is tuple of 
        (<minimum radius>, <maximum radius>)
    '''
    if len(coordinates)==0:
        coordinates.append([atom_name, 0.0, 0.0, 0.0])
    elif len(coordinates)== 1:
        coordinates.append([atom_name, 0.0, 0.0, get_rand_bond_length(atom_name, coordinates[0][0], radius_dict)])
    elif len(coordinates) > 1:
        failures = 0
        max_appending_failures = 100
        while failures <= max_appending_failures:
            appendee = random.randint(0, len(coordinates)-1)
            bond_length = get_rand_bond_length(atom_name, coordinates[appendee][0], radius_dict)
            new_atom = np.random.randn(3)
            new_atom = new_atom*bond_length/np.linalg.norm(new_atom) + np.array(coordinates[appendee][1:4])
            x, y, z = new_atom[0], new_atom[1], new_atom[2]
            successfully_appended = True
            for atom in coordinates:
                if np.linalg.norm(atom[1:4]-new_atom) < radius_dict[atom[0]][0] + radius_dict[atom_name][0]:
                    failures += 1
                    successfully_appended = False
                    break
            if successfully_appended == True:
                coordinates.append([atom_name, x, y, z])
                break
            if failures > max_appending_failures:
                print('FAILURE')
    return coordinates

def create_radius_dict(atom_set, csv_name = 'atomic_information.csv', min_mult=1.0, max_mult=1.0):
    '''
    Returns a dictionary object with keys of atomic symbols and values of a 
    tuple of (minimum atomic radius*min_mult, maximum atomic radius*max_mult). 
    This is all done so that the global minimum search can use empirical 
    paramaters to easily rule out structures that are unrealistic or yield 
    divergent calculations due to either structural fragmentation or too small
    internuclear distances.
    atom_set: set
        The atoms we need to build this dictionary to include, whcich are the 
        list of elements in the molecule we're searching for.
    min_mult: float
        Determines the minimum bond length upon multiplication with the sum of 
        the minimum atomic radii of each atom. Smaller values mean smaller 
        possible bonds, which leads to less bias. But beware! This also raises 
        the likelihood of a calculation failing.
    max_mult: float
        Determines the maximum bond length upon multiplication with the sum of 
        the maximum atomic radii of each atom. Larger values mean larger 
        possible bonds, which leads to less bias. But beware! This also raises 
        the likelihood of fragmentation ruining the search.
    '''
    radius_dict = dict({})
    with open(csv_name) as radius_csv:
        reader = csv.reader(radius_csv)
        for row in reader:
            if row[1] in atom_set:
                radius_dict.update({row[1] : ((min_mult*float(row[3]), max_mult*float(row[2])))})
    print(radius_dict)
    return radius_dict

def gencoords(atom_list, radius_dict):
    '''
    Parameters
    ----------
    atom_list : list of strings
        A list of all atoms in the cluster/molecule.
    radius_dict : dictionary
        The radius dictionary. Key is atomic symbol, element is tuple of 
        (<minimum radius>, <maximum radius>)

    Returns
    -------
    coordinates : List of lists
        Each element has the form 
        [<atomic symbol (string)>, <x (float)>, <y (float)>, <z (float)>].
    '''
    alist = copy.copy(atom_list)
    random.shuffle(alist)
    coordinates = []
    for atom in alist:
        coordinates = append_atom(coordinates, atom, radius_dict)
    return coordinates


def far_enough(atom, x, y, z, c2, radius_dict):
    '''
    Checks to see if two atoms are closer together than they should be, which
    is the minimum bond length for the two atoms (determined by their atomic 
    radii times a paramater that lets us lower the minimum bond length)
    '''
    if (x-c2[1])**2 + (y-c2[2])**2 + (z-c2[3])**2 > radius_dict[atom][0]+radius_dict[c2[0]][0]:
        return True
    return False

def hopcoords(coords, max_p, radius_dict, max_errors = 10000):
    '''
    This method perturbs each coordinate of each atom by a random amount. This 
    is done by adding a uniformly random number between + and - max_p to each 
    coordinate of each atom. The perutbartion of each atom is accepted if it is
    a minimum distance from all other atoms. Note that this is done sequentially,
    with each atom pertubed one-at-a-time, rather than all of them hopping at 
    once.
    I'll try a comparison of this sequential hop with a simeltaneous hop, when 
    I do a study on B20. I wonder which is better, though I suspect it's the 
    simeltaneous hop. That is more true to the spirit of Wales and Doye...
    '''
    for coord in coords:
        for n in range(0, max_errors):
            rx = coord[1] + random.uniform(-max_p, max_p)
            ry = coord[2] + random.uniform(-max_p, max_p)
            rz = coord[3] + random.uniform(-max_p, max_p)
            acceptable = True
            for coord2 in coords:
                if coord2 != coord:
                    if not far_enough(coord[0], rx, ry, rz, coord2, radius_dict):
                        acceptable = False
                        break
            if acceptable:
                coord[1] = rx
                coord[2] = ry
                coord[3] = rz
                break
    return coords


def write_data(name, output, num_calc, num_hop, kT, max_p, accept):
    '''
    This writes a line in _data.txt summarizing the results of a calculation. 
    This line consists of the file name, energy, acceptance, temperature, and 
    maximum perturbation.
    accept: string
        Either 'T' or 'F' that says if the last hop was accepted
    '''
    with open(name + '/_data.txt', 'a') as data:
        data.write(outname(name, num_calc, num_hop-1) + ';Energy=' + str(output[0]))
        data.write(';Accept=' + accept + ';kT=' + str(kT) + ';max_p=' + str(max_p) + '\n')


def hop(name, num_calc, num_hop, atom_list, radius_dict, data=True, max_p=0.5, kT=0.000001):
    '''
    This compares the two most recent hop calculations. It will accept the most
    recent if either it is lower in energy than the second most recent or 
    math.exp((output2[0]-output1[0])/kT is greater than a random number between
    0 and 1. If kT=0, then the latter never happens.
    Additionally, this will write the result of the previous hop calculation in
    the _data.txt file of the calculation, using the write_data method.
    kT: float
        Temperature of the metropolis acceptance criterion. kT=0 is monotonic
    max_p:
        maximum perturbation experienced by each of the coordinates. Larger 
        values cause more deformation, but smaller values may not allow the 
        structure to hop basins. This value will be tailored in the future.
    '''
    if num_hop == 0:
        coords = gencoords(atom_list, radius_dict)
    elif num_hop == 1:
        output1 = nwchem_interface.readopt(outname(name + '/' + name, num_calc, num_hop-1))
        if output1[1] == 'nc':
            coords = gencoords(atom_list, radius_dict)
            write_data(name, output1, num_calc, num_hop, kT, max_p, 'F')
        else:
            coords = hopcoords(output1[1], max_p, radius_dict)
            write_data(name, output1, num_calc, num_hop, kT, max_p, 'T')
    else:
        output1 = nwchem_interface.readopt(outname(name + '/' + name, num_calc, num_hop-1))
        output2 = nwchem_interface.readopt(outname(name + '/' + name, num_calc, num_hop-2))
        if output1[1] == 'nc':
            write_data(name, output1, num_calc, num_hop, kT, max_p, 'F')
            if output2[1] == 'nc':
                coords = gencoords(atom_list, radius_dict)
            else:
                coords = hopcoords(output2[1], max_p, radius_dict)
 #           coords = gencoords(atom_list, radius_dict)
        elif output1[0] < output2[0] or random.uniform(0,1) > math.exp((output2[0]-output1[0])/kT):
            coords = hopcoords(output1[1], max_p, radius_dict)
            write_data(name, output1, num_calc, num_hop, kT, max_p, 'T')
        else:
            coords = hopcoords(output2[1], max_p, radius_dict)
            write_data(name, output1, num_calc, num_hop, kT, max_p, 'F')
    for c in coords:
        c[1], c[2], c[3] = round(c[1],5),round(c[2],5),round(c[3],5)
    return coords

def calculate(name, num_calc, num_hop, max_hops, max_calcs, multiplicity):
    '''
    Submits an nwchem6.8 calculation on a linux cluster with the SLURM workload
    manager. A calculation for the next basin hop (nwbh.py) is also submitted 
    with a dependency on the nwchem calculation such that it only begins after 
    the nwchem calculation has completed.
    '''
    j1 = subprocess.Popen('sbatch nwsubmit.sh ' + name + '/' + name + '_' + str(num_calc) + '_' + str(num_hop), shell = True, stdout=subprocess.PIPE)
    jid1 = ''.join(list(filter(str.isdigit, j1.communicate()[0].decode('utf-8'))))
    subprocess.call('sbatch --dependency=afterany:'+str(jid1)+' nwbh.sh '+str(num_calc)+' '+str(num_hop + 1)+' '+str(max_hops)+' '+str(max_calcs)+' '+str(multiplicity),shell=True)

def get_atom_list(name):
    '''
    Finds the list of atoms from an atomic formula string. The string formula 
    should have all elemental symbols followed by a number (even if there is 
    only one atom of that element), the formula should terminate with an
    underscore, and all elements should have proper capitalization.
    Good formulas: C2H6_t, C1H4_methane, C1O1_monoxide
    Bad formulas: C2H6, CH4_methane, c1o1_monoxide
    
    For example, the formula 'Re2B3_test' gives the list 
    ['Re', 'Re', 'B', 'B', 'B'] 
    
    The flexibility of this will be improved soon, I promise.
    '''
    l = 0
    word_list, atom_list = [], []
    for k in range(1, len(name)):
        if name[k].isupper() or name[k] == '_':
            word_list.append(name[l:k])
            l=k
    for w in word_list:
        for k in range(0, len(w)):
            if w[k].islower():
                if w[k+1:] == '':
                    atom_list.append(w)
                else:
                    for m in range(0, int(w[k+1:])):
                        atom_list.append(w[0:k+1])
                break
            elif w[k].isdigit():
                for m in range(0,int(w[k:])):
                    atom_list.append(w[0:k])
                break
            elif len(w) == 1:
                atom_list.append(w)
    return atom_list

search_specs = 'search_specs.txt'
min_dist_multiplier = 0.8

with open(search_specs) as specs:
    for line in specs:
        if line.find('title') >= 0:
            t = re.search('\"(.+?)\"',line)
            if t:
                name = t.group(1)
            break

num_calc = int(sys.argv[1])
num_hop = int(sys.argv[2])
max_hops = int(sys.argv[3])
max_calcs = int(sys.argv[4])
multiplicity = int(sys.argv[5])
name = name + str(multiplicity)

atom_list = get_atom_list(name)

seed = num_calc*1000000 + num_hop
random.seed(seed)
np.random.seed(seed)

radius_dict = create_radius_dict(set(atom_list), min_mult=min_dist_multiplier)
coords = hop(name, num_calc, num_hop, atom_list, radius_dict)
print(num_hop)
print(max_hops)
if num_hop < max_hops:
    nwchem_interface.makeinp(inpname(name+'/'+name, num_calc, num_hop), coords, 'search_specs.txt', multiplicity)
    calculate(name, num_calc, num_hop, max_hops, max_calcs, multiplicity)
else:
    done= False
    with open(name + '/_data.txt') as data:
        if sum(1 for line in data) >= max_calcs*max_hops:
            done = True
    if done == True:
        print('sbatch nwsifter.sh ' + name)
        j = subprocess.call('sbatch nwsifter.sh ' + name + ' ' + str(multiplicity), shell = True, stdout=subprocess.PIPE)
        with open(name + '/_data.txt', 'a') as data:
            data.write('DONE!!\n')
