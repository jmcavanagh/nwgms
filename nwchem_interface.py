'''
NWChem interface methods. 
Handles reading output files and writing input files for nwchem6.8-openmpi.
'''


'''
Reading output files
'''

def readopt(outputfile):
    '''
    Reads an nwchem optimization output file. This finds the lowest energy 
    found and the set of atomic coordinates that correspond to this energy in 
    tuple form. If the calculation fails to converge, it returns the tuple 
    (0, 'nc').

    Parameters
    ----------
    outputfile : string
        The path of the nwchem output file

    Returns
    -------
    energy : float
        The lowest energy found in the nwchem output file
    coordinates / no coordinates : list of lists / string
        atomic coordinates of the cluster, if found. If not, the string 'nc'.
        coordinates are all in angstroms, and their form is
        [<atomic symbol (string)>, <x (float)>, <y (float)>, <z (float)>].

    '''
    with open(outputfile) as out:
        energy = 0
        coords = []
        writecoords=False
        while True:
            line = out.readline()
            if line.find('@') >= 0 and line.find('--') < 0 and line.find('Step') < 0:
                if float(line.split()[2]) < energy:
                    energy = float(line.split()[2])
                    writecoords = True
                    coords = []
            elif writecoords == True and line.find('Output coordinates in angstroms') >= 0:
                line = out.readline()
                line = out.readline()
                line = out.readline()
                while True:
                    line = out.readline()
                    raw = line.split()
                    if raw == []:
                        break
                    coords.append([raw[1], float(raw[3]), float(raw[4]), float(raw[5])])
                    writecoords = False
            elif not line:
                if coords == []:
                    return (energy, 'nc')
                return (energy, coords)







'''
Writing input files
'''

def makeinp(inputfile, coords, template, multiplicity):
    '''
    Writes an input file for nwchem 6.8 that performs the desired calculation
    as specified by the template, on the atomic coordinates specified
    in coords.

    Parameters
    ----------
    inputfile : string
        The path/name of the nwchem input file to be created
    coords : list of lists
        Atomic coordinates of the cluster, in angstroms. Their form is
        [<atomic symbol (string)>, <x (float)>, <y (float)>, <z (float)>]
    template : string
        The path/name of the specifications of the calculation. This template 
        file contains all of the information (basis set, approximation method, 
        multiplicity, etc.) needed other than the atomic coordinates.

    Returns
    -------
    None. Creates a file.
    '''
    with open(inputfile, 'w+') as inp:
        inp.write('geometry\n')
        for c in coords:
            inp.write(c[0] + ' ' + str(c[1]) + ' ' + str(c[2]) + ' ' + str(c[3]) + '\n')
        inp.write('end\n')
        with open('search_specs.txt') as template:
            for t in template:
                if t.find('title') < 0:
                    inp.write(t)
                    if t.find('xc') >= 0:
                        print("found xc")
                        inp.write(' mult ' + str(multiplicity) + '\n')
    return


                
                
            

