'''
Module for parsing nwchem 6.8 output files, particularly for finding the lowest 
energy and corresponding coordinates.

This may have additional functionality in the future, but, as of now, it only 
finds the coordinates and energy.
'''


class nwoutput:
    def __init__(self, filepath):
        with open(filepath) as out:
            energy = 0
            coords = []
            counter = 0
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
                        self.coords = 'nc'
                    else:
                        self.coords = coords
                    self.energy = energy