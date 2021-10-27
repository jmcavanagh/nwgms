import sys
import numpy as np
import nwchem_interface as nw
import matplotlib.pyplot as plt
from pyscf.data import nist
from pyscf import gto, scf, dft, tddft
from pyscf.geomopt.berny_solver import optimize


def anionic_dft(mol, xc, anion_mult):
    mol.charge = -1
    mol.spin = anion_mult - 1
    mf = dft.UKS(mol)
    mf.xc = xc
    energy = mf.kernel()
    return energy

def neutral_tddft(mol, xc, neutral_mult):
    mol.charge = 0
    mol.spin = neutral_mult - 1
    mf = dft.UKS(mol)
    mf.xc = xc
    gs = mf.kernel()
    mytd = tddft.TDDFT(mf)
    mytd.kernel(nstates=3)
    mytd.analyze()
    e_ev = np.asarray(mytd.e) * nist.HARTREE2EV
    return gs, e_ev, mytd


def build_mol(output_file, basis_set, xc, anion_mult):
    energy, coord_list = nw.readopt(output_file)
    print(energy)
    coord_string = '; '.join(' '.join([str(i) for i in c]) for c in coord_list)
    print(coord_string)
    mol = gto.Mole()
    mol.build(
        atom = coord_string,
        basis = basis_set,
        symmetry = True,
        charge = -1,
        spin = anion_mult - 1,
        verbose = 0,
        ecp = basis_set
    )
    return mol
    
try:
    nwoutput = sys.argv[1]
    basis_set = sys.argv[2]
    xc = sys.argv[3]
    anion_mult = int(sys.argv[4])
    neutral_mult = int(sys.argv[5])
    mol = build_mol(nwoutput, basis_set, xc, anion_mult)
    ags = anionic_dft(mol, xc, anion_mult)
    ngs, e_ev, mytd = neutral_tddft(mol, xc, neutral_mult)
    
    vde1 = (ngs - ags)*nist.HARTREE2EV
    e_ev = np.append(e_ev,0)
    e_ev += vde1
    
    x = np.arange(1, 7, 0.01)
    y = np.zeros_like(x)
    for ex in e_ev:
        y += np.exp(-40*(x-ex)**2)
    plt.plot(x,y, color= 'b')
    
    y_single = np.exp(-40*(x-e_ev[0])**2)
    xy = mytd.xy
    for k in range(0, len(xy)):
        exc_mat = abs(xy[k][0])
        orbital_transition = np.where(exc_mat == np.max(exc_mat))
        orbital_transition = [int(h) + 1 for h in orbital_transition]
        if orbital_transition[1] == 0:
            y_single += np.exp(-40*(x-e_ev[k+1])**2)
        else:
            print('MULTIELECTRON EXCITATION: ' + str(k))
    plt.plot(x, y_single, color='r')
    plt.savefig('B3_spec.png')
except:
    print("sbatch testspec.sh nwoutput basis_set xc anion_mult neutral_mult")
    
    
    
    
    
    
