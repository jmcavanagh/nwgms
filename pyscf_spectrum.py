import numpy
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
    e_ev = numpy.asarray(mytd.e) * nist.HARTREE2EV
    return gs, e_ev


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
        verbose = 5
    )
    return mol
    

nwoutput = 'B3_a1.out'
basis_set = '6-311+g'
xc = 'b3lyp'
anion_mult = 1
neutral_mult = 2

mol = build_mol(nwoutput, basis_set, xc, anion_mult)


ags = anionic_dft(mol, xc, anion_mult)
ngs, e_ev = neutral_tddft(mol, xc, neutral_mult)

vde1 = (ngs - ags)*nist.HARTREE2EV
e_ev = numpy.append(e_ev,0)
e_ev += vde1

x = numpy.arange(1, 7, 0.01)
y = numpy.zeros_like(x)
for ex in e_ev:
    y += numpy.exp(-40*(x-ex)**2)
plt.plot(x,y)




