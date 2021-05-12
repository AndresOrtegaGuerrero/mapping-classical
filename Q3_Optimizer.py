# This code is part of Qiskit Quantum Electronic Stucture Qalculator (Q3).
#
# (C) Copyright Q3-Team 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

''' Optimizer Module'''
# Import Qiskit useful libraries
#=========================================================
import sys
import os
import time
from pyscf import gto, scf, dft
from pyscf import lo
from pyscf.geomopt.geometric_solver import optimize as geometric
from pyscf.geomopt.berny_solver import optimize as berny
import numpy
#=========================================================

# Modules for the Molecule, working directory and pre-run energies

from QuEst import working_dir, electronic_set
from ..molecule import Molecule

mol = Molecule
method = electronic_set{method}
shell = electronic_set{shell}


conv_params = { # These are the default settings
    'convergence_energy': 1e-6,  # Eh
    'convergence_grms': 3e-4,    # Eh/Bohr
    'convergence_gmax': 4.5e-4,  # Eh/Bohr
    'convergence_drms': 1.2e-3,  # Angstrom
    'convergence_dmax': 1.8e-3,  # Angstrom
}

def Q3OptimizerDriver(pre_energy,mol,method,shell,**conv_params):
    if pre_energy ==True:
        opt_pre = optimize(pre_energy, **conv_params)
        opt_geom = opt_pre.tofile(working_dir+'/opt-geom.xyz')
    else:
        if method==KS:
            if shell == O:
                SPE =scf.UKS(mol).run()
            else shell == R:
                SPE =scf.RKS(mol).run()
        else method == HF:
            if shell == O:
                SPE =scf.UHF(mol).run()
            else shell == R:
                SPE =scf.RHF(mol).run()
        SPE_opt = optimize(SPE, **conv_params)
        opt_geom = mol_eq.tofile(working_dir+'opt-geom.xyz')
