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

import sys
import time
from pyscf import gto, scf
from pyscf import lo
from pyscf.tools import molden, cubegen
from pyscf.geomopt.geometric_solver import optimize
#from QuEst import mol, working_dir, 
import numpy

"""pre_energy, mol , "HF", "R", **conv_params
Args:
    pre_energy: Previous single point energy calcualtion generated previously
    mol : Molecule build by the gto module
    "HF" : HF standard theory for optimizing simple coordinates, "KS" for DFT theory
    "R" : Restricted shell of the calcualtion, "O" for open shell calcualtions
    conv_params: Specific paramenters for the optimization convergence
    
Raises:
    Q3_Error: Invalid Input
"""

conv_params = { # These are the default settings
    'convergence_energy': 1e-6,  # Eh
    'convergence_grms': 3e-4,    # Eh/Bohr
    'convergence_gmax': 4.5e-4,  # Eh/Bohr
    'convergence_drms': 1.2e-3,  # Angstrom
    'convergence_dmax': 1.8e-3,  # Angstrom
}

def Q3OptimizerDriver(pre_energy, mol , HF, R, conv_params):
    if pre_energy ==True:
        opt_pre = optimize(pre_energy, **conv_params)
        opt_geom = opt_pre.tofile(working_dir+'/opt-geom.xyz')
    else:
        if theory==KS:
            if shell == O:
                SPE =scf.UKS(mol).run()
            else shell == R:
                SPE =scf.RKS(mol).run()
        else theory == HF:
            if shell == O:
                SPE =scf.UHF(mol).run()
            else shell == R:
                SPE =scf.RHF(mol).run()

        SPE_opt = optimize(SPE, **conv_params)
        opt_geom = mol_eq.tofile(working_dir+'opt-geom.xyz')
