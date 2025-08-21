# from part import *
# from material import *
# from section import *
# from assembly import *
# from step import *
# from interaction import *
# from load import *
# from mesh import *
# from optimization import *
# from job import *
# from sketch import *
# from visualization import *
# from connectorBehavior import *
from abaqus import *
from abaqusConstants import *
import numpy as np
import os
from PostProcessing import *
from ZhangEtAlValidation import *
import time

# mdb.close()
yield_stress_sweep = [
    200.0,
    300.0,
    500.0,
    600.0,
    800.0,
    900.0,
    1100.0,
    1200.0,
    1400.0,
    1500.0,
    1700.0,
    1800.0,
    2000.0,
]
strain_hardening_sweep = np.arange(0.1, 0.6, 0.2)
depth = -50e-3
friction_coefficient = (
    0.0  # does nothing currently as tangential contact is frictionless
)
E_modulus = 200000.0
density = 7.8e-9
poisson = 0.3
jobName = "ExplicitRigidIndenterScratch"
for sigma_y in yield_stress_sweep:
    for n in strain_hardening_sweep:
        print(sigma_y, n)
        # Call model and run model
        ScratchModelParameterSweep(
            jobName=jobName,
            sigma_y=sigma_y,
            strain_hardening_index=n,
            depth=depth,
            friction_coefficient=friction_coefficient,
            E_modulus=E_modulus,
            density=density,
            poisson=poisson,
        )
        # Do post processing
        PostProcess(
            sigma_y=sigma_y,
            strain_hardening_index=n,
            depth=depth,
            friction_coefficient=friction_coefficient,
            E_modulus=E_modulus,
            density=density,
            poisson=poisson,
        )
