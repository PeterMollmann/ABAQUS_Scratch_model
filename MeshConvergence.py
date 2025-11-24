from abaqus import *
from abaqusConstants import *
from ProgressiveLoadScratch.PostProcessing import *
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from ProgressiveLoadScratch.helpers import run_job_and_wait

# from material_parameters import parameters
import numpy as np
import shutil
from cleanup import cleanupAbaqusJunk
import ProgressiveLoadScratch.Constants as C


### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "MeshConvergence"

# Array of shape (#mesh combinations, 3). Rows are mesh combinations, coloums are [meshSizeX, meshSizeY, meshSizeZ],
meshSizes = np.array(
    [
        # [0.040, 0.040, 0.040],
        # [0.030, 0.030, 0.030],
        # [0.020, 0.020, 0.020],
        # [0.015, 0.015, 0.015],
        # [0.010, 0.010, 0.010],
        # [0.008, 0.010, 0.010],
        # [0.010, 0.010, 0.008],
        # [0.008, 0.010, 0.008],
        # [0.008, 0.008, 0.008],
        # [0.006, 0.006, 0.006],
        [0.005, 0.005, 0.005],
        [0.004, 0.004, 0.004],
    ]
)


# Parallelisation
num_cpus = 6
num_domains = num_cpus

# Change abaqus working directory
rundir = os.path.join("runs", "MeshConvergence")
if not os.path.exists(rundir):
    os.makedirs(rundir)
os.chdir(rundir)


rho = 8.0e-09
E = 200000.0
nu = 0.3
A = 700.0
B = 700.0
n = 0.5
D1 = 1.0
D2 = 0.5
D3 = -0.25
uts = 1000.0
kc = 2000
mu = 0.1
kappa = 1e-4

# make to dict
material_params = {
    "rho": rho,
    "E": E,
    "nu": nu,
    "A": A,
    "B": B,
    "n": n,
    "D1": D1,
    "D2": D2,
    "D3": D3,
    "uts": uts,
    "kc": kc,
    "mu": mu,
    "kappa": kappa,
}

include_wear = False

for id, meshSize in enumerate(meshSizes):

    fileName = (
        "new_mesh" + str(meshSize[0]) + "_" + str(meshSize[1]) + "_" + str(meshSize[2])
    )

    ScratchModel, SubstratePart = ScratchModelSetup(
        SubstrateSizeX=meshSize[1],
        SubstrateSizeY=meshSize[0],
        SubstrateSizeZ=meshSize[2],
        mass_scale=10e4,
        use_ALE=True,
        include_wear=include_wear,
    )

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        rho=rho,
        youngs_modulus=E,
        poisson_ratio=nu,
    )

    material.JohnsonCookHardening(A=A, B=B, n=n)
    material.JohnsonCookDamage(d1=D1, d2=D2, d3=D3)
    material.DamageEvolution(kc=kc, uts=uts, E=E, nu=nu)
    material.SectionAssignment()
    material.UpdateFrictionAndWear(mu)

    run_job_and_wait(jobName)
    PostProcess(jobName, fileName, material_params)

    mdb.close()

    sta_file = jobName + ".sta"
    target_dir = "SimDataOutputs/"
    new_name = fileName + ".sta"
    dst_file = os.path.join(target_dir, new_name)
    shutil.move(sta_file, dst_file)

    odb_file = jobName + ".odb"
    target_dir = "SimDataOutputs/"
    new_name = fileName + ".odb"
    dst_file = os.path.join(target_dir, new_name)
    shutil.move(odb_file, dst_file)

cleanupAbaqusJunk()
