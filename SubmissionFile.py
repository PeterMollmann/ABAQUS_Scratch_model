from abaqus import *
from abaqusConstants import *
import ProgressiveLoadScratch.Constants as C
from ProgressiveLoadScratch.PostProcessing import PostProcess
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from material_parameters import parameters
from cleanup import cleanupAbaqusJunk
from ProgressiveLoadScratch.helpers import run_job_and_wait


### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "MaterialSweep"


meshSize = [0.040, 0.030, 0.020, 0.010, 0.008, 0.006, 0.004]
meshSizeIdx = 3


# Change abaqus working directory
rundir = os.path.join("runs", jobName)
if not os.path.exists(rundir):
    os.makedirs(rundir)
os.chdir(rundir)

include_wear = False
# Setup scratch model. Only needs to be called once
ScratchModel, SubstratePart = ScratchModelSetup(
    SubstrateSizeY=meshSize[meshSizeIdx],
    SubstrateSizeX=meshSize[meshSizeIdx],
    SubstrateSizeZ=meshSize[meshSizeIdx],
    mass_scale=10e4,
    use_ALE=False,
    include_wear=include_wear,
)

start_from_sim_id = 1  # Set to desired starting ID to skip completed simulations


for arg in parameters:
    run_id = arg["id"]
    if int(run_id) < start_from_sim_id:
        continue
    rho = float(arg["rho"])
    E = float(arg["E"])
    nu = float(arg["nu"])
    A = float(arg["A"])
    B = float(arg["B"])
    n = float(arg["n"])
    # D1 = float(arg["D1"])
    # D2 = float(arg["D2"])
    # D3 = float(arg["D3"])
    # uts = float(arg["uts"])
    # kc = float(arg["kc"])
    # u_pl_f = float(arg["u_pl_f"])
    mu = float(arg["mu"])
    # kappa = float(arg["kappa"])
    # H = float(arg["Hardness"])

    fileName = "sim" + str(run_id)

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        rho=rho,
        youngs_modulus=E,
        poisson_ratio=nu,
    )

    material.JohnsonCookHardening(A=A, B=B, n=n)
    # material.JohnsonCookDamage(d1=D1, d2=D2, d3=D3)
    # material.DamageEvolution(kc=kc, uts=uts, E=E, nu=nu)
    # material.DamageEvolution(u_pl_f)
    material.SectionAssignment()
    # material.UpdateFrictionAndWear(mu, include_wear, kappa)
    material.UpdateFrictionAndWear(mu)

    run_job_and_wait(fileName, "Model-1")

    PostProcess(jobName, fileName, arg)

mdb.close()
cleanupAbaqusJunk()
