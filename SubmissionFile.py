from abaqus import *
from abaqusConstants import *
import ProgressiveLoadScratch.Constants as C
from ProgressiveLoadScratch.PostProcessing import PostProcess
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from cleanup import cleanupAbaqusJunk
from ProgressiveLoadScratch.helpers import run_job_and_wait

from material_parameters.halton_discrete_material_parameter_sweep import parameters

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

start_from_sim_id = 37  # Set to desired starting ID to skip completed simulations
stop_at_id = 100  # Set to desired stopping ID. Runs this ID simulation

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
    mu = float(arg["mu"])

    fileName = "sim" + str(run_id)

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        rho=rho,
        youngs_modulus=E,
        poisson_ratio=nu,
    )

    material.JohnsonCookHardening(A=A, B=B, n=n)
    material.SectionAssignment()
    material.UpdateFrictionAndWear(mu)

    run_job_and_wait(jobName)

    PostProcess(jobName, fileName, arg)

    if int(run_id) == stop_at_id:
        break

mdb.close()
cleanupAbaqusJunk()
