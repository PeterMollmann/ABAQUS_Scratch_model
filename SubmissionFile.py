from abaqus import *
from abaqusConstants import *
import ProgressiveLoadScratch.Constants as C
from ProgressiveLoadScratch.PostProcessing import PostProcess
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from cleanup import cleanupAbaqusJunk
from ProgressiveLoadScratch.helpers import run_job_and_wait
import shutil

from material_parameters.halton_discrete_material_parameter_sweep import parameters

### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "MaterialSweepNew"
# jobName = "OutlierInvestigation"


meshSize = [0.030, 0.020, 0.010, 0.008, 0.006, 0.004, 0.002]
meshSizeIdx = 4


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
    mass_scale=5e5,
    use_ALE=True,
    include_wear=include_wear,
)

# 34
start_from_sim_id = 7  # Set to desired starting ID to skip completed simulations
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
    # fileName = (
    #     "Hard_contact_default_ALE_sweep1per20_sweep1per200" + "_sim" + str(run_id)
    # )
    # fileName = "scale_factor_cssf2.0_issf1.0_ocf2.0" + "_sim" + str(run_id)
    # fileName = "test" + "_sim" + str(run_id)

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
