from abaqus import *
from abaqusConstants import *
from ProgressiveLoadScratch.PostProcessing import *
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os

# from material_parameters import parameters
import numpy as np
import shutil


### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "MeshConvergence"

# Array of shape (#mesh combinations, 3). Rows are mesh combinations, coloums are [meshSizeX, meshSizeY, meshSizeZ],
meshSizes = np.array(
    [
        [0.040, 0.040, 0.040],
        [0.030, 0.030, 0.030],
        [0.020, 0.020, 0.020],
        [0.010, 0.020, 0.020],
        [0.010, 0.020, 0.010],
        [0.010, 0.010, 0.010],
        # [0.005, 0.010, 0.010],
        # [0.010, 0.010, 0.005],
        # [0.010, 0.005, 0.010],
        # [0.005, 0.010, 0.005],
        # [0.010, 0.005, 0.005],
        # [0.005, 0.005, 0.010],
        # [0.005, 0.005, 0.005], # Takes too long
        # [0.002, 0.002, 0.002], # Takes too long
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


depth = -25e-3
rho = 7.8e-09
E = 210000.0
nu = 0.3
A = 1430.0
B = 2545.0
n = 0.7
D1 = 0.028
D2 = 1.0
D3 = -0.916
uts = 1600.0
kc = 2000
mu = 0.2


for id, meshSize in enumerate(meshSizes):

    fileName = (
        "mesh" + str(meshSize[0]) + "_" + str(meshSize[1]) + "_" + str(meshSize[2])
    )
    ScratchModel, SubstratePart, SubstrateSet = ScratchModelSetup(
        depth=depth,
        SubstrateSizeY=meshSize[0],
        SubstrateSizeX=meshSize[1],
        SubstrateSizeZ=meshSize[2],
        target_time_increment=0.0,
    )

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        SubstrateSet,
        rho=rho,
        youngs_modulus=E,
        poisson_ratio=nu,
    )

    material.JohnsonCookHardening(A=A, B=B, n=n)
    material.JohnsonCookDamage(d1=D1, d2=D2, d3=D3)
    material.DamageEvolution(kc=kc, uts=uts, E=E, nu=nu)
    material.SectionAssignment()
    material.UpdateFriction(mu)

    #### ------------------------------ ####
    #           Create Job
    #### ------------------------------ ####
    mdb.Job(
        activateLoadBalancing=False,
        atTime=None,
        contactPrint=OFF,
        description="",
        echoPrint=OFF,
        explicitPrecision=SINGLE,
        historyPrint=OFF,
        memory=90,
        memoryUnits=PERCENTAGE,
        model="Model-1",
        modelPrint=OFF,
        multiprocessingMode=MPI,
        name=jobName,
        nodalOutputPrecision=SINGLE,
        numCpus=num_cpus,
        numDomains=num_domains,
        parallelizationMethodExplicit=DOMAIN,
        queue=None,
        resultsFormat=ODB,
        scratch="",
        type=ANALYSIS,
        userSubroutine="",
        waitHours=0,
        waitMinutes=0,
    )

    #### ------------------------------ ####
    #             Submit Job
    #### ------------------------------ ####
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()

    PostProcess(jobName, fileName)

    mdb.close()

    sta_file = jobName + ".sta"
    target_dir = "SimDataOutputs/"
    new_name = fileName + ".sta"
    dst_file = os.path.join(target_dir, new_name)
    shutil.move(sta_file, dst_file)
