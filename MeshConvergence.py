from abaqus import *
from abaqusConstants import *
from ProgressiveLoadScratch.PostProcessing import *
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from material_parameters import parameters
import numpy as np

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
        [0.010, 0.010, 0.010],
        [0.005, 0.005, 0.005],
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


depth = -50e-3
arg = parameters[0]
run_id = arg["id"]
rho = float(arg["rho"])
E = float(arg["E"])
nu = float(arg["nu"])
A = float(arg["A"])
B = float(arg["B"])
n = float(arg["n"])
D1 = float(arg["D1"])
D2 = float(arg["D2"])
D3 = float(arg["D3"])
uts = float(arg["uts"])
kc = float(arg["kc"])
mu = float(arg["mu"])


for id, meshSize in enumerate(meshSizes):

    fileName = "sim" + str(id)
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
    material.DamageEvolution(kc=kc, E=E, nu=nu)
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
