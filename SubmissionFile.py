from abaqus import *
from abaqusConstants import *
from PostProcessing import *
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os
from material_parameters import parameters

# import pandas as pd

### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "ProgressiveLoadScratchTest"
indenter = "RockwellIndenter"
materialModel = "JohnsonCook"
damageModel = False

depth = -100e-3
friction_coefficient = 0.0

# density = 7.8e-9

# Parallelisation
num_cpus = 6
num_domains = num_cpus


# Change abaqus working directory
rundir = os.path.join("runs", jobName)
if not os.path.exists(rundir):
    os.makedirs(rundir)
os.chdir(rundir)


# Setup scratch model. Only needs to be called once
ScratchModel, SubstratePart, SubstrateSet = ScratchModelSetup(
    depth=depth, IndenterToUse=indenter
)


for arg in parameters:
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

    fileName = "sim" + str(run_id)

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
    material.SectionAssignment()

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
