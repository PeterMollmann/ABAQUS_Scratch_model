from abaqus import *
from abaqusConstants import *
from PostProcessing import *
from ProgressiveLoading.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoading.SubstrateMaterial import SubstrateMaterialAssignment
from itertools import product


depth = -50e-3
friction_coefficient = 0.0

# Material elastic prpperties
E_modulus = 200000.0
density = 7.8e-9

# Material density
poisson = 0.3

# Material hardening properties - Isotropic hardening
# yield_stress_sweep = [200.0]
# strain_hardening_sweep = [0.2]
# materialIterationProduct = product(yield_stress_sweep, strain_hardening_sweep)

# Material hardening properties - Johnson-Cook
A = [200.0]
B = [100.0]
n = [0.2]
m = 0.0
Tm = 0.0
Tt = 0.0

# Material damage model - Johnson-Cook damage initiation
d1 = 0.028
d2 = 1
d3 = -0.916
d4 = 0
d5 = 0
Sr = 0


materialIterationProduct = product(A, B, n)

# Parallelisation
num_cpus = 6
num_domains = num_cpus

ScratchModel, SubstratePart, SubstrateSet = ScratchModelSetup(
    depth=depth,
    IndenterToUse="RockwellIndenter",
)

for arg in materialIterationProduct:

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        SubstrateSet,
        rho=density,
        youngs_modulus=E_modulus,
        poisson_ratio=poisson,
    )
    # material.IsotrpopicHardening(yield_strength=arg[0], n=arg[1])
    material.JohnsonCookHardening(A=arg[0], B=arg[1], n=arg[2], m=m, Tm=Tm, Tt=Tt)
    # material.JohnsonCookDamage(d1=d1, d2=d2, d3=d3, d4=d4, d5=d5, Tm=Tm, Tt=Tt, Sr=Sr)
    material.SectionAssignment()

    jobName = "ProgressiveLoadScratchTest"

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

    # PostProcess(
    #     jobName,
    #     sigma_y=sigma_y,
    #     strain_hardening_index=n,
    #     depth=depth,
    #     friction_coefficient=friction_coefficient,
    #     E_modulus=E_modulus,
    #     density=density,
    #     poisson=poisson,
    # )

mdb.close()
