from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
from abaqus import *
from abaqusConstants import *
import os

import sys


def test_progressive_load_scratch():
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    jobName = "test_progressive_load_scratch"
    rundir = os.path.join("test_runs", jobName)
    if not os.path.exists(rundir):
        os.makedirs(rundir)

    # Change current working directory
    os.chdir(rundir)
    num_cpus = 6
    num_domains = num_cpus

    E_modulus = 200000.0
    density = 7.8e-9
    poisson = 0.3

    SubstrateMeshSizeX = 0.04
    SubstrateMeshSizeY = 0.04
    SubstrateMeshSizeZ = 0.04
    depth = -50e3

    ScratchModel, SubstratePart, SubstrateSet = ScratchModelSetup(
        depth=depth,
        IndenterToUse="RockwellIndenter",
        SubstrateMeshSizeX=SubstrateMeshSizeX,
        SubstrateMeshSizeY=SubstrateMeshSizeY,
        SubstrateMeshSizeZ=SubstrateMeshSizeZ,
    )

    material = SubstrateMaterialAssignment(
        ScratchModel,
        SubstratePart,
        SubstrateSet,
        rho=density,
        youngs_modulus=E_modulus,
        poisson_ratio=poisson,
    )

    material.IsotrpopicHardening(yield_strength=600.0, n=0.2)
    material.SectionAssignment()

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

    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()
    mdb.close()


if __name__ == "__main__":
    test_progressive_load_scratch()
