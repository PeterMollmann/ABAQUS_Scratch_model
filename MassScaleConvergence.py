from abaqus import *
from abaqusConstants import *
from ProgressiveLoadScratch.PostProcessing import *
from ProgressiveLoadScratch.ProgressiveLoadScratchTest import ScratchModelSetup
from ProgressiveLoadScratch.SubstrateMaterial import SubstrateMaterialAssignment
import os

# from material_parameters import parameters
import shutil
from cleanup import cleanupAbaqusJunk
import ProgressiveLoadScratch.Constants as C

### ---------------- ###
# SETTINGS
### ---------------- ###
jobName = "MassScaleConvergence"

# Array of shape (#mesh combinations, 3). Rows are mesh combinations, coloums are [meshSizeX, meshSizeY, meshSizeZ],
# massScales = [1e5, 5e4, 1e4, 1e3, 1e2, 1e1]
# massScales = [1e3, 1e2, 1e1]
massScales = [4e5]


meshSize = [0.040, 0.030, 0.020, 0.010, 0.005]
meshSizeIdx = 3


# Change abaqus working directory
rundir = os.path.join("runs", "MassScaleConvergence")
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

include_wear = True

for id, massScale in enumerate(massScales):

    fileName = "MassScale" + str(massScale)
    ScratchModel, SubstratePart = ScratchModelSetup(
        SubstrateSizeY=meshSize[meshSizeIdx],
        SubstrateSizeX=meshSize[meshSizeIdx],
        SubstrateSizeZ=meshSize[meshSizeIdx],
        mass_scale=massScale,
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
    # material.DamageEvolution(u_pl_f)
    material.SectionAssignment()
    material.UpdateFrictionAndWear(mu, include_wear, kappa)

    #### ------------------------------ ####
    #           Create Job
    #### ------------------------------ ####
    job = mdb.Job(
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
        numCpus=C.num_cpus,
        numDomains=C.num_domains,
        parallelizationMethodExplicit=DOMAIN,
        queue=None,
        resultsFormat=ODB,
        scratch="",
        type=ANALYSIS,
        userSubroutine="",
        waitHours=0,
        waitMinutes=0,
    )

    job.writeInput(consistencyChecking=OFF)

    if include_wear:
        # Go into .inp file and locate:
        #
        # *Surface Property, name=IntProp-2
        # *WEAR SURFACE PROPERTIES, FRIC COEF DEPENDENT=YES, UNITLESS WEAR COEF=NO
        # 0.0001, , , , **
        #
        # There is a bug here, it should be:
        #
        # 0.0001, , , ,
        # **
        #
        # This is fixed here.
        trigger_found = False
        inpFile = jobName + ".inp"
        with open(inpFile, "r") as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if "*WEAR SURFACE PROPERTIES" in line:
                trigger_found = True
                # next line is the one we want to change
                if i + 1 < len(lines):
                    old_line = lines[i + 1]

                    new_line = f" {kappa}, , , , \n**\n"

                    lines[i + 1] = new_line

                break

        if not trigger_found:
            raise ValueError(
                f"Could not find the trigger line in the inp file. Old line: {old_line}"
            )

        # Write the modified content back
        # inpFileBackup = jobName + "Backup" + ".inp"
        with open(inpFile, "w") as f:
            f.writelines(lines)

    job = mdb.JobFromInputFile(
        name=jobName,
        inputFileName=jobName + ".inp",
        activateLoadBalancing=False,
        atTime=None,
        # contactPrint=OFF,
        # description="",
        # echoPrint=OFF,
        explicitPrecision=SINGLE,
        # historyPrint=OFF,
        memory=90,
        memoryUnits=PERCENTAGE,
        # model="Model-1",
        # modelPrint=OFF,
        multiprocessingMode=MPI,
        nodalOutputPrecision=SINGLE,
        numCpus=C.num_cpus,
        numDomains=C.num_domains,
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
    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()

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
