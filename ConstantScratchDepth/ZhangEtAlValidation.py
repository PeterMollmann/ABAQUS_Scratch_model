from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import *
import os
import numpy as np


def ScratchModelParameterSweep(
    jobName="ExplicitRigidIndenterScratch",
    sigma_y=500.0,
    strain_hardening_index=0.2,
    depth=-50e-3,
    friction_coefficient=0.0,
    E_modulus=200000.0,
    density=7.8e-9,
    poisson=0.3,
):
    ### Remember to use consistent units! This document uses SI(mm) units ###

    # Make selecting geometries easier to understand. Makes for more generalisable code
    session.journalOptions.setValues(
        replayGeometry=COORDINATE, recoverGeometry=COORDINATE
    )

    ScratchModel = mdb.models["Model-1"]

    # Naming - Not important
    SubstrateName = "TestSpecimen"
    IndenterName = "Indenter"
    MaterialName = "MyMaterial"

    # Parameters for sketching
    sheet_size = 10

    # Geometry coordinates for specimen
    xs1 = 0.0  # x coordinate of first point
    ys1 = 0.0  # y coordinate of first point
    xs2 = 0.96  # x coordinate of second point
    ys2 = 0.64  # y coordinate of second point
    zs1 = 0.0  # z coordinate of first point
    zs2 = 2.88  # z coordinate of second point - extrude depth

    # Geometry coordinates for indenter
    xi1 = -0.2
    yi1 = xi1 + ys2
    xi2 = -xi1
    yi2 = xi2 + ys2
    indenter_angle = -60  # should be negative to get correct extrude
    z_factor = (
        1
        / cos((90 + indenter_angle) * pi / 180)
        * sin((90 + indenter_angle) * pi / 180)
    )
    zi = (
        xi2 * z_factor
    )  # is the sum of 90 and the angle as the angel has to be negative

    # Datum plane offsets
    dpo_x = ys2 / 2.0
    dpo_y = 0.10

    # Meshing parameters
    IndenterMinSize = 0.001
    IndenterMaxSize = 0.025

    SubstrateSizeY = 0.025
    SubstrateSizeX = 0.02
    SubstrateSizeZ = 0.04

    CoarseMeshSize0 = 0.05
    CoarseMeshSize1 = 0.15
    CoarseMeshSize2 = 0.3

    # Material characteristics
    rho = density  # density [tonne/mm^3]
    youngs_modulus = E_modulus  # [MPa]
    poisson_ratio = poisson
    yield_strength = sigma_y  # [MPa]
    mu_s = friction_coefficient
    n = strain_hardening_index  # strain hardening index

    # Analysis time
    indentation_time = 0.0001  # [s]
    scratch_time = 0.001  # [s]

    # Scratch parameters
    scratch_depth = depth  # [mm]
    scratch_length = 2  # [mm]

    sample_frequency_indentation = indentation_time / 2.0
    sample_frequency_scratching = scratch_time / 2.0
    sample_force_frequency_scratching = scratch_time / 40.0
    mass_scale = 1e4

    # Parallelisation
    num_cpus = 6
    num_domains = num_cpus

    #### ------------------------------ ####
    #         Substrate geometry
    #### ------------------------------ ####
    ScratchModel.ConstrainedSketch(
        name="__profile__", sheetSize=sheet_size
    )  # Make sketching sheet
    Sketch = ScratchModel.sketches["__profile__"]
    Sketch.rectangle(
        point1=(xs1, ys1), point2=(xs2, ys2)
    )  # Create rectangle based on coordinates defined for the substrate
    ScratchModel.Part(
        dimensionality=THREE_D, name=SubstrateName, type=DEFORMABLE_BODY
    )  # Substrate is modelled as a deformable body
    ScratchModel.parts[SubstrateName].BaseSolidExtrude(
        depth=zs2, sketch=Sketch
    )  # Extrude sketch by "zs2" amount
    del Sketch
    SubstratePart = ScratchModel.parts[SubstrateName]

    # Create datum plane for partition
    SubstratePart.DatumPlaneByPrincipalPlane(
        offset=dpo_x, principalPlane=YZPLANE
    )  # Datum plane for z axis partition
    SubstratePart.DatumPlaneByPrincipalPlane(
        offset=(ys2 - dpo_y), principalPlane=XZPLANE
    )  # Datum plane for y axis partition

    # Make partition of substrate for mesh refinement along scratch path
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((xs2, ys2, zs2),),
        ),
        datumPlane=SubstratePart.datums[2],
    )  # partition along z axis
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs2, ys2, zs2),),
        ),
        datumPlane=SubstratePart.datums[3],
    )  # partition along y axis

    #### ------------------------------ ####
    #          Indenter geometry
    #### ------------------------------ ####
    ScratchModel.ConstrainedSketch(name="__profile__", sheetSize=sheet_size)
    Sketch = ScratchModel.sketches["__profile__"]
    Sketch.rectangle(
        point1=(xi1, yi1),
        point2=(xi2, yi2),
    )
    ScratchModel.Part(
        dimensionality=THREE_D, name=IndenterName, type=DISCRETE_RIGID_SURFACE
    )
    ScratchModel.parts[IndenterName].BaseSolidExtrude(
        depth=1.0, draftAngle=indenter_angle, sketch=Sketch
    )
    del Sketch
    IndenterPart = ScratchModel.parts[IndenterName]
    IndenterPart.RemoveCells(
        cellList=(
            IndenterPart.cells.findAt(
                (xi1, yi1, 0.0),
            ),
        )
    )

    #### ------------------------------ ####
    #         Materials Indenter
    #### ------------------------------ ####
    # Creating reference point and inertia for indenter
    indenter_set = IndenterName + "Set"
    IndenterPart.ReferencePoint(
        point=IndenterPart.vertices.findAt(
            ((xi1 + xi2) / 2.0, (yi1 + yi2) / 2.0, zi),
        )
    )
    IndenterPart.Set(
        name=indenter_set, referencePoints=(IndenterPart.referencePoints[3],)
    )
    IndenterPart.engineeringFeatures.PointMassInertia(
        alpha=0.0,
        composite=0.0,
        i11=0.1,
        i22=0.1,
        i33=0.1,
        mass=1.0,
        name="IndenterInertia",
        region=IndenterPart.sets[indenter_set],
    )

    #### ------------------------------ ####
    #         Materials Substrate
    #### ------------------------------ ####
    # Creating the material for the substrate
    ScratchModel.Material(name=MaterialName)
    ScratchModel.materials[MaterialName].Density(table=((rho,),))
    ScratchModel.materials[MaterialName].Elastic(
        table=((youngs_modulus, poisson_ratio),)
    )

    total_strain_data = np.append(np.arange(0.0001, 0.1, 0.01), np.linspace(0.1, 2, 20))
    plastic_behaviour = []
    yield_strain = yield_strength / youngs_modulus
    K = youngs_modulus * yield_strain ** (1 - n)
    plastic_behaviour.append((yield_strength, 0))
    for total_strain in total_strain_data:
        if total_strain > yield_strain:
            yield_strength_ = round(K * total_strain**n, 5)
            plastic_strain = round(total_strain - yield_strength_ / youngs_modulus, 5)
            plastic_behaviour.append(
                (yield_strength_, plastic_strain)
            )  # (Yield strength, plastic strain) [MPa, -]

    plastic_behaviour = tuple(plastic_behaviour)
    ScratchModel.materials[MaterialName].Plastic(table=plastic_behaviour)

    # Creating section
    SubstrateSet = SubstrateName + "Set"
    ScratchModel.HomogeneousSolidSection(
        material=MaterialName, name="SubstrateSection", thickness=None
    )
    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs1, ys2, zs1),),
            ((xs2, ys2, zs1),),
            ((xs2, ys1, zs1),),
        ),
        name=SubstrateSet,
    )
    SubstratePart.SectionAssignment(
        offset=0.0,
        offsetField="",
        offsetType=MIDDLE_SURFACE,
        region=SubstratePart.sets[SubstrateSet],
        sectionName="SubstrateSection",
        thicknessAssignment=FROM_SECTION,
    )

    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys2, zs1),),
        ),
        name="RefinedArea",
    )

    #### ------------------------------ ####
    #               Meshing
    #### ------------------------------ ####
    # Meshing of indenter
    eps = 0.075  # Small offset to help select edges
    IndenterPart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end2Edges=IndenterPart.edges.findAt(
            ((xi1 + eps, yi1 + eps, eps * z_factor),),
            ((xi1 + eps, yi2 - eps, eps * z_factor),),
            ((xi2 - eps, yi2 - eps, eps * z_factor),),
            ((xi2 - eps, yi1 + eps, eps * z_factor),),
        ),
        maxSize=IndenterMaxSize,
        minSize=IndenterMinSize,
    )
    IndenterPart.generateMesh()

    # Meshing of substrate
    SubstratePart.setMeshControls(
        elemShape=TET,
        regions=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs1, ys2, zs1),),
            ((xs2, ys1, zs1),),
            ((xs2, ys2, zs1),),
        ),
        technique=FREE,
    )
    SubstratePart.setElementType(
        elemTypes=(
            ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT),
            ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
            ElemType(
                elemCode=C3D4,
                elemLibrary=EXPLICIT,
                secondOrderAccuracy=OFF,
                distortionControl=DEFAULT,
            ),
        ),
        regions=(
            SubstratePart.cells.findAt(
                ((xs1, ys1, zs1),),
                ((xs1, ys2, zs1),),
                ((xs2, ys1, zs1),),
                ((xs2, ys2, zs1),),
            ),
        ),
    )

    ### Refined area meshing ###
    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys2 - dpo_y / 2.0, zs1),),
            ((dpo_x, ys2 - dpo_y / 2.0, zs1),),
            ((xs1, ys2 - dpo_y / 2.0, zs2),),
            ((dpo_x, ys2 - dpo_y / 2.0, zs2),),
        ),
        size=SubstrateSizeY,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((dpo_x / 2.0, ys2 - dpo_y, zs1),),
            ((dpo_x / 2.0, ys2 - dpo_y, zs2),),
            ((dpo_x / 2.0, ys2, zs2),),
            ((dpo_x / 2.0, ys2, zs1),),
        ),
        size=SubstrateSizeX,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys2 - dpo_y, (zs2 + zs1) / 2.0),),
            ((xs1, ys2, (zs2 + zs1) / 2.0),),
            ((dpo_x, ys2, (zs2 + zs1) / 2.0),),
        ),
        size=SubstrateSizeZ,
    )

    ### Meshing the rest ###
    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys1, (zs2 + zs1) / 2.0),),
            ((dpo_x, ys1, (zs2 + zs1) / 2.0),),
            ((xs2, ys2, (zs2 + zs1) / 2.0),),
            ((xs2, ys2 - dpo_y, (zs2 + zs1) / 2.0),),
            ((dpo_x / 2.0, ys1, zs1),),
            ((dpo_x / 2.0, ys1, zs2),),
        ),
        size=CoarseMeshSize1,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs2, ys2 - dpo_y / 2.0, zs1),),
            ((xs2, ys2 - dpo_y / 2.0, zs2),),
        ),
        size=dpo_y / 2.0,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((dpo_x, (ys2 - dpo_y) / 2.0, zs2),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((dpo_x, (ys2 - dpo_y) / 2.0, zs1),),
        ),
        maxSize=CoarseMeshSize1,
        minSize=CoarseMeshSize0,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((dpo_x, ys2 - dpo_y, (zs2 + zs1) / 2.0),),
        ),
        size=CoarseMeshSize0,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FIXED,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs2, ys1, (zs2 + zs1) / 2.0),),
        ),
        size=CoarseMeshSize2,
    )

    SubstratePart.generateMesh()

    #### ------------------------------ ####
    #            Assembly
    #### ------------------------------ ####
    indenterInstanceName = IndenterName + "Inst"
    specimenInstanceName = SubstrateName + "Inst"
    ScratchModelAssembly = ScratchModel.rootAssembly
    ScratchModelAssembly.DatumCsysByDefault(CARTESIAN)
    ScratchModelAssembly.Instance(
        dependent=ON, name=indenterInstanceName, part=IndenterPart
    )
    IndenterInstance = ScratchModelAssembly.instances[indenterInstanceName]
    ScratchModelAssembly.Instance(
        dependent=ON, name=specimenInstanceName, part=SubstratePart
    )
    SpecimenInstance = ScratchModelAssembly.instances[specimenInstanceName]
    ScratchModelAssembly.translate(
        instanceList=(indenterInstanceName,), vector=(0.0, 0.0, -zi)
    )
    ScratchModelAssembly.rotate(
        angle=90.0,
        axisDirection=(10.0, 0.0, 0.0),
        axisPoint=(xs1, ys2, zs1),
        instanceList=(indenterInstanceName,),
    )
    ScratchModelAssembly.translate(
        instanceList=(indenterInstanceName,),
        vector=(0.0, 0.0, (zs2 - scratch_length) / 2.0),
    )

    #### ------------------------------ ####
    #           Indentation Step
    #### ------------------------------ ####
    # Create explicit dynamics step with mass scaling
    ScratchModel.ExplicitDynamicsStep(
        improvedDtMethod=ON,
        massScaling=(
            (
                SEMI_AUTOMATIC,
                MODEL,
                AT_BEGINNING,
                mass_scale,
                0.0,
                None,
                0,
                0,
                0.0,
                0.0,
                0,
                None,
            ),
        ),
        name="IndentationStep",
        previous="Initial",
        timePeriod=indentation_time,
    )

    # Fixed constraint on the bottom two faces
    fixedBCSet = "fixedBCSet"
    ScratchModelAssembly.Set(
        faces=SpecimenInstance.faces.findAt(
            ((dpo_x / 2.0, ys1, (zs2 + zs1) / 2.0),),
            (((dpo_x + xs2) / 2.0, ys1, (zs2 + zs1) / 2.0),),
        ),
        name=fixedBCSet,
    )
    ScratchModel.EncastreBC(
        createStepName="IndentationStep",
        localCsys=None,
        name="Fixed_constraint",
        region=ScratchModelAssembly.sets[fixedBCSet],
    )

    # Symmetry constraint
    symmetryBCSet = "symmetryBCSet"
    ScratchModelAssembly.Set(
        faces=SpecimenInstance.faces.findAt(
            ((xs1, (ys2 + ys1) / 2.0, (zs2 + zs1) / 2.0),),
            ((xs1, ys2 - dpo_y / 2.0, (zs2 + zs1) / 2.0),),
        ),
        name=symmetryBCSet,
    )
    ScratchModel.XsymmBC(
        createStepName="IndentationStep",
        localCsys=None,
        name="x_axis_symmetry",
        region=ScratchModelAssembly.sets[symmetryBCSet],
    )

    # Defining the indenter movement bc
    ScratchModelAssembly.Set(
        name=indenter_set, referencePoints=(IndenterInstance.referencePoints[3],)
    )
    ScratchModel.TabularAmplitude(
        data=(
            (0.0, 0.0),
            (indentation_time, 1.0),
            (indentation_time + scratch_time, 1.0),
            (2 * indentation_time + scratch_time, 0.0),
        ),
        name="Amp-1",
        smooth=SOLVER_DEFAULT,
        timeSpan=TOTAL,
    )

    ScratchModel.DisplacementBC(
        amplitude="Amp-1",
        createStepName="IndentationStep",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="IndenterIndentation",
        region=ScratchModelAssembly.sets[indenter_set],
        u1=UNSET,
        u2=scratch_depth,
        u3=UNSET,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
    )

    ScratchModel.DisplacementBC(
        amplitude=UNSET,
        createStepName="IndentationStep",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="IndenterConstraint",
        region=ScratchModelAssembly.sets[indenter_set],
        u1=SET,
        u2=UNSET,
        u3=SET,
        ur1=SET,
        ur2=SET,
        ur3=SET,
    )

    ### History and field output requests ###
    del ScratchModel.fieldOutputRequests["F-Output-1"]
    del ScratchModel.historyOutputRequests["H-Output-1"]

    ScratchModel.HistoryOutputRequest(
        createStepName="IndentationStep",
        name="ReactionForces",
        rebar=EXCLUDE,
        region=ScratchModelAssembly.sets[indenter_set],
        sectionPoints=DEFAULT,
        timeInterval=sample_frequency_indentation,
        variables=("RF1", "RF2", "RF3"),
    )

    ScratchModel.FieldOutputRequest(
        createStepName="IndentationStep",
        name="FieldOutput",
        timeInterval=sample_frequency_indentation,
        variables=PRESELECT,
    )
    ScratchModel.FieldOutputRequest(
        createStepName="IndentationStep",
        name="ContactForce",
        timeInterval=sample_frequency_indentation,
        variables=("CFORCE",),
    )
    #### ------------------------------ ####
    #           Scratching Step
    #### ------------------------------ ####
    ScratchModel.ExplicitDynamicsStep(
        improvedDtMethod=ON,
        name="ScratchingStep",
        previous="IndentationStep",
        timePeriod=scratch_time,
    )

    ScratchModel.TabularAmplitude(
        data=((indentation_time, 0.0), (indentation_time + scratch_time, 1.0)),
        name="Amp-2",
        smooth=SOLVER_DEFAULT,
        timeSpan=TOTAL,
    )

    ScratchModel.boundaryConditions["IndenterConstraint"].setValuesInStep(
        stepName="ScratchingStep", u3=FREED
    )

    ScratchModel.DisplacementBC(
        amplitude="Amp-2",
        createStepName="ScratchingStep",
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="IndenterScratching",
        region=ScratchModelAssembly.sets[indenter_set],
        u1=UNSET,
        u2=UNSET,
        u3=scratch_length,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
    )

    ScratchModel.fieldOutputRequests["FieldOutput"].setValuesInStep(
        stepName="ScratchingStep", timeInterval=sample_frequency_scratching
    )
    ScratchModel.fieldOutputRequests["ContactForce"].setValuesInStep(
        stepName="ScratchingStep", timeInterval=sample_frequency_scratching
    )

    ScratchModel.historyOutputRequests["ReactionForces"].setValuesInStep(
        stepName="ScratchingStep", timeInterval=sample_force_frequency_scratching
    )

    #### ------------------------------ ####
    #           Unloading Step
    #### ------------------------------ ####
    ScratchModel.ExplicitDynamicsStep(
        improvedDtMethod=ON,
        name="UnloadingStep",
        previous="ScratchingStep",
        timePeriod=indentation_time,
    )
    ScratchModel.TabularAmplitude(
        data=(
            (indentation_time + scratch_time, 0.0),
            (2 * indentation_time + scratch_time, 1.0),
        ),
        name="Amp-3",
        smooth=SOLVER_DEFAULT,
        timeSpan=TOTAL,
    )

    ScratchModel.fieldOutputRequests["FieldOutput"].setValuesInStep(
        stepName="UnloadingStep", timeInterval=sample_frequency_indentation
    )
    ScratchModel.fieldOutputRequests["ContactForce"].setValuesInStep(
        stepName="UnloadingStep", timeInterval=sample_frequency_indentation
    )

    ScratchModel.historyOutputRequests["ReactionForces"].setValuesInStep(
        stepName="UnloadingStep", timeInterval=sample_frequency_indentation
    )

    #### ------------------------------ ####
    #          Contact modelling
    #### ------------------------------ ####
    ScratchModel.ContactProperty("IntProp-1")
    ScratchModel.interactionProperties["IntProp-1"].TangentialBehavior(
        formulation=FRICTIONLESS
    )
    ScratchModel.interactionProperties["IntProp-1"].NormalBehavior(
        allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT,
        pressureOverclosure=HARD,
    )

    # Use abaqus general contact
    ScratchModel.ContactExp(createStepName="Initial", name="Int-1")
    ScratchModel.interactions["Int-1"].includedPairs.setValuesInStep(
        stepName="Initial", useAllstar=ON
    )
    ScratchModel.interactions["Int-1"].contactPropertyAssignments.appendInStep(
        assignments=((GLOBAL, SELF, "IntProp-1"),), stepName="Initial"
    )

    #### ------------------------------ ####
    #         Final touches
    #### ------------------------------ ####
    # Rotate indenter 45 degress around x axis
    ScratchModelAssembly.rotate(
        angle=45.0,
        axisDirection=(xs1, ys2, zs1),
        axisPoint=(xs1, ys2, (zs2 - scratch_length) / 2.0),
        instanceList=(indenterInstanceName,),
    )

    # create node set for post processing ease
    contactSurfaceSet = "contactSurfaceSet"
    ScratchModelAssembly.Set(
        faces=SpecimenInstance.faces.findAt(((dpo_x / 2.0, ys2, (zs2 + zs1) / 2.0),)),
        name=contactSurfaceSet,
    )

    #### ------------------------------ ####
    #           Create Job
    #### ------------------------------ ####
    # jobName = "SY"+str(sigma_y)+"_n"+str(strain_hardening_index)+"_d"+str(depth)+"_E"+str(E_modulus)+"_Friction"+str(friction_coefficient) + "_Poisson"+str(poisson)
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
    mdb.close()
