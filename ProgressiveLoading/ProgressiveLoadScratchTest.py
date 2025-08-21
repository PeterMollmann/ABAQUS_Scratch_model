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
from PyramidIndenter import PyramidIndenter
from SubstrateGeneration import (
    SubstrateGeneration,
    SubstrateMeshing,
)


def ScratchModelSetup(
    depth=-50e-3,
):
    ### Remember to use consistent units! This document uses SI(mm) units ###

    # Make selecting geometries easier to understand. Makes for more generalisable code
    session.journalOptions.setValues(
        replayGeometry=COORDINATE, recoverGeometry=COORDINATE
    )

    ScratchModel = mdb.models["Model-1"]

    # Parameters for sketching
    sheet_size = 10

    # Geometry coordinates for specimen
    xs1 = 0.0  # x coordinate of first point
    ys1 = 0.0  # y coordinate of first point
    xs2 = 0.96  # x coordinate of second point
    ys2 = 0.64  # y coordinate of second point
    zs1 = 0.0  # z coordinate of first point
    zs2 = 6.0  # z coordinate of second point - extrude depth

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

    # Analysis time
    indentation_time = 0.0001  # [s]
    scratch_time = 0.001  # [s]

    # Scratch parameters
    scratch_depth = depth  # [mm]
    scratch_length = 5  # [mm]

    sample_frequency_indentation = indentation_time / 1.0
    sample_frequency_scratching = scratch_time / 1.0
    sample_force_frequency_scratching = scratch_time / 100.0
    mass_scale = 1e4

    # Create and mesh substate
    ScratchModel, SubstratePart, _, SubstrateSet = SubstrateGeneration(
        ScratchModel,
        xs1,
        ys1,
        xs2,
        ys2,
        zs1,
        zs2,
        dpo_x,
        dpo_y,
        sheet_size,
    )

    SubstratePart = SubstrateMeshing(
        SubstratePart,
        xs1,
        ys1,
        zs1,
        xs2,
        ys2,
        zs2,
        dpo_x,
        dpo_y,
        CoarseMeshSize0,
        CoarseMeshSize1,
        CoarseMeshSize2,
        SubstrateSizeX,
        SubstrateSizeY,
        SubstrateSizeZ,
    )

    # Create and mesh indenter
    ScratchModel, IndenterPart, indenter_set, IndenterName, indenterHeight = (
        PyramidIndenter(
            ScratchModel,
            xs2,
            ys2,
            IndenterMaxSize,
            IndenterMinSize,
            sheet_size=sheet_size,
        )
    )

    #### ------------------------------ ####
    #            Assembly
    #### ------------------------------ ####
    indenterInstanceName = IndenterName + "Inst"
    specimenInstanceName = "SubstrateInst"
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
        instanceList=(indenterInstanceName,), vector=(0.0, 0.0, -indenterHeight)
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
    #           Scratching Step
    #### ------------------------------ ####
    # Create explicit dynamics step with mass scaling
    StepName1 = "ProgressiveScratchStep"
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
        name=StepName1,
        previous="Initial",
        timePeriod=scratch_time,
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
        createStepName=StepName1,
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
        createStepName=StepName1,
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
            # (indentation_time, 1.0),
            (scratch_time, 1.0),
            (indentation_time + scratch_time, 0.0),
        ),
        name="Amp-1",
        smooth=SOLVER_DEFAULT,
        timeSpan=TOTAL,
    )

    ScratchModel.DisplacementBC(
        amplitude="Amp-1",
        createStepName=StepName1,
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="IndenterScratching",
        region=ScratchModelAssembly.sets[indenter_set],
        u1=UNSET,
        u2=scratch_depth,
        u3=scratch_length,
        ur1=UNSET,
        ur2=UNSET,
        ur3=UNSET,
    )

    ScratchModel.DisplacementBC(
        amplitude=UNSET,
        createStepName=StepName1,
        distributionType=UNIFORM,
        fieldName="",
        fixed=OFF,
        localCsys=None,
        name="IndenterConstraint",
        region=ScratchModelAssembly.sets[indenter_set],
        u1=SET,
        u2=UNSET,
        u3=UNSET,
        ur1=SET,
        ur2=SET,
        ur3=SET,
    )

    ### History and field output requests ###
    del ScratchModel.fieldOutputRequests["F-Output-1"]
    del ScratchModel.historyOutputRequests["H-Output-1"]

    ScratchModel.HistoryOutputRequest(
        createStepName=StepName1,
        name="ReactionForces",
        rebar=EXCLUDE,
        region=ScratchModelAssembly.sets[indenter_set],
        sectionPoints=DEFAULT,
        timeInterval=sample_force_frequency_scratching,
        variables=("RF1", "RF2", "RF3"),
    )

    ScratchModel.FieldOutputRequest(
        createStepName=StepName1,
        name="FieldOutput",
        timeInterval=sample_frequency_scratching,
        variables=PRESELECT,
    )
    # ScratchModel.FieldOutputRequest(
    #     createStepName=StepName1,
    #     name="ContactForce",
    #     timeInterval=sample_frequency_scratching,
    #     variables=("CFORCE",),
    # )

    #### ------------------------------ ####
    #           Unloading Step
    #### ------------------------------ ####
    StepName2 = "UnloadingStep"
    ScratchModel.ExplicitDynamicsStep(
        improvedDtMethod=ON,
        name=StepName2,
        previous=StepName1,
        timePeriod=indentation_time,
    )

    ScratchModel.boundaryConditions["IndenterConstraint"].setValuesInStep(
        stepName=StepName2, u3=SET
    )

    ScratchModel.fieldOutputRequests["FieldOutput"].setValuesInStep(
        stepName=StepName2, timeInterval=sample_frequency_indentation
    )
    # ScratchModel.fieldOutputRequests["ContactForce"].setValuesInStep(
    #     stepName=StepName2, timeInterval=sample_frequency_indentation
    # )

    ScratchModel.historyOutputRequests["ReactionForces"].setValuesInStep(
        stepName=StepName2, timeInterval=sample_frequency_indentation
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

    return (ScratchModel, SubstratePart, SubstrateSet)
