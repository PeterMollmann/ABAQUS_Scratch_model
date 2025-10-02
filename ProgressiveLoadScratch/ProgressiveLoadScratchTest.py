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
from .RockwellIndenter import RockwellIndenter
from .SubstrateGeneration import SubstrateGeneration, SubstrateMeshing
from . import Constants as C

# from SubstratePartitionPattern import FullPartitionOfFace


def ScratchModelSetup(
    SubstrateSizeY=0.020,  # Very coarse mesh for fast simulations
    SubstrateSizeX=0.020,
    SubstrateSizeZ=0.020,
    target_time_increment=0.0,  # if zero, applies mass scaling factor only
    mass_scale=1e4,
):
    """
    Sets up the scratch model with substrate and indenter parts, assembly,
    steps, boundary conditions, loading conditions, contact interactions, and output requests.

    Call once to setup the model.
    Assigning of the material properties is done separately to not construct the model multiple times.
    Run the job after assigning the material properties.

    Remember to use consistent units as abaqus does not check for unit consistency! This document uses SI(mm) units.

    Args:
        SubstrateSizeX: The mesh size of the refined area in the x-direction
        SubstrateSizeY: The mesh size of the refined area in the y-direction
        SubstrateSizeZ: The mesh size of the refined area in the z-direction
        target_time_increment: Target stable time increment for variable mass scaling. If 0.0, variable mass scaling is not used.
        mass_scale: Fixed mass scaling factor. Is not used if target_time_increment is not 0.0.

    Returns:
        ScratchModel: The Abaqus model object with the complete scratch test setup.
        SubstratePart: The Abaqus part object of the substrate.
    """

    # Set the replay options
    session.journalOptions.setValues(
        replayGeometry=COORDINATE, recoverGeometry=COORDINATE
    )

    ScratchModel = mdb.models["Model-1"]

    # Create and mesh substate
    SubstratePart = SubstrateGeneration(
        ScratchModel,
    )

    SubstrateMeshing(
        SubstratePart,
        SubstrateSizeX,
        SubstrateSizeY,
        SubstrateSizeZ,
    )

    # Create and mesh indenter
    IndenterPart = RockwellIndenter(ScratchModel)

    #### ------------------------------ ####
    #            Assembly
    #### ------------------------------ ####
    # indenterInstanceName = IndenterName + "Inst"
    # C.substrate_instance_name = "SubstrateInst"
    ScratchModelAssembly = ScratchModel.rootAssembly
    ScratchModelAssembly.DatumCsysByDefault(CARTESIAN)
    ScratchModelAssembly.Instance(
        dependent=ON, name=C.indenter_instance_name, part=IndenterPart
    )
    IndenterInstance = ScratchModelAssembly.instances[C.indenter_instance_name]
    ScratchModelAssembly.Instance(
        dependent=ON, name=C.substrate_instance_name, part=SubstratePart
    )
    SubstrateInstance = ScratchModelAssembly.instances[C.substrate_instance_name]

    eps = 0.00
    ScratchModelAssembly.translate(
        instanceList=(C.indenter_instance_name,), vector=(0.0, C.ys2 + eps, 0.0)
    )

    ScratchModelAssembly.translate(
        instanceList=(C.indenter_instance_name,),
        vector=(0.0, 0.0, C.dpo_z),
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
                AT_BEGINNING if target_time_increment == 0.0 else THROUGHOUT_STEP,
                mass_scale if target_time_increment == 0.0 else 0.0,
                target_time_increment,
                BELOW_MIN if target_time_increment != 0.0 else None,
                0,
                10,  # Update mass scaling every # equally spaced intervals
                0.0,
                0.0,
                0,
                None,
            ),
        ),
        name=StepName1,
        previous="Initial",
        timePeriod=C.scratch_time,
        nlgeom=ON,
    )

    # Fixed constraint on the bottom two faces
    fixedBCSet = "FIXEDBCSET"
    ScratchModelAssembly.Set(
        faces=SubstrateInstance.faces.findAt(
            ((C.xs1 + C.dpo_x / 2.0, C.ys1, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs1 + C.dpo_x / 2.0, C.ys1, (C.zs2 + C.zs1) / 2.0),),
            ((C.xs1 + C.dpo_x / 2.0, C.ys1, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs2 - C.dpo_x / 2.0, C.ys1, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs2 - C.dpo_x / 2.0, C.ys1, (C.zs2 + C.zs1) / 2.0),),
            ((C.xs2 - C.dpo_x / 2.0, C.ys1, C.zs1 + C.dpo_z / 2.0),),
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
    XsymmetryBCSet = "XsymmetryBCSet"
    ScratchModelAssembly.Set(
        faces=SubstrateInstance.faces.findAt(
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs1, C.ys1 + C.dpo_y / 2.0, (C.zs2 + C.zs1) / 2.0),),
            ((C.xs1, C.ys2 - C.dpo_y / 2.0, (C.zs2 + C.zs1) / 2.0),),
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs2 - C.dpo_z / 2.0),),
        ),
        name=XsymmetryBCSet,
    )
    ScratchModel.XsymmBC(
        createStepName=StepName1,
        localCsys=None,
        name="x_axis_symmetry",
        region=ScratchModelAssembly.sets[XsymmetryBCSet],
    )

    # Z axis symmetry constraint
    ZsymmetryBCSet = "ZsymmetryBCSet"
    ScratchModelAssembly.Set(
        faces=SubstrateInstance.faces.findAt(
            ((C.xs1 + C.dpo_x / 2.0, (C.ys2 + C.ys1) / 2.0, C.zs1),),
            ((C.xs2 - C.dpo_x / 2.0, (C.ys2 + C.ys1) / 2.0, C.zs1),),
            ((C.xs1 + C.dpo_x / 2.0, (C.ys2 + C.ys1) / 2.0, C.zs2),),
            ((C.xs2 - C.dpo_x / 2.0, (C.ys2 + C.ys1) / 2.0, C.zs2),),
        ),
        name=ZsymmetryBCSet,
    )
    ScratchModel.ZsymmBC(
        createStepName=StepName1,
        localCsys=None,
        name="z_axis_symmetry",
        region=ScratchModelAssembly.sets[ZsymmetryBCSet],
    )

    ScratchModel.TabularAmplitude(
        data=(
            (0.0, 0.0),
            (C.scratch_time, 1.0),
            (C.unload_time + C.scratch_time, 0.0),
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
        region=IndenterInstance.sets[C.indenter_set_name],
        u1=UNSET,
        u2=C.scratch_depth,
        u3=C.scratch_length,
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
        region=IndenterInstance.sets[C.indenter_set_name],
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
        region=IndenterInstance.sets[C.indenter_set_name],
        sectionPoints=DEFAULT,
        timeInterval=C.sample_force_frequency,
        variables=("RF1", "RF2", "RF3"),
    )
    ScratchModel.HistoryOutputRequest(
        createStepName=StepName1,
        name="Energy",
        region=SubstrateInstance.sets["SubstrateSet"],
        timeInterval=C.sample_force_frequency,
        variables=("ALLKE", "ALLIE"),
    )

    ScratchModel.FieldOutputRequest(
        createStepName=StepName1,
        name="FieldOutput",
        region=SubstrateInstance.sets["SubstrateSet"],
        timeInterval=C.sample_frequency_scratching,
        variables=(
            "MISES",
            "TRIAX",
            # "A",
            "CSTRESS",
            "EVF",
            "LE",
            "PE",
            "PEEQ",
            # "PEEQVAVG",
            # "PEVAVG",
            "RF",
            "S",
            "SVAVG",
            "U",
            # "V",
            # "DUCTCRT",
            # "JCCRT",
            "DAMAGEC",
            "DAMAGET",
            "DMICRT",
            "STATUS",
            "SDEG",
            "EMSF",
            # "DMASS",
        ),
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
        timePeriod=C.unload_time,
    )

    ScratchModel.boundaryConditions["IndenterConstraint"].setValuesInStep(
        stepName=StepName2, u3=SET
    )

    ScratchModel.fieldOutputRequests["FieldOutput"].setValuesInStep(
        stepName=StepName2, timeInterval=C.sample_frequency_unloading
    )

    # ScratchModel.fieldOutputRequests["ContactForce"].setValuesInStep(
    #     stepName=StepName2, timeInterval=sample_frequency_indentation
    # )

    ScratchModel.historyOutputRequests["ReactionForces"].setValuesInStep(
        stepName=StepName2, timeInterval=C.sample_frequency_unloading
    )

    ScratchModel.historyOutputRequests["Energy"].setValuesInStep(
        stepName=StepName2, timeInterval=C.sample_frequency_unloading
    )

    #### ------------------------------ ####
    #          Contact modelling
    #### ------------------------------ ####
    # Initiate tangential behaviour with zero friction. Friction is updated in another file.
    ScratchModel.ContactProperty("IntProp-1")
    ScratchModel.interactionProperties["IntProp-1"].TangentialBehavior(
        formulation=PENALTY, table=((0.0,),), fraction=0.005
    )
    mdb.models["Model-1"].interactionProperties["IntProp-1"].NormalBehavior(
        allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT,
        pressureOverclosure=HARD,
    )
    # ScratchModel.interactionProperties["IntProp-1"].NormalBehavior(
    #     constraintEnforcementMethod=DEFAULT,
    #     contactStiffnessScaleFactor=2.0,
    #     initialStiffnessScaleFactor=1.0,
    #     overclosureFactor=1.0,
    #     overclosureMeasure=0.0,
    #     pressureOverclosure=SCALE_FACTOR,
    # )

    ScratchModelAssembly.Surface(
        name=C.master_surface_name,
        side1Faces=IndenterInstance.faces.findAt(
            ((C.xs1, C.ys2 + eps, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.xl2, C.ys2 + C.yl2 + eps, C.zs1 + C.dpo_z),),
        ),
    )

    # Get all relevant contact elements of the substrate (+- 0.1 for ensureing all elements are selected)
    elemsInBox = SubstrateInstance.elements.getByBoundingBox(
        C.xs1 - 0.1,
        (C.ys2 - C.dpo_y) - 0.1,
        C.zs1 + C.dpo_z - 0.1,
        C.xs1 + C.dpo_x + 0.1,
        C.ys2 + 0.1,
        (C.zs2 - C.dpo_z) + 0.1,
    )

    ScratchModelAssembly.Surface(
        name=C.slave_surface_name,
        face1Elements=elemsInBox,
        face2Elements=elemsInBox,
        face3Elements=elemsInBox,
        face4Elements=elemsInBox,
        face5Elements=elemsInBox,
        face6Elements=elemsInBox,
    )

    # Use abaqus general contact
    ScratchModel.ContactExp(createStepName="Initial", name="Int-1")
    ScratchModel.interactions["Int-1"].includedPairs.setValuesInStep(
        addPairs=(
            (
                ScratchModelAssembly.surfaces[C.master_surface_name],
                ScratchModelAssembly.surfaces[C.slave_surface_name],
            ),
        ),
        stepName="Initial",
        useAllstar=OFF,
    )
    ScratchModel.interactions["Int-1"].contactPropertyAssignments.appendInStep(
        assignments=((GLOBAL, SELF, "IntProp-1"),), stepName="Initial"
    )

    ScratchModelAssembly.Set(
        name=C.contact_region_nodes_name,
        nodes=ScratchModelAssembly.allSurfaces[C.slave_surface_name].nodes,
    )

    # ScratchModel.AdaptiveMeshControl(name="Ada-1")
    # ScratchModel.DisplacementAdaptiveMeshConstraint(
    #     amplitude=UNSET,
    #     createStepName=StepName1,
    #     localCsys=None,
    #     motionType=FOLLOW,
    #     name="Ada-Cons-1",
    #     region=SubstrateInstance.faces.findAt(
    #         (C.xs1, C.ys2 - C.dpo_y / 2.0, (C.zs2 + C.zs1) / 2.0),
    #     ),
    #     u1=0.0,
    #     u2=UNSET,
    #     u3=UNSET,
    #     ur1=UNSET,
    #     ur2=UNSET,
    #     ur3=UNSET,
    # )
    # ScratchModel.steps[StepName1].AdaptiveMeshDomain(
    #     controls="Ada-1",
    #     meshSweeps=5,
    #     frequency=20,
    #     initialMeshSweeps=5,
    #     region=Region(
    #         cells=SubstrateInstance.cells.findAt(
    #             ((C.xs1, C.ys1, C.zs1),),
    #             ((C.xs2, C.ys1, C.zs1),),
    #             ((C.xs1, C.ys1, C.zs2),),
    #             ((C.xs2, C.ys1, C.zs2),),
    #             ((C.xs1, C.ys1, (C.zs1 + C.zs2) / 2.0),),
    #             ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
    #             ((C.xs2, C.ys1, (C.zs1 + C.zs2) / 2.0),),
    #         )
    #     ),
    # )
    # ScratchModel.steps[StepName2].AdaptiveMeshDomain(
    #     controls="Ada-1",
    #     meshSweeps=5,
    #     frequency=20,
    #     initialMeshSweeps=5,
    #     region=Region(
    #         cells=SubstrateInstance.cells.findAt(
    #             ((C.xs1, C.ys1, C.zs1),),
    #             ((C.xs2, C.ys1, C.zs1),),
    #             ((C.xs1, C.ys1, C.zs2),),
    #             ((C.xs2, C.ys1, C.zs2),),
    #             ((C.xs1, C.ys1, (C.zs1 + C.zs2) / 2.0),),
    #             ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
    #             ((C.xs2, C.ys1, (C.zs1 + C.zs2) / 2.0),),
    #         )
    #     ),
    # )

    return ScratchModel, SubstratePart
