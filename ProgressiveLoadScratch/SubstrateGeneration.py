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
from . import Constants as C


def SubstrateGeneration(
    ScratchModel,
):
    """
    Makes the substrate part in the given Model. First a 2D sketch is made and then extruded to create a 3D part.
    The substrate is modelled as a deformable body. Several partitions are made to optimise the mesh.


    Args:
        ScratchModel: The Abaqus model object where the substrate will be created.

    Returns:
        ScratchModel: The Abaqus model object with the substrate added.
        SubstratePart: The created substrate part.
    """

    #### ------------------------------ ####
    #         Substrate geometry
    #### ------------------------------ ####
    ScratchModel.ConstrainedSketch(
        name="__profile__", sheetSize=C.sheet_size
    )  # Make sketching sheet
    Sketch = ScratchModel.sketches["__profile__"]
    Sketch.rectangle(
        point1=(C.xs1, C.ys1), point2=(C.xs2, C.ys2)
    )  # Create rectangle based on coordinates defined for the substrate
    ScratchModel.Part(
        dimensionality=THREE_D, name=C.substrate_name, type=DEFORMABLE_BODY
    )  # Substrate is modelled as a deformable body
    ScratchModel.parts[C.substrate_name].BaseSolidExtrude(
        depth=C.zs2, sketch=Sketch
    )  # Extrude sketch by "C.zs2" amount
    del Sketch
    SubstratePart = ScratchModel.parts[C.substrate_name]

    # Create datum plane for partition
    SubstratePart.DatumPlaneByPrincipalPlane(offset=C.dpo_x, principalPlane=YZPLANE)
    SubstratePart.DatumPlaneByPrincipalPlane(
        offset=C.zs1 + C.dpo_z, principalPlane=XYPLANE
    )
    SubstratePart.DatumPlaneByPrincipalPlane(
        offset=C.zs2 - C.dpo_z, principalPlane=XYPLANE
    )
    SubstratePart.DatumPlaneByPrincipalPlane(
        offset=C.ys2 - C.dpo_y, principalPlane=XZPLANE
    )

    # Make partition of substrate for mesh refinement along scratch path
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((C.xs2, C.ys2, C.zs2),),
        ),
        datumPlane=SubstratePart.datums[2],
    )  # partition along z axis

    # partitions along x axis for mesh optimization
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((C.xs1, C.ys1, C.zs1),),
            ((C.xs2, C.ys1, C.zs1),),
        ),
        datumPlane=SubstratePart.datums[3],
    )
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((C.xs1, C.ys1, C.zs2),),
            ((C.xs2, C.ys1, C.zs2),),
        ),
        datumPlane=SubstratePart.datums[4],
    )

    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((C.xs1, C.ys1, (C.zs2 + C.zs1) / 2.0),),
        ),
        datumPlane=SubstratePart.datums[5],
    )

    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((C.xs1, C.ys1, C.zs1),),
            ((C.xs2, C.ys1, C.zs1),),
            ((C.xs1, C.ys1, C.zs2),),
            ((C.xs2, C.ys1, C.zs2),),
            ((C.xs1, C.ys1, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs2, C.ys1, (C.zs1 + C.zs2) / 2.0),),
        ),
        name=C.substrate_set_name,
    )

    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
        ),
        name="RefinedArea",
    )
    return SubstratePart


def SubstrateMeshing(
    SubstratePart,
    SubstrateSizeX,
    SubstrateSizeY,
    SubstrateSizeZ,
):
    """
    The substrate is meshing with a structured hex mesh.

    Args:
        SubstratePart: The Abaqus part object of the substrate to be meshed.

    Returns:
        SubstratePart: The Abaqus part object of the substrate with the mesh.
    """

    SubstratePart.setMeshControls(
        elemShape=HEX,
        regions=SubstratePart.cells.findAt(
            ((C.xs1, C.ys1, C.zs1),),
            ((C.xs2, C.ys1, C.zs1),),
            ((C.xs1, C.ys1, C.zs2),),
            ((C.xs2, C.ys1, C.zs2),),
            ((C.xs1, C.ys1, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs2, C.ys1, (C.zs1 + C.zs2) / 2.0),),
        ),
        technique=STRUCTURED,
    )

    SubstratePart.setElementType(
        elemTypes=(
            ElemType(
                elemCode=C3D8R,
                elemLibrary=EXPLICIT,
                secondOrderAccuracy=OFF,
                distortionControl=DEFAULT,
                hourglassControl=DEFAULT,
                elemDeletion=OFF,
                maxDegradation=0.8,
                # particleConversion=STRAIN,
                # particleConversionThreshold=0.200000002980232,
                # particleConversionPPD=1,
                # particleConversionKernel=CUBIC,
            ),
            ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
            ElemType(
                elemCode=C3D4,
                elemLibrary=EXPLICIT,
            ),
        ),
        regions=(
            SubstratePart.cells.findAt(
                ((C.xs1, C.ys1, C.zs1),),
                ((C.xs2, C.ys1, C.zs1),),
                ((C.xs1, C.ys1, C.zs2),),
                ((C.xs2, C.ys1, C.zs2),),
                ((C.xs1, C.ys1, (C.zs1 + C.zs2) / 2.0),),
                ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
                ((C.xs2, C.ys1, (C.zs1 + C.zs2) / 2.0),),
            ),
        ),
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((C.xs1, C.ys2, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs1, C.ys2 - C.dpo_y, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs1 + C.dpo_x, C.ys2, (C.zs1 + C.zs2) / 2.0),),
            ((C.xs1 + C.dpo_x, C.ys2 - C.dpo_y, (C.zs1 + C.zs2) / 2.0),),
        ),
        size=SubstrateSizeZ,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((C.xs1 + C.dpo_x / 2.0, C.ys2, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.dpo_x / 2.0, C.ys2 - C.dpo_y, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.dpo_x / 2.0, C.ys2, C.zs2 - C.dpo_z),),
            ((C.xs1 + C.dpo_x / 2.0, C.ys2 - C.dpo_y, C.zs2 - C.dpo_z),),
        ),
        size=SubstrateSizeX,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((C.xs1, C.ys2 - C.dpo_y / 2.0, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.dpo_x, C.ys2 - C.dpo_y / 2.0, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.dpo_x, C.ys2 - C.dpo_y / 2.0, C.zs2 - C.dpo_z),),
            ((C.xs1 + C.dpo_x, C.ys2 - C.dpo_y / 2.0, C.zs2 - C.dpo_z),),
        ),
        size=SubstrateSizeY,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs1 + C.dpo_z),),
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs2 - C.dpo_z),),
            ((C.xs2, (C.ys1 + C.ys2) / 2.0, C.zs1),),
            ((C.xs2, (C.ys1 + C.ys2) / 2.0, C.zs2),),
            ((C.xs1 + C.dpo_x, (C.ys1 + C.ys2) / 2.0, C.zs2),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((C.xs2, (C.ys1 + C.ys2) / 2.0, C.zs1 + C.dpo_z),),
            ((C.xs2, (C.ys1 + C.ys2) / 2.0, C.zs2 - C.dpo_z),),
            ((C.xs1 + C.dpo_x, (C.ys1 + C.ys2) / 2.0, C.zs1 + C.dpo_z),),
            ((C.xs1 + C.dpo_x, (C.ys1 + C.ys2) / 2.0, C.zs2 - C.dpo_z),),
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs1),),
            ((C.xs1, (C.ys1 + C.ys2) / 2.0, C.zs2),),
            ((C.xs1 + C.dpo_x, (C.ys1 + C.ys2) / 2.0, C.zs1),),
        ),
        maxSize=C.coarse_mesh_size_2,
        minSize=SubstrateSizeY,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((C.xs1 + C.dpo_x, C.ys2, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs1, C.ys2, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs2, C.ys2, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs1, C.ys1, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs1 + C.dpo_x, C.ys1, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs2, C.ys1, C.zs2 - C.dpo_z / 2.0),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((C.xs1 + C.dpo_x, C.ys2, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs1, C.ys2, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs2, C.ys2, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs1, C.ys1, C.zs1 + C.dpo_z / 2.0),),
            ((C.xs1 + C.dpo_x, C.ys1, C.zs2 - C.dpo_z / 2.0),),
            ((C.xs2, C.ys1, C.zs1 + C.dpo_z / 2.0),),
        ),
        maxSize=C.coarse_mesh_size_1,
        minSize=SubstrateSizeZ,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys2, C.zs1),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys2, C.zs2),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys1, C.zs1 + C.dpo_z),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys1, C.zs2 - C.dpo_z),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys2, C.zs1 + C.dpo_z),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys2, C.zs2 - C.dpo_z),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys1, C.zs1),),
            ((C.dpo_x + (C.xs2 - C.dpo_x) / 2.0, C.ys1, C.zs2),),
        ),
        maxSize=C.coarse_mesh_size_2,
        minSize=SubstrateSizeX,
    )

    SubstratePart.generateMesh()

    # return SubstratePart
