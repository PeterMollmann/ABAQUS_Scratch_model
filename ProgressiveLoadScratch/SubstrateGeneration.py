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


def SubstrateGeneration(
    ScratchModel,
    xs1,
    ys1,
    xs2,
    ys2,
    zs1,
    zs2,
    dpo_x,
    dpo_z,
    sheet_size,
):
    """
    Makes the substrate part in the given Model. First a 2D sketch is made and then extruded to create a 3D part.
    The substrate is modelled as a deformable body. Several partitions are made to optimise the mesh.


    Args:
        ScratchModel: The Abaqus model object where the substrate will be created.
        xs1 (float): x coordinate of one corner of the substrate rectangle.
        ys1 (float): y coordinate of one corner of the substrate rectangle.
        xs2 (float): x coordinate of the opposite corner of the substrate rectangle.
        ys2 (float): y coordinate of the opposite corner of the substrate rectangle.
        zs1 (float): The extrusion depth start value for the substrate.
        zs2 (float): The extrusion depth end value for the substrate.
        dpo_x (float): Datum plane offset in the x direction for partitioning and meshing.
        dpo_z (float): Datum plane offset in the z direction for partitioning and meshing.
        sheet_size (float): Size of the sketching sheet.

    Returns:
        ScratchModel: The Abaqus model object with the substrate added.
        SubstratePart: The created substrate part.
        SubstrateName: The name of the substrate part.
        SubstrateSet: The name of the set containing all the substrate cells.
    """

    SubstrateName = "Substrate"
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
    SubstratePart.DatumPlaneByPrincipalPlane(offset=dpo_x, principalPlane=YZPLANE)
    SubstratePart.DatumPlaneByPrincipalPlane(offset=zs1 + dpo_z, principalPlane=XYPLANE)
    SubstratePart.DatumPlaneByPrincipalPlane(offset=zs2 - dpo_z, principalPlane=XYPLANE)

    # Make partition of substrate for mesh refinement along scratch path
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((xs2, ys2, zs2),),
        ),
        datumPlane=SubstratePart.datums[2],
    )  # partition along z axis

    # partitions along x axis for mesh optimization
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs2, ys1, zs1),),
        ),
        datumPlane=SubstratePart.datums[3],
    )
    SubstratePart.PartitionCellByDatumPlane(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs2),),
            ((xs2, ys1, zs2),),
        ),
        datumPlane=SubstratePart.datums[4],
    )

    SubstrateSet = SubstrateName + "Set"
    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs2, ys1, zs1),),
            ((xs1, ys1, zs2),),
            ((xs2, ys1, zs2),),
            ((xs1, ys1, (zs1 + zs2) / 2.0),),
            ((xs2, ys1, (zs1 + zs2) / 2.0),),
        ),
        name=SubstrateSet,
    )

    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys2, (zs1 + zs2) / 2.0),),
        ),
        name="RefinedArea",
    )
    return ScratchModel, SubstratePart, SubstrateName, SubstrateSet


def SubstrateMeshing(
    SubstratePart,
    xs1,
    ys1,
    zs1,
    xs2,
    ys2,
    zs2,
    dpo_x,
    dpo_z,
    CoarseMeshSize0,
    CoarseMeshSize1,
    CoarseMeshSize2,
    SubstrateSizeX,
    SubstrateSizeY,
    SubstrateSizeZ,
):
    """
    The substrate is meshing with a structured hex mesh.

    Args:
        SubstratePart: The Abaqus part object of the substrate to be meshed.
        xs1 (float): x coordinate of one corner of the substrate rectangle.
        ys1 (float): y coordinate of one corner of the substrate rectangle.
        xs2 (float): x coordinate of the opposite corner of the substrate rectangle.
        ys2 (float): y coordinate of the opposite corner of the substrate rectangle.
        zs1 (float): The extrusion depth start value for the substrate.
        zs2 (float): The extrusion depth end value for the substrate.
        dpo_x (float): Datum plane offset in the x direction for partitioning and meshing.
        dpo_z (float): Datum plane offset in the z direction for partitioning and meshing.


    Returns:
        SubstratePart: The Abaqus part object of the substrate with the mesh.
    """

    SubstratePart.setMeshControls(
        elemShape=HEX,
        regions=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs2, ys1, zs1),),
            ((xs1, ys1, zs2),),
            ((xs2, ys1, zs2),),
            ((xs1, ys1, (zs1 + zs2) / 2.0),),
            ((xs2, ys1, (zs1 + zs2) / 2.0),),
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
                hourglassControl=ENHANCED,
            ),
            ElemType(elemCode=UNKNOWN_WEDGE, elemLibrary=EXPLICIT),
            ElemType(
                elemCode=UNKNOWN_TET,
                elemLibrary=EXPLICIT,
            ),
        ),
        regions=(
            SubstratePart.cells.findAt(
                ((xs1, ys1, zs1),),
                ((xs2, ys1, zs1),),
                ((xs1, ys1, zs2),),
                ((xs2, ys1, zs2),),
                ((xs1, ys1, (zs1 + zs2) / 2.0),),
                ((xs2, ys1, (zs1 + zs2) / 2.0),),
            ),
        ),
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys2, (zs1 + zs2) / 2.0),),
            ((xs1 + dpo_x, ys2, (zs1 + zs2) / 2.0),),
        ),
        size=SubstrateSizeZ,
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1 + dpo_x / 2.0, ys2, zs1 + dpo_z),),
            ((xs1 + dpo_x / 2.0, ys2, zs2 - dpo_z),),
        ),
        size=SubstrateSizeX,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((xs1, (ys1 + ys2) / 2.0, zs1 + dpo_z),),
            ((xs1, (ys1 + ys2) / 2.0, zs2 - dpo_z),),
            ((xs2, (ys1 + ys2) / 2.0, zs1),),
            ((xs2, (ys1 + ys2) / 2.0, zs2),),
            ((xs1 + dpo_x, (ys1 + ys2) / 2.0, zs2),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((xs2, (ys1 + ys2) / 2.0, zs1 + dpo_z),),
            ((xs2, (ys1 + ys2) / 2.0, zs2 - dpo_z),),
            ((xs1, (ys1 + ys2) / 2.0, zs1),),
            ((xs1, (ys1 + ys2) / 2.0, zs2),),
            ((xs1 + dpo_x, (ys1 + ys2) / 2.0, zs1),),
        ),
        maxSize=CoarseMeshSize2,
        minSize=SubstrateSizeY,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((xs1 + dpo_x, ys2, zs2 - dpo_z / 2.0),),
            ((xs1, ys2, zs2 - dpo_z / 2.0),),
            ((xs2, ys2, zs2 - dpo_z / 2.0),),
            ((xs1, ys1, zs2 - dpo_z / 2.0),),
            ((xs1 + dpo_x, ys1, zs1 + dpo_z / 2.0),),
            ((xs2, ys1, zs2 - dpo_z / 2.0),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((xs1 + dpo_x, ys2, zs1 + dpo_z / 2.0),),
            ((xs1, ys2, zs1 + dpo_z / 2.0),),
            ((xs2, ys2, zs1 + dpo_z / 2.0),),
            ((xs1, ys1, zs1 + dpo_z / 2.0),),
            ((xs1 + dpo_x, ys1, zs2 - dpo_z / 2.0),),
            ((xs2, ys1, zs1 + dpo_z / 2.0),),
        ),
        maxSize=CoarseMeshSize1,
        minSize=SubstrateSizeZ,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs1),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs2),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs1 + dpo_z),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs2 - dpo_z),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs1 + dpo_z),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs2 - dpo_z),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs1),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs2),),
        ),
        maxSize=CoarseMeshSize2,
        minSize=SubstrateSizeX,
    )

    SubstratePart.generateMesh()

    return SubstratePart
