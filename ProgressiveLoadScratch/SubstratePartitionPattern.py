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


def FullPartitionOfFace(
    SubstratePart,
    xs1,
    xs2,
    ys1,
    ys2,
    zs1,
    zs2,
    dpo_x,
    dpo_z,
    partitionAreaLength,
    meshSizeX,
    meshSizeZ,
):
    """
    Adds hierarchical partitions (vertical, horizontal, roof-style)
    to an existing substrate part for scratch testing.

    Args:
        SubstratePart: Abaqus Part object (already extruded substrate).
        xs1 (float): x coordinate of one corner of the substrate rectangle.
        ys1 (float): y coordinate of one corner of the substrate rectangle.
        xs2 (float): x coordinate of the opposite corner of the substrate rectangle.
        ys2 (float): y coordinate of the opposite corner of the substrate rectangle.
        zs1 (float): The extrusion depth start value for the substrate.
        zs2 (float): The extrusion depth end value for the substrate.
        dpo_x (float): Datum plane offset in the x direction for partitioning and meshing.
        dpo_z (float): Datum plane offset in the z direction for partitioning and meshing.
        partitionAreaLength (float): length (in the direction of desired partition) of the area that is to be partitioned
        meshSizeX (float): mesh size in the x-direction
        meshSizeZ (float): mesh size in the z-direction

    Returns:
        None (modifies SubstratePart in-place).
    """

    top_face = SubstratePart.faces.getByBoundingBox(
        xs1 + dpo_x, ys2 - 1e-3, zs1 + dpo_z, xs2, ys2 + 1e-3, zs2 - dpo_z
    )

    sk = SubstratePart.MakeSketchTransform(
        # sketchPlane=SubstratePart.faces.findAt(((0.5, 1.0, 2.0),)),
        sketchPlane=top_face[0],
        # sketchPlane=SubstratePart.faces.findAt(
        #     (((xs2 + (xs1 + dpo_x)) / 2.0, ys2, (zs2 + zs1) / 2.0),)
        # ),
        sketchPlaneSide=SIDE1,
        sketchUpEdge=SubstratePart.edges.findAt(
            ((xs2 + (xs1 + dpo_x)) / 2.0, ys2, zs2 - dpo_z),
        ),
        origin=(0.0, ys2, 0.0),
        sketchOrientation=RIGHT,
    )

    s = mdb.models[SubstratePart.modelName].ConstrainedSketch(
        name="__profile__",
        sheetSize=10,
        transform=sk,
    )
    SubstratePart.projectReferencesOntoSketch(
        filter=COPLANAR_EDGES, sketch=mdb.models["Model-1"].sketches["__profile__"]
    )

    # n_partitions = 1
    n_partitions = int((partitionAreaLength / meshSizeZ) / 2)

    for n in range(n_partitions):
        x1 = xs1 + dpo_x
        x2 = x1 + 2 * meshSizeX
        z1 = zs1 + dpo_z + n * 2 * meshSizeZ
        z2 = z1 + 2 * meshSizeZ
        if n % 2 == 0:
            SketchBlock(s=s, x1=x1, x2=x2, z1=z1, z2=z2)
        else:
            SketchBlock(s=s, x1=x1, x2=x2, z1=z2, z2=z1)

    SubstratePart.PartitionCellBySketch(
        cells=SubstratePart.cells.findAt(
            (((xs2 + (xs1 + dpo_x)) / 2.0, ys2, (zs2 + zs1) / 2.0),)
        ),
        sketch=s,
        sketchPlane=top_face[0],
        sketchUpEdge=SubstratePart.edges.findAt(
            ((xs2 + (xs1 + dpo_x)) / 2.0, ys2, zs2 - dpo_z),
        ),
    )
    del s


def SketchBlock(s, x1, x2, z1, z2):
    """
    Sketches one block.
    x is the y-dimension in the sketch, and z is the x-dimension in the sketch

    Args:
        s: Abaqus sketch object.
        x1, x2 (float): x-bounds of segment.
        z1, z2 (float): z-bounds of segment.
        ys2 (float): The y-coordinate of the top surface of the substrate.
    """

    # Middle point
    x_mid = (x2 + x1) / 2.0
    z_mid = (z2 + z1) / 2.0

    # rectangle outline
    s.rectangle(point1=(z1, x1), point2=(z2, x2))

    # Slanted line from corner to middle
    s.Line(point1=(z1, x2), point2=(z_mid, x_mid))

    # Vertical line from bottom to middle
    s.Line(point1=(z_mid, x1), point2=(z_mid, x_mid))

    # Horizontal line from middle to side
    s.Line(point1=(z_mid, x_mid), point2=(z2, x_mid))
