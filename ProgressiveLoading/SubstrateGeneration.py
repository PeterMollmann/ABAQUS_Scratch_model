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
    dpo_y,
    sheet_size,
):
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

    SubstrateSet = SubstrateName + "Set"
    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs1, ys2, zs1),),
            ((xs2, ys2, zs1),),
            ((xs2, ys1, zs1),),
        ),
        name=SubstrateSet,
    )

    SubstratePart.Set(
        cells=SubstratePart.cells.findAt(
            ((xs1, ys2, zs1),),
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
    dpo_y,
    CoarseMeshSize0,
    CoarseMeshSize1,
    CoarseMeshSize2,
    SubstrateSizeX,
    SubstrateSizeY,
    SubstrateSizeZ,
):

    SubstratePart.setMeshControls(
        elemShape=HEX,
        regions=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs1, ys2, zs1),),
            ((xs2, ys1, zs1),),
            ((xs2, ys2, zs1),),
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
                ((xs1, ys2, zs1),),
                ((xs2, ys1, zs1),),
                ((xs2, ys2, zs1),),
            ),
        ),
    )

    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys2, zs2 / 2.0),),
            ((xs1, ys2 - dpo_y, zs2 / 2.0),),
            ((xs1 + dpo_x, ys2, zs2 / 2.0),),
            ((xs1 + dpo_x, ys2 - dpo_y, zs2 / 2.0),),
        ),
        size=SubstrateSizeZ,
    )
    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1 + dpo_x / 2.0, ys2, zs1),),
            ((xs1 + dpo_x / 2.0, ys2 - dpo_y, zs1),),
            ((xs1 + dpo_x / 2.0, ys2, zs2),),
            ((xs1 + dpo_x / 2.0, ys2 - dpo_y, zs2),),
        ),
        size=SubstrateSizeX,
    )
    SubstratePart.seedEdgeBySize(
        constraint=FINER,
        deviationFactor=0.1,
        edges=SubstratePart.edges.findAt(
            ((xs1, ys2 - dpo_y / 2.0, zs1),),
            ((xs1, ys2 - dpo_y / 2.0, zs1),),
        ),
        size=SubstrateSizeY,
    )

    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((xs1 + dpo_x, (ys2 - dpo_y) / 2.0, zs2),),
            ((xs2, (ys2 - dpo_y) / 2.0, zs2),),
            ((xs2, (ys2 - dpo_y) / 2.0, zs1),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((xs1, (ys2 - dpo_y) / 2.0, zs1),),
            ((xs1, (ys2 - dpo_y) / 2.0, zs2),),
            ((xs1 + dpo_x, (ys2 - dpo_y) / 2.0, zs1),),
        ),
        maxSize=CoarseMeshSize2,
        minSize=SubstrateSizeX,
    )
    SubstratePart.seedEdgeByBias(
        biasMethod=SINGLE,
        constraint=FINER,
        end1Edges=SubstratePart.edges.findAt(
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs1),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2, zs2),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2 - dpo_y, zs2),),
        ),
        end2Edges=SubstratePart.edges.findAt(
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys2 - dpo_y, zs1),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs1),),
            ((dpo_x + (xs2 - dpo_x) / 2.0, ys1, zs2),),
        ),
        maxSize=CoarseMeshSize2,
        minSize=SubstrateSizeX,
    )

    SubstratePart.generateMesh()

    return SubstratePart
