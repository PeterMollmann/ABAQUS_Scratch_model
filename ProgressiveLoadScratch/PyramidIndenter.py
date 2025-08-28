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


def PyramidIndenter(
    Model,
    xs2,
    ys2,
    meshMaxSize,
    meshMinSize,
    sheet_size,
):
    IndenterName = "PyramidIndenter"
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
    indenterHeight = (
        xi2 * z_factor
    )  # is the sum of 90 and the angle as the angel has to be negative

    #### ------------------------------ ####
    #          Indenter geometry
    #### ------------------------------ ####
    Model.ConstrainedSketch(name="__profile__", sheetSize=sheet_size)
    Sketch = Model.sketches["__profile__"]
    Sketch.rectangle(
        point1=(xi1, yi1),
        point2=(xi2, yi2),
    )
    Model.Part(dimensionality=THREE_D, name=IndenterName, type=DISCRETE_RIGID_SURFACE)
    Model.parts[IndenterName].BaseSolidExtrude(
        depth=1.0, draftAngle=indenter_angle, sketch=Sketch
    )
    del Sketch
    IndenterPart = Model.parts[IndenterName]
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
            ((xi1 + xi2) / 2.0, (yi1 + yi2) / 2.0, indenterHeight),
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
    #             Meshing
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
        maxSize=meshMaxSize,
        minSize=meshMinSize,
    )
    IndenterPart.generateMesh()

    return Model, IndenterPart, indenter_set, IndenterName, indenterHeight
