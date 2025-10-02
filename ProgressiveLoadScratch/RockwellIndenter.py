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


def RockwellIndenter(Model):
    """
    Makes the Rockwell indenter part in the given Model. First a 2D sketch is made and then revolved to create a 3D part.
    The indenter is modelled as a rigid body with a reference point (and mass) at the center of the spherical tip.
    The indenter is meshed with a sweep mesh.

    Args:
        Model: The Abaqus model object where the indenter will be created.

    Returns:
        Model: The Abaqus model object with the indenter added.
        IndenterPart: The created indenter part.
    """

    # IndenterName = "RockwellIndenter"
    Model.ConstrainedSketch(name="__profile__", sheetSize=C.sheet_size)
    Sketch = Model.sketches["__profile__"]

    # Geometry coordinates for indenter
    # xc1 = 0.0  # x coordinate of center of spherical tip
    # yc1 = 0.0  # y coordinate of center of spherical tip
    # xc2 = R * cos(-theta * pi / 180)  # x coordinate of sphere end point
    # yc2 = R + R * sin(-theta * pi / 180)  # y coordinate of sphere end point
    # xc3 = R * cos(
    #     (-theta - (90 - theta) / 2.0) * pi / 180
    # )  # x coordinate of shpere mid point
    # yc3 = R + R * sin(
    #     (-theta - (90 - theta) / 2.0) * pi / 180
    # )  # y coordinate of sphere mid point

    # xl1 = xc2  # line x coorcinate 1
    # yl1 = yc2  # line y coordinate 1

    # xl2 = xl1 + 0.5 * cos((90 - theta) * pi / 180)  # line x coordinate 2
    # yl2 = yl1 + 0.5 * sin((90 - theta) * pi / 180)  # line y coordinate 2

    # Make sketch
    Sketch.ConstructionLine(point1=(0.0, -5.0), point2=(0.0, 5.0))
    Sketch.FixedConstraint(
        entity=Sketch.geometry.findAt(
            (0.0, 0.0),
        )
    )

    # Make construction line for horizontal axis
    Sketch.ConstructionLine(point1=(0.0, 0.0), point2=(1.0, 0.0))
    Sketch.HorizontalConstraint(
        addUndoState=False,
        entity=Sketch.geometry.findAt(
            (0.5, 0.0),
        ),
    )
    Sketch.FixedConstraint(
        entity=Sketch.geometry.findAt(
            (1.0, 0.0),
        )
    )

    # make arc for spherical tip
    Sketch.ArcByCenterEnds(
        center=(0.0, C.tip_radius), point1=(C.xc1, C.yc1), point2=(C.xc2, C.yc2)
    )
    Sketch.CoincidentConstraint(
        entity1=Sketch.vertices.findAt(
            (C.xc1, C.yc1),
        ),
        entity2=Sketch.geometry.findAt(
            (0.5, 0.0),
        ),
    )
    Sketch.CoincidentConstraint(
        entity1=Sketch.vertices.findAt(
            (C.xc1, C.yc1),
        ),
        entity2=Sketch.geometry.findAt(
            (0.0, 1.0),
        ),
    )

    # make lines for conical part
    Sketch.Line(point1=(C.xl1, C.yl1), point2=(C.xl2, C.yl2))

    # Make conical line tangent to spherical tip
    Sketch.TangentConstraint(
        entity1=Sketch.geometry.findAt(
            (C.xl2, C.yl2),
        ),
        entity2=Sketch.geometry.findAt(
            (C.xc3, C.yc3),
        ),
    )

    # Make conical line coincident with spherical tip
    Sketch.CoincidentConstraint(
        entity1=Sketch.vertices.findAt(
            (C.xl1, C.yl1),
        ),
        entity2=Sketch.vertices.findAt(
            (C.xc2, C.yc2),
        ),
    )

    # Revolve to make 3D part
    Sketch.sketchOptions.setValues(constructionGeometry=ON)
    Sketch.assignCenterline(
        line=Sketch.geometry.findAt(
            (0.0, 1.0),
        )
    )
    Model.Part(
        dimensionality=THREE_D, name=C.indenter_name, type=ANALYTIC_RIGID_SURFACE
    )
    Model.parts[C.indenter_name].AnalyticRigidSurfRevolve(
        sketch=Model.sketches["__profile__"]
    )
    # Model.parts[IndenterName].BaseShellRevolve(
    #     angle=360.0,
    #     flipRevolveDirection=OFF,
    #     sketch=Sketch,
    # )

    # Sketch.Line(point1=(xl2, yl2), point2=(0.0, yl2))
    # Sketch.Line(point1=(0.0, yl2), point2=(0.0, 0.0))
    # Model.Part(dimensionality=THREE_D, name=IndenterName, type=DEFORMABLE_BODY)
    # Model.parts[IndenterName].BaseSolidRevolve(
    #     angle=360.0,
    #     flipRevolveDirection=OFF,
    #     sketch=Sketch,
    # )
    del Sketch
    IndenterPart = Model.parts[C.indenter_name]

    #### ------------------------------ ####
    #         Materials Indenter
    #### ------------------------------ ####
    # Creating reference point and inertia for indenter
    # indenter_set = C.indenter_name + "Set"
    IndenterPart.ReferencePoint(
        point=IndenterPart.vertices.findAt(
            (C.xc1, C.yc1, 0.0),
        )
    )
    IndenterPart.Set(
        name=C.indenter_set_name, referencePoints=(IndenterPart.referencePoints[2],)
    )

    IndenterPart.engineeringFeatures.PointMassInertia(
        alpha=0.0,
        composite=0.0,
        i11=0.0,
        i22=0.0,
        i33=0.0,
        mass=1.0,
        name="IndenterInertia",
        region=IndenterPart.sets[C.indenter_set_name],
    )

    # mat = Model.Material(name="IndenterMaterial")
    # mat.Density(table=((3.5e-3,),))
    # mat.Elastic(table=((1050e3, 0.2),))
    # Model.HomogeneousSolidSection(
    #     material="IndenterMaterial", name="IndenterSection", thickness=None
    # )

    # IndenterPart.SectionAssignment(
    #     offset=0.0,
    #     offsetField="",
    #     offsetType=MIDDLE_SURFACE,
    #     region=IndenterPart.sets[indenter_set],
    #     sectionName="IndenterSection",
    #     thicknessAssignment=FROM_SECTION,
    # )

    #### ------------------------------ ####
    #             Meshing
    #### ------------------------------ ####
    # IndenterPart.setMeshControls(
    #     regions=IndenterPart.faces.findAt(
    #         ((xc1, yc1, 0.0),),
    #         ((xl2, yl2, 0.0),),
    #     ),
    #     technique=SWEEP,
    # )
    # IndenterPart.seedEdgeByBias(
    #     biasMethod=SINGLE,
    #     constraint=FINER,
    #     end2Edges=IndenterPart.edges.findAt(((xc1, yc1, 0.0),)),
    #     maxSize=meshMaxSize,
    #     minSize=meshMinSize,
    # )

    # IndenterPart.seedEdgeBySize(
    #     constraint=FINER,
    #     deviationFactor=0.1,
    #     edges=IndenterPart.edges.findAt((((xl1 + xl2) / 2.0, (yl1 + yl2) / 2.0, 0.0),)),
    #     size=meshMaxSize,
    # )

    # IndenterPart.generateMesh()

    return IndenterPart
