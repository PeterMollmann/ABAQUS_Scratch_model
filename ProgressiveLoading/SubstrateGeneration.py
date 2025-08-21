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
import numpy as np


class Substrate:
    def __init__(self, model, xs1, ys1, xs2, ys2, zs1, zs2, dpo_x, dpo_y, sheet_size):
        self.model = model
        self.xs1, self.ys1, self.xs2, self.ys2 = xs1, ys1, xs2, ys2
        self.zs1, self.zs2 = zs1, zs2
        self.dpo_x, self.dpo_y = dpo_x, dpo_y
        self.sheet_size = sheet_size

        self.name = "Substrate"
        self.part = None

    # ------------------------------
    # Geometry
    # ------------------------------
    def generate_geometry(self):
        self.model.ConstrainedSketch(name="__profile__", sheetSize=self.sheet_size)
        sketch = self.model.sketches["__profile__"]
        sketch.rectangle(point1=(self.xs1, self.ys1), point2=(self.xs2, self.ys2))

        self.model.Part(
            dimensionality=THREE_D,
            name=self.name,
            type=DEFORMABLE_BODY,
        )
        self.model.parts[self.name].BaseSolidExtrude(depth=self.zs2, sketch=sketch)
        del sketch

        self.part = self.model.parts[self.name]

        # Datum planes
        self.part.DatumPlaneByPrincipalPlane(offset=self.dpo_x, principalPlane=YZPLANE)
        self.part.DatumPlaneByPrincipalPlane(
            offset=(self.ys2 - self.dpo_y), principalPlane=XZPLANE
        )

        # Partitions
        self.part.PartitionCellByDatumPlane(
            cells=self.part.cells.findAt(((self.xs2, self.ys2, self.zs2),)),
            datumPlane=self.part.datums[2],
        )
        self.part.PartitionCellByDatumPlane(
            cells=self.part.cells.findAt(
                ((self.xs1, self.ys1, self.zs1),),
                ((self.xs2, self.ys2, self.zs2),),
            ),
            datumPlane=self.part.datums[3],
        )
        return self.part

    # ------------------------------
    # Material + Section assignment
    # ------------------------------
    def assign_material(self, rho, E, nu, yield_strength, n):
        mat_name = "SubstrateMaterial"

        self.model.Material(name=mat_name)
        self.model.materials[mat_name].Density(table=((rho,),))
        self.model.materials[mat_name].Elastic(table=((E, nu),))

        # Plasticity
        total_strain_data = np.append(
            np.arange(0.0001, 0.1, 0.01), np.linspace(0.1, 2, 20)
        )
        plastic_behaviour = []
        yield_strain = yield_strength / E
        K = E * yield_strain ** (1 - n)
        plastic_behaviour.append((yield_strength, 0))
        for total_strain in total_strain_data:
            if total_strain > yield_strain:
                yield_strength_ = round(K * total_strain**n, 5)
                plastic_strain = round(total_strain - yield_strength_ / E, 5)
                plastic_behaviour.append((yield_strength_, plastic_strain))

        self.model.materials[mat_name].Plastic(table=tuple(plastic_behaviour))

        # Section
        self.model.HomogeneousSolidSection(
            material=mat_name, name="SubstrateSection", thickness=None
        )
        set_name = self.name + "Set"
        self.part.Set(
            cells=self.part.cells.findAt(
                ((self.xs1, self.ys1, self.zs1),),
                ((self.xs1, self.ys2, self.zs1),),
                ((self.xs2, self.ys2, self.zs1),),
                ((self.xs2, self.ys1, self.zs1),),
            ),
            name=set_name,
        )
        self.part.SectionAssignment(
            offset=0.0,
            offsetField="",
            offsetType=MIDDLE_SURFACE,
            region=self.part.sets[set_name],
            sectionName="SubstrateSection",
            thicknessAssignment=FROM_SECTION,
        )

        # Refined area set
        self.part.Set(
            cells=self.part.cells.findAt(((self.xs1, self.ys2, self.zs1),)),
            name="RefinedArea",
        )

        return self.part

    # ------------------------------
    # Meshing
    # ------------------------------
    def mesh(
        self,
        CoarseMeshSize0,
        CoarseMeshSize1,
        CoarseMeshSize2,
        SubstrateSizeX,
        SubstrateSizeY,
        SubstrateSizeZ,
    ):

        self.part.setMeshControls(
            elemShape=TET,
            regions=self.part.cells.findAt(
                ((self.xs1, self.ys1, self.zs1),),
                ((self.xs1, self.ys2, self.zs1),),
                ((self.xs2, self.ys1, self.zs1),),
                ((self.xs2, self.ys2, self.zs1),),
            ),
            technique=FREE,
        )
        self.part.setElementType(
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
                self.part.cells.findAt(
                    ((self.xs1, self.ys1, self.zs1),),
                    ((self.xs1, self.ys2, self.zs1),),
                    ((self.xs2, self.ys1, self.zs1),),
                    ((self.xs2, self.ys2, self.zs1),),
                ),
            ),
        )

        ### Refined area meshing ###
        self.part.seedEdgeBySize(
            constraint=FINER,
            edges=self.part.edges.findAt(
                ((self.xs1, self.ys2 - self.dpo_y / 2.0, self.zs1),),
                ((self.dpo_x, self.ys2 - self.dpo_y / 2.0, self.zs1),),
                ((self.xs1, self.ys2 - self.dpo_y / 2.0, self.zs2),),
                ((self.dpo_x, self.ys2 - self.dpo_y / 2.0, self.zs2),),
            ),
            size=SubstrateSizeY,
        )

        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.dpo_x / 2.0, self.ys2 - self.dpo_y, self.zs1),),
                ((self.dpo_x / 2.0, self.ys2 - self.dpo_y, self.zs2),),
                ((self.dpo_x / 2.0, self.ys2, self.zs2),),
                ((self.dpo_x / 2.0, self.ys2, self.zs1),),
            ),
            size=SubstrateSizeX,
        )

        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.xs1, self.ys2 - self.dpo_y, (self.zs2 + self.zs1) / 2.0),),
                ((self.xs1, self.ys2, (self.zs2 + self.zs1) / 2.0),),
                ((self.dpo_x, self.ys2, (self.zs2 + self.zs1) / 2.0),),
            ),
            size=SubstrateSizeZ,
        )

        ### Meshing the rest ###
        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.xs1, self.ys1, (self.zs2 + self.zs1) / 2.0),),
                ((self.dpo_x, self.ys1, (self.zs2 + self.zs1) / 2.0),),
                ((self.xs2, self.ys2, (self.zs2 + self.zs1) / 2.0),),
                ((self.xs2, self.ys2 - self.dpo_y, (self.zs2 + self.zs1) / 2.0),),
                ((self.dpo_x / 2.0, self.ys1, self.zs1),),
                ((self.dpo_x / 2.0, self.ys1, self.zs2),),
            ),
            size=CoarseMeshSize1,
        )

        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.xs2, self.ys2 - self.dpo_y / 2.0, self.zs1),),
                ((self.xs2, self.ys2 - self.dpo_y / 2.0, self.zs2),),
            ),
            size=self.dpo_y / 2.0,
        )

        self.part.seedEdgeByBias(
            biasMethod=SINGLE,
            constraint=FINER,
            end1Edges=self.part.edges.findAt(
                ((self.dpo_x, (self.ys2 - self.dpo_y) / 2.0, self.zs2),),
            ),
            end2Edges=self.part.edges.findAt(
                ((self.dpo_x, (self.ys2 - self.dpo_y) / 2.0, self.zs1),),
            ),
            maxSize=CoarseMeshSize1,
            minSize=CoarseMeshSize0,
        )

        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.dpo_x, self.ys2 - self.dpo_y, (self.zs2 + self.zs1) / 2.0),),
            ),
            size=CoarseMeshSize0,
        )

        self.part.seedEdgeBySize(
            constraint=FIXED,
            deviationFactor=0.1,
            edges=self.part.edges.findAt(
                ((self.xs2, self.ys1, (self.zs2 + self.zs1) / 2.0),),
            ),
            size=CoarseMeshSize2,
        )

        self.part.generateMesh()
        return self.part


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

    #### ------------------------------ ####
    #               Meshing
    #### ------------------------------ ####

    # Meshing of substrate
    SubstratePart.setMeshControls(
        elemShape=TET,  # TET
        regions=SubstratePart.cells.findAt(
            ((xs1, ys1, zs1),),
            ((xs1, ys2, zs1),),
            ((xs2, ys1, zs1),),
            ((xs2, ys2, zs1),),
        ),
        technique=FREE,  # FREE
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

    return SubstratePart
