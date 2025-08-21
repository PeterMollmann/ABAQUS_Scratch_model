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


class SubstrateMaterialAssignment:
    def __init__(
        self,
        ScratchModel,
        SubstratePart,
        SubstrateSet,
        rho,
        youngs_modulus,
        poisson_ratio,
    ):
        self.ScratchModel = ScratchModel
        self.SubstratePart = SubstratePart
        self.SubstrateSet = SubstrateSet
        self.youngs_modulus = youngs_modulus

        self.MaterialName = "SubstrateMaterial"

        if self.MaterialName in self.ScratchModel.materials.keys():
            del self.ScratchModel.materials[self.MaterialName]

        # Creating the material for the substrate
        self.mat = self.ScratchModel.Material(name=self.MaterialName)
        self.mat.Density(table=((rho,),))
        self.mat.Elastic(table=((youngs_modulus, poisson_ratio),))

    def IsotrpopicHardening(self, yield_strength, n):
        total_strain_data = np.arange(0.0001, 0.5, 0.01)
        plastic_behaviour = []
        yield_strain = yield_strength / self.youngs_modulus
        K = self.youngs_modulus * yield_strain ** (1 - n)
        plastic_behaviour.append((yield_strength, 0))
        for total_strain in total_strain_data:
            if total_strain > yield_strain:
                yield_strength_ = round(K * total_strain**n, 5)
                plastic_strain = round(
                    total_strain - yield_strength_ / self.youngs_modulus, 5
                )
                plastic_behaviour.append((yield_strength_, plastic_strain))

        plastic_behaviour = tuple(plastic_behaviour)
        self.mat.Plastic(table=plastic_behaviour)
        return self.mat

    def JohnsonCookHardening(self, A, B, n, m, Tm, Tt):
        self.mat.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, m, Tm, Tt),))

        return self.mat

    def JohnsonCookDamage(self, d1, d2, d3, d4, d5, Tm, Tt, Sr):
        self.mat.JohnsonCookDamageInitiation(table=((d1, d2, d3, d4, d5, Tm, Tt, Sr),))
        return self.mat

    def SectionAssignment(self):
        # Creating section
        self.ScratchModel.HomogeneousSolidSection(
            material=self.MaterialName, name="SubstrateSection", thickness=None
        )

        self.SubstratePart.SectionAssignment(
            offset=0.0,
            offsetField="",
            offsetType=MIDDLE_SURFACE,
            region=self.SubstratePart.sets[self.SubstrateSet],
            sectionName="SubstrateSection",
            thicknessAssignment=FROM_SECTION,
        )
        return self.ScratchModel, self.SubstratePart
