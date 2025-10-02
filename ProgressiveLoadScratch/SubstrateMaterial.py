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
import Constants as C


class SubstrateMaterialAssignment:
    """A class to handle the creation and assignment of material properties to the substrate part in an Abaqus model."""

    def __init__(
        self,
        ScratchModel,
        SubstratePart,
        rho,
        youngs_modulus,
        poisson_ratio,
    ):
        """
        Initializes the SubstrateMaterialAssignment class with the given parameters.

        Args:
            ScratchModel: The Abaqus model object.
            SubstratePart: The Abaqus part object of the substrate.
            rho (float): Density of the substrate material.
            youngs_modulus (float): Young's modulus of the substrate material.
            poisson_ratio (float): Poisson's ratio of the substrate material.
        """

        self.ScratchModel = ScratchModel
        self.SubstratePart = SubstratePart
        self.youngs_modulus = youngs_modulus

        self.MaterialName = "SubstrateMaterial"

        if self.MaterialName in self.ScratchModel.materials.keys():
            del self.ScratchModel.materials[self.MaterialName]

        # Creating the material for the substrate
        self.mat = self.ScratchModel.Material(name=self.MaterialName)
        self.mat.Density(table=((rho,),))
        self.mat.Elastic(table=((youngs_modulus, poisson_ratio),))

    def IsotrpopicHardening(self, yield_strength, n):
        """
        Applies isotropic hardening to the material up to 50% total strain

        Args:
            yield_strength (float): The yield strength of the material.
            n (float): The strain hardening exponent.
        Returns:
            mat: The Abaqus material object with isotropic hardening applied.
        """

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

    def JohnsonCookHardening(self, A, B, n, m=0.0, Tm=0.0, Tt=0.0):
        """
        Applies Johnson-Cook hardening to the material.
        Args:
            A (float): Yield stress at reference conditions.
            B (float): Hardening modulus.
            n (float): Strain hardening exponent.
            m (float): Thermal softening exponent.
            Tm (float): Melting temperature of the material.
            Tt (float): Reference temperature.
        Returns:
            mat: The Abaqus material object with Johnson-Cook hardening applied.
        """

        self.mat.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, m, Tm, Tt),))

        return self.mat

    def JohnsonCookDamage(self, d1, d2, d3, d4=0.0, d5=0.0, Tm=0.0, Tt=0.0, Sr=0.0):
        """
        Applies Johnson-Cook damage initiation to the material.

        Args:
            d1 (float): Damage parameter 1.
            d2 (float): Damage parameter 2.
            d3 (float): Damage parameter 3.
            d4 (float): Damage parameter 4.
            d5 (float): Damage parameter 5.
            Tm (float): Melting temperature of the material.
            Tt (float): Reference temperature.
            Sr (float): Reference strain rate.
        Returns:
            mat: The Abaqus material object with Johnson-Cook damage initiation applied.

        """
        self.mat.JohnsonCookDamageInitiation(table=((d1, d2, d3, d4, d5, Tm, Tt, Sr),))
        # eps_f = d1 + d2*exp(-d3*)
        return self.mat

    def DamageEvolution(self, kc, uts, E, nu):
        """
        Assigns damage evolution to the model.

        Args:
            kc (float): Fracture toughness.
            E (float): Young's modulus.
            nu (float): Poisson's ratio
        """
        E_star = E / (1 - nu**2)  # [MPa]
        fractureEnergy = kc**2 / E_star  # [N/mm] = [MJ/m2]
        delta_f = 2 * fractureEnergy / uts  # [mm]
        self.mat.johnsonCookDamageInitiation.DamageEvolution(
            table=((delta_f,),), type=DISPLACEMENT, softening=LINEAR
        )
        return self.mat

    def SectionAssignment(self):
        """
        Assigns the created material to the substrate part by creating and assigning a section.

        Returns:
            ScratchModel: The Abaqus model object with the material and section assigned.
            SubstratePart: The Abaqus part object of the substrate with the section assigned.
        """
        # Creating section
        self.ScratchModel.HomogeneousSolidSection(
            material=self.MaterialName, name="SubstrateSection", thickness=None
        )

        self.SubstratePart.SectionAssignment(
            offset=0.0,
            offsetField="",
            offsetType=MIDDLE_SURFACE,
            region=self.SubstratePart.sets[C.substrate_set_name],
            sectionName="SubstrateSection",
            thicknessAssignment=FROM_SECTION,
        )
        return self.ScratchModel, self.SubstratePart

    def UpdateFriction(self, mu):
        """
        Updates the interfacial coefficient of friction in the tangential behaviour of the contact properties.

        Args:
            mu (float): Interfacial coefficient of friction.

        Returns:
            ScratchModel: The Abaqus model object with updated friction coefficient.

        """
        self.ScratchModel.interactionProperties[
            "IntProp-1"
        ].tangentialBehavior.setValues(table=((mu,),))
        return self.ScratchModel
