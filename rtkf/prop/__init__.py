#!/usr/bin/env python

from CoolProp.CoolProp import PropsSI, FluidsList
from .units import PhysicalQuantity, PhysicalQuantitiesHolder,\
                   Pressure, Temperature, Enthalpy, Entropy, Density,\
                   SpecificGasConstant, HeatCapacityRatio,\
                   ThermalConductivity, DynamicViscosity,\
                   KinematicViscosity, Speed


IDEALGAS = "****IDEALGAS****"

_USE_REFPROP = False
_REFPROP_PATH = ""



class Properties(PhysicalQuantitiesHolder):
    
    def __init__(self, medium):
        self.__medium = medium
        self.__useREFPROP = False

    def __checkQueryInputRealgas(self, v):
        """Check an input variable for a new query. Assure the input is a
        `PhysicalQuantity` instance, assure it can be used to construct an input, 
        determine its representation, and convert it to the suitable unit."""
        if not isinstance(v, PhysicalQuantity):
            raise ValueError("Input must be an instance of `PhysicalQuantity`")

        table = {
            "P": (Pressure, "Pa"),
            "H": (Enthalpy, "J/kg"),
            "T": (Temperature, "K"),
            "S": (Entropy, "J/kg/K"),
        }
        for rep in table:
            if isinstance(v, table[rep][0]):
                return (rep, v.to(table[rep][1]))
        raise ValueError("Given input cannot be used for a query.")

    def useREFPROP(self, path=""):
        """Tell the program to use REFPROP for calculation."""
        self.__useREFPROP = True

    def query(self, var1, var2, *argv):
        mediumName = self.__medium
        if mediumName == IDEALGAS:
            return self.__loadIdealgas(var1, var2, *argv)

        if self.__useREFPROP: mediumName = "REFPROP::" + mediumName

        var1name, var1value = self.__checkQueryInputRealgas(var1)
        var2name, var2value = self.__checkQueryInputRealgas(var2)

        doQuery = lambda var: PropsSI(\
            var, var1name, var1value, var2name, var2value, mediumName)

        self.a = Speed(doQuery("A")).units("m/s")
        self.h = Enthalpy(doQuery('H')).units("J/kg")
        self.s = Entropy(doQuery('S')).units("J/kg/K")
        self.t = Temperature(doQuery('T')).units("K")
        self.T = self.t
        self.p = Pressure(doQuery('P')).units("Pa")
        self.rho = Density(doQuery('D')).units("kg/m^3")
        
        self.Q = doQuery("Q")
        if not (0 <= self.Q <= 1): self.Q = 1
        

        self.mu = DynamicViscosity(doQuery("V")).units("Pa.s")
        self.nu = KinematicViscosity(float(self.mu) / float(self.rho)).units()

        return self

def setREFPROP(use, path=""):
    global _USE_REFPROP, _REFPROP_PATH
    _USE_REFPROP = False
    return # Disable usage of REFPROP
    _USE_REFPROP = bool(use)
    _REFPROP_PATH = path

def usingREFPROP():
    return False # Disable usage of REFPROP
    global _USE_REFPROP, _REFPROP_PATH
    return bool(_USE_REFPROP)

def getPropertiesQuerier(medium):
    global _USE_REFPROP, _REFPROP_PATH
    r = Properties(medium)
    if _USE_REFPROP:
        r.useREFPROP(_REFPROP_PATH)
    return r

def listFluids():
    return FluidsList()
