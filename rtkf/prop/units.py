#!/usr/bin/env python3

"""
Defines a class representing a physical value(a value with units).

An abstraction class is named `PhysicalQuantity`, classes like `Pressure`,
`Temperature` are inherited from this.

Usage:
    Pressure(101325).units('Pa')       creates an instance for `101325 Pa`
    Pressure(101325).units('Pa').to('MPa')
                                       returns `0.101325`
"""
import numpy as np

def isNumber(i):
    return type(i) in [
        int,
        float,
        np.float64,
    ]

class PhysicalQuantitiesHolder:
    pass

class PhysicalQuantity:
    def __init__(self, value):
        self.__v = value

    def denote(self, unitName=None):
        """Returns a human-readable unit symbol representation."""
        if unitName == None:
            unitName = self.defaultUnit
        if hasattr(self, "displayTable") and unitName in self.displayTable:
            return self.displayTable[unitName].encode("utf-8")
        else:
            return unitName
    
    def units(self, unitName=None):
        if unitName == None:
            unitName = self.defaultUnit
        else:
            if not unitName in self.inputTable:
                raise ValueError(
                    "Unsupported input unit `%s` for type `%s`" %
                    (unitName, self.__class__.__name__)
                )
        self.__unit = unitName
        return self

    def to(self, unitName):
        if not unitName in self.outputTable:
            raise ValueError(
                "Unsupported output unit `%s` for type `%s`" %
                (unitName, self.__class__.__name__)
            )
        return float(self.outputTable[unitName](
            self.inputTable[self.__unit](self.__v)
        ))

    def getInternalValue(self):
        return self.__v

    def __float__(self):
        return self.to(self.defaultUnit)

    def __repr__(self):
        return "<%s> %E %s" % (
            self.__class__.__name__,
            self.__v,
            self.__unit
        )

    def __eq__(self, you):
        if self.__class__.__name__ != you.__class__.__name__: return False
        return self.getInternalValue() == you.getInternalValue()

    def __add__(self, you):
        if self.__class__.__name__ != you.__class__.__name__:
            if type(you) != float:
                raise Exception(
                    "Cannot add 2 values of different physical quantities.")
        newValue = float(self) + float(you)
        if type(you) == float: return newValue
        return self.__class__(newValue).units()

    def __sub__(self, you):
        if self.__class__.__name__ != you.__class__.__name__:
            if type(you) != float:
                raise Exception(
                    "Cannot substract 2 values of different physical quantities.")
        newValue = float(self) - float(you)
        if type(you) == float: return newValue
        return self.__class__(newValue).units()

    def __truediv__(self, you):
        # e.g.: self / another physical quantity
        if isNumber(you):
            you = float(you)
            return self.__class__(float(self) / you).units()
        elif self.__class__.__name__ == you.__class__.__name__:
            return float(self) / float(you)
        else:
            raise Exception("Cannot calculate division.")

    def __div__(self, you):
        return self.__truediv__(you)
        
    def __mul__(self, you):
        # e.g.: self * 1.2
        if isNumber(you):
            you = float(you)
            return self.__class__(float(self) * you).units()
        else:
            print(type(you))
            raise Exception("Cannot calculate multiplication.")
        
    def __rmul__(self, you):
        # e.g.: 1.2 * self
        if isNumber(you):
            you = float(you)
            return self.__class__(float(self) * you).units()
        else:
            print(type(you))
            raise Exception("Cannot calculate multiplication.")



class Speed(PhysicalQuantity):
    defaultUnit = "m/s"
    inputTable = {
        "m/s": lambda i: i,
        "km/h": lambda i: i / 3.6,
    }
    outputTable = {
        "m/s": lambda i: i,
        "km/h": lambda i: i * 3.6,
    }

class Pressure(PhysicalQuantity):
    defaultUnit = "Pa"
    inputTable = {
        "Pa":  lambda i: i,
        "MPa": lambda i: i * 1.0e6,
        "kPa": lambda i: i * 1.0e3,
        "bar": lambda i: i * 1.0e5,
        "atm": lambda i: i * 101325.0,
    }
    outputTable = {
        "Pa": lambda v: v,
        "MPa": lambda v: v / 1.0e6,
        "kPa": lambda v: v / 1.0e3,
        "bar": lambda v: v / 1.0e5,
        "atm": lambda v: v / 101325.0
    }

class Temperature(PhysicalQuantity):
    defaultUnit = "K"
    inputTable = {
        "K": lambda i: i,
        "C": lambda i: i + 273.15,
    }
    outputTable = {
        "K": lambda v: v,
        "C": lambda v: v - 273.15,
    }
    displayTable = {
        "C": u"\u2103",
    }

class Density(PhysicalQuantity):
    defaultUnit = "kg/m^3"
    inputTable = {
        "kg/m^3": lambda i: i,
        "g/cm^3": lambda i: i * 1000.0,
    }
    outputTable = {
        "kg/m^3": lambda i: i,
        "g/cm^3": lambda i: i / 1000.0,
    }
    displayTable = {
        "kg/m^3": u"kg/m\u00B3",
        "g/cm^3": u"g/cm\u00B3",
    }

class Entropy(PhysicalQuantity):
    defaultUnit = "J/kg/K"
    inputTable = {
        "J/kg/K": lambda i: i,
        "kJ/kg/K": lambda i: i * 1000.0,
    }
    outputTable = {
        "J/kg/K": lambda v: v,
        "kJ/kg/K": lambda v: v / 1000.0,
    }

class Enthalpy(PhysicalQuantity):
    defaultUnit = "J/kg"
    inputTable = {
        "J/kg": lambda i: i,
        "kJ/kg": lambda i: i * 1000.0,
    }
    outputTable = {
        "J/kg": lambda v: v,
        "kJ/kg": lambda v: v / 1000.0,
    }

class SpecificGasConstant(PhysicalQuantity):
    defaultUnit = "J/kg/K"
    inputTable = {
        "J/kg/K": lambda i: i,
    }
    outputTable = {
        "J/kg/K": lambda v: v,
    }

class HeatCapacityRatio(PhysicalQuantity):
    defaultUnit = "kappa"
    inputTable = { "kappa": lambda i: i, "K": lambda i: i * 1.0 / (i - 1.0) }
    outputTable = { "kappa": lambda o: o, "K": lambda o: o * 1.0 / (o - 1.0) }


class ThermalConductivity(PhysicalQuantity):
    defaultUnit = "W/m/K"
    inputTable = { "W/m/K": lambda i: i }
    outputTable = { "W/m/K": lambda o: o }

class DynamicViscosity(PhysicalQuantity):
    defaultUnit = "Pa.s"
    inputTable = { "Pa.s": lambda i: i }
    outputTable = { "Pa.s": lambda o: o }

class KinematicViscosity(PhysicalQuantity):
    defaultUnit = "m^2/s"
    inputTable = { "m^2/s": lambda i: i }
    outputTable = { "m^2/s": lambda o: o }
