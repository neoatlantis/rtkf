#!/usr/bin/env python3

from ..prop import Properties
from ..prop.units import *
from math import pi

from ._report import CalculationReport


class WorkingPoint(CalculationReport):

    """Calculations for the volute part.

    Input
    :calc - An instance of calc.Calculation, the main calculation object.
    """

    def __init__(self, calc, variation):
        self.parent = calc
        self.designPoint = calc.config.designPoint

        CalculationReport.__init__(self)

        self.variation = variation
        self.medium = self.designPoint.medium

        self.T_total_inlet = self.designPoint.T_total_inlet
        self.p_total_inlet = \
            self.designPoint.p_total_inlet * variation.ratioPressure
        self.m_dot = self.designPoint.m_dot * variation.ratioMassFlow
        self.n = self.designPoint.nRot * variation.ratioRotationSpeed
        self.omega = self.n * 2 * pi 

        self.set(
            "properties", 
            Properties(self.medium).query(
                Temperature(self.T_total_inlet).units(),
                Pressure(self.p_total_inlet).units()),
            skip=True
        )



