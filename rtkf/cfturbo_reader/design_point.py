#!/usr/bin/env python3

from ._types import CFTurboObject
from ..prop.units import *


class DesignPoint(CFTurboObject):

    def __init__(self, xml):

        CFTurboObject.__init__(self, xml, defaultType="String")

        pi_t_s = getattr(self, "PressureRatio-t-s")

        self.medium = list(self.FluidGas.keys())[0]
        self.T_total_inlet = Temperature(self.TtInlet).units("K")
        self.p_static_outlet = Pressure(self.pOutlet).units("Pa")
        self.p_total_inlet = self.p_static_outlet * pi_t_s

        self.m_dot = float(self.mFlow)
        assert self.FlowMode == "Mass flow"