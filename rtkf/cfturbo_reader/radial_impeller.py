#!/usr/bin/env python3

from ._types import CFTurboObject

class RadialImpeller(CFTurboObject):

    def __init__(self, xmlNode, geometry):

        CFTurboObject.__init__(self, xmlNode)

        rotorGeometry = geometry.Component["Rotor"]
        self.geometry = rotorGeometry
