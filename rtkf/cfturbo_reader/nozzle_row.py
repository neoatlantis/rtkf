#!/usr/bin/env python3

from ._types import CFTurboObject

class NozzleRow(CFTurboObject):


    def __init__(self, xmlNode, geometry):
        CFTurboObject.__init__(self, xmlNode)

        statorGeometry = geometry.Component["Stator"]
        self.geometry = statorGeometry

