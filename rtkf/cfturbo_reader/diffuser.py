#!/usr/bin/env python3

from ._types import CFTurboObject

class Diffuser(CFTurboObject):

    def __init__(self, xmlNode, geometry):
        CFTurboObject.__init__(self, xmlNode)

        diffuserGeometry = geometry.Component["Pipe out"]
        self.geometry = diffuserGeometry