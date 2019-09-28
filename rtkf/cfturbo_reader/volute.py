#!/usr/bin/env python3

from ._types import CFTurboObject


class Volute(CFTurboObject):

    def __init__(self, xml, geometry):

        CFTurboObject.__init__(self, xml)

        try:
            assert self.PrimaryInlet and not self.PrimaryOutlet
        except:
            print("WARNING: Volute seems not configured as main inlet!")

#        self.design = CFTurboObject(xml.find("SpiralCasingSection"))
#        self.design = CFTurboObject(xml.find("ModelSetup_Spiral"))
#        self.design = CFTurboObject(xml)

        self.geometry = geometry.Component["Volute"]