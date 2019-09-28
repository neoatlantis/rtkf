#!/usr/bin/env python3

import xml.etree.ElementTree as ET

from ..cfturbo_reader._types import CFTurboObject

from ..calc.parametric_curve_distance import *
from ..calc.curve import DiscreteCurve, BezierCurve, MiddleOf2Curves
from numpy import *

import os



class ThroatFinder:

    def __init__(self, inputFile):

        """Read in a CFTurbo project file. Assure it contains all the 4
        components required in our calculation: volute, nozzle ring, radial
        impeller and exhaust. Besides, it must contain a design point profile.
        After that, call sub modules to convert each component from XML to
        some more friendly internal format.
        """

        if not inputFile.endswith(".geo-xml"):
            raise Exception("Input file must be .geo-xml file.")

        baseName = os.path.splitext(inputFile)[0]
        geometryFile = baseName + ".geo-xml"

        if not os.path.isfile(geometryFile):
            raise Exception("Requires geometry file: %s" % geometryFile)

        self.geometryTree = ET.parse(geometryFile)
        print("Read: %s" % geometryFile)

        root2 = self.geometryTree.getroot()

        self.nodeGeometry = root2
        if root2.tag != "CFturboGeometryXML":
            raise Exception("Invalid CFTurbo geometry input file.")

        self.geometry = CFTurboObject(
            self.nodeGeometry, checkType=False, defaultType="Object")


    def __call__(self, name=None, N=14):
        if not name:
            print(self.geometry.Component.keys())
            exit()

        MainBladeMeanlines = \
            self.geometry.Component[name].Blade["Main"].Meanlines

        mainCurve3dHub = MainBladeMeanlines[0]
        mainCurve3dMid = MainBladeMeanlines[int(len(MainBladeMeanlines)/2)]
        mainCurve3dShr = MainBladeMeanlines[-1]
        o, r_th, b_th = self.__findThroat(
            mainCurve3dHub, mainCurve3dMid, mainCurve3dShr,
            N=N
        )
        return o, r_th, b_th



    def __findThroat(self, curve3dHub, curve3dMid, curve3dShr, N=14):
        delta = 2 * pi / N

        i1, i2, omid = findThroatFromCurve3d(delta, curve3dMid)

        # TODO substract omid with blade thickness!

        M_th = curve3dMid[i1].M
        r_th = (curve3dMid[i1].r + curve3dMid[i2].r)/2
        t = M_th / curve3dMid[-1].M

        hubRtoM = DiscreteCurve([np.array([e.M, e.r]) for e in curve3dHub])
        hubZtoM = DiscreteCurve([np.array([e.M, e.z]) for e in curve3dHub])
        shrRtoM = DiscreteCurve([np.array([e.M, e.r]) for e in curve3dShr])
        shrZtoM = DiscreteCurve([np.array([e.M, e.z]) for e in curve3dShr])

        hubR, hubZ = hubRtoM(t)[1], hubZtoM(t)[1]
        shrR, shrZ = shrRtoM(t)[1], shrZtoM(t)[1]

        bmid = sqrt((shrR-hubR)**2 + (shrZ-hubZ)**2)
        return omid, r_th, bmid
