#!/usr/bin/env python3

import xml.etree.ElementTree as ET

from ._types import CFTurboObject, YAMLObject

from .design_point import DesignPoint
from .volute import Volute
from .nozzle_row import NozzleRow
from .radial_impeller import RadialImpeller
from .diffuser import Diffuser

import os
import yaml


class CFTurboReader:

    def __init__(self, inputFile):

        """Read in a CFTurbo project file. Assure it contains all the 4
        components required in our calculation: volute, nozzle ring, radial
        impeller and exhaust. Besides, it must contain a design point profile.
        After that, call sub modules to convert each component from XML to
        some more friendly internal format.
        """

        if not (
            inputFile.endswith(".cft") or 
            inputFile.endswith(".geo-xml") or
            inputFile.endswith(".yaml")
        ):
            raise Exception("Input file must be .cft/.geo-xml/.yaml file.")

        self.overriding = None
        if inputFile.endswith(".yaml"):
            self.overriding = yaml.load(open(inputFile, "r").read())
            self.__initComponents()
            return

        baseName = os.path.splitext(inputFile)[0]
        projectFile = baseName + ".cft"
        geometryFile = baseName + ".geo-xml"

        if not os.path.isfile(projectFile):
            raise Exception("Requires project file: %s" % projectFile)

        if not os.path.isfile(geometryFile):
            raise Exception("Requires geometry file: %s" % geometryFile)

        self.projectTree = ET.parse(projectFile)
        print("Read: %s" % projectFile)

        self.geometryTree = ET.parse(geometryFile)
        print("Read: %s" % geometryFile)

        root = self.projectTree.getroot()
        root2 = self.geometryTree.getroot()

        self.nodeProject = root.find("CFturboProject")
        if not self.nodeProject:
            raise Exception("Invalid CFTurbo project input file.")

        self.nodeGeometry = root2
        if root2.tag != "CFturboGeometryXML":
            raise Exception("Invalid CFTurbo geometry input file.")

        # ---- Design point

        self.nodeDesignPoint = self.nodeProject.find("DesignPoint_Turb")
        if not self.nodeDesignPoint:
            raise Exception(
                "Input CFTurbo project doesn't contain turbine designs.")

        namespace = lambda x: "CFturboDesign_%s" % x

        # ---- Volute

        self.nodeVolute = self.nodeProject.find(namespace("Volute"))
        if not self.nodeVolute:
            raise Exception("Project must contain a volute.")

        # ---- Impeller

        self.nodeRadialImpeller = self.nodeProject.find(
            namespace("RadialImpeller"))
        if not self.nodeRadialImpeller:
            raise Exception("Project must contain a radial impeller.")

        stators = self.nodeProject.findall(namespace("Stator"))
        if len(stators) != 2:
            raise Excetpion(
                "Project must contain a nozzle ring and an exhaust.")
        for stator in stators:
            if stator.find("BladeProfiles"):
                self.nodeNozzleRow = stator
            else:
                self.nodeDiffuser = stator
        assert self.nodeNozzleRow != None and self.nodeDiffuser != None

        self.__initComponents()


    def __initComponents(self):

        """Convert XML based data into instances."""

        if not self.overriding:
            self.geometry = CFTurboObject(
                self.nodeGeometry, checkType=False, defaultType="Object")

            self.designPoint = DesignPoint(self.nodeDesignPoint)
            self.volute = Volute(self.nodeVolute, self.geometry)
            self.nozzleRow = NozzleRow(self.nodeNozzleRow, self.geometry) 
            self.rotor = RadialImpeller(self.nodeRadialImpeller, self.geometry)
            self.diffuser = Diffuser(self.nodeDiffuser, self.geometry)

        else:

            print("Overriding geometry using YAML.")

            self.designPoint = YAMLObject(self.overriding["design_point"])
            self.volute =      YAMLObject(self.overriding["volute"])
            self.nozzleRow =   YAMLObject(self.overriding["nozzle_row"])
            self.rotor =       YAMLObject(self.overriding["rotor"])
            self.diffuser =    YAMLObject(self.overriding["diffuser"])
