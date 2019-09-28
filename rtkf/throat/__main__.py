#!/usr/bin/env python3

from . import ThroatFinder

import os
import sys
import argparse

parser = argparse.ArgumentParser(
    description="Find the throat parameters from a given geometry file."
)

parser.add_argument(
    "input",
    metavar="CFTURBO_GEOMETRY_FILE",
    action="store",
    help="The CFTurbo geometry export file, *.geo-xml."
)

parser.add_argument(
    "component",
    help="Name of the component with blades."
)

parser.add_argument(
    "N",
    type=int,
    help="Number of blades."
)


args = parser.parse_args()

x = ThroatFinder(args.input)


o, r_th, b_th = x(args.component, N=args.N)
print("o=%f\nr_th=%f\nb_th=%f\n" % (o, r_th, b_th))
