#!/usr/bin/env python3

import sys
import os
import argparse
import logging

from .cfturbo_reader import CFTurboReader
from .calc import Calculation, CalculationVariation, PerformanceAnalyzer
from .report import serve


# Boot up and self check
assert sys.version_info > (3, 5)


# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(levelname)s: %(message)s"
)

# Settings of argparser

parser = argparse.ArgumentParser(
    prog="python3 -m rtkf",
    description="""Calculate the performance of a given radial inflow turbine.
    """
)

parser.add_argument(
    "input",
    metavar="INPUT_FILE",
    action="store",
    help="The input file, *.cft, *.geo-xml or *.yaml"
)

parser.add_argument(
    "action",
    choices=["point", "massflow", "pressure"]
)

parser.add_argument(
    "--output", "-o",
    action="store",
    help="Write calculation result in given *.txt file."
)


parser.add_argument("--pvar", type=float, default=1.0)
parser.add_argument("--nvar", type=float, default=1.0)
parser.add_argument("--mvar", type=float, default=1.0)

args = parser.parse_args()

# Load project file and read in

cfturbo = CFTurboReader(args.input)

calculator = Calculation(cfturbo)
analyzer = PerformanceAnalyzer(cfturbo)

serve(
    CalculationVariation,
    calculator,
    analyzer,
    action=args.action.upper(),
    output=args.output,
    pVar=args.pvar,
    NVar=args.nvar,
    mVar=args.mvar
)