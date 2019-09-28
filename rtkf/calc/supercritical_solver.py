#!/usr/bin/env python3


import numpy as np
from logging import *
import time

from .choked import ChokedException



class SupercriticalSolver:

    RELATIVE_ERROR = 1e-3

    def __init__(self, cfturbo, variation, Calculation):
        self.cfturbo = cfturbo
        self.variation = variation
        self.Calculation = Calculation

    def newCalculationInstance(self, p_after):
        c = self.Calculation(self.cfturbo)
        return c.calculate(variation=self.variation, p_after=p_after)

    def solve(self, p_dis_target):
        """Solve for target discharge static pressure p_dis, and returns the
        result."""

        p_afters = {}

        while True:
            choked1 = False
            try:
                calculation1 = self.newCalculationInstance(p_afters)
                break # if nothing happens here, nothing chokes, just exit.
            except ChokedException as chokedException1:
                choked1 = True
                where1 = chokedException1.where
                assert where1 in ["rotor", "nozzle_row"]
                p_range1 = [
                    float(chokedException1.p_min),
                    float(chokedException1.p_max)
                ]
                if where1 == "rotor":
                    # is the last component chokable
                    if p_range1[0] < p_dis_target:
                        p_range1[0] = p_dis_target
                        p_afters[where1] = sum(p_range1) / 2.0
                    else:
                        # p_min > p_dis_target
                        p_afters[where1] = p_range1[0]
                        print("Expansion exit > required discharge pressure.")
                        exit()
                else:
                    p_afters[where1] = sum(p_range1) / 2.0
                assert p_range1[0] <= p_range1[1]

            # Below are procedures for the component detected above.

            secondChokeEverSensed = False

            while choked1:
                choked2 = False
                try:
                    calculation2 = self.newCalculationInstance(p_afters)
                except ChokedException as chokedException2:
                    choked2 = True
                    where2 = chokedException2.where
                    assert where2 != where1


                if not choked2:
                    p_dis = float(calculation2.diffuser.p2)
                    if abs(p_dis / float(p_dis_target) - 1) < self.RELATIVE_ERROR:
                        break
                    assert p_range1[0] <= p_range1[1]
                    print(p_dis, p_dis_target, p_range1)
                    time.sleep(1)
                    if p_dis < p_dis_target:
                        p_range1[0] = p_afters[where1]
                    else:
                        p_range1[1] = p_afters[where1]
                    assert p_range1[0] <= p_range1[1]

                else:
                    secondChokeEverSensed = True

                    # Otherwise: downstream is thought to be chokable. we will
                    # adjust p_afters[where1] (the back pressure for component 1
                    # and the input pressure for component 2) to such a value,
                    # that the second component is just at the edge of not being
                    # choked(but still choked).

                    warn("When analyzing %s, downstream %s also choked." % (
                        where1, where2
                    ))

                    p_range1[0] = p_afters[where1]
                    print(">>>>>>>>>>", p_range1)
                    time.sleep(1)

                    
                if abs(p_range1[1]/p_range1[0] - 1) < self.RELATIVE_ERROR:
                    break
                p_afters[where1] = sum(p_range1) / 2.0

        # ---- end of loop

        try:
            calc = calculation2
        except:
            calc = calculation1

        p_dis = float(calc.diffuser.p2)

        print(p_afters)
        print(p_dis)


        return p_dis, calc



