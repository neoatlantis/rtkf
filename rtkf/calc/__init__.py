#!/usr/bin/env python3

import numpy as np
import time

from ..cfturbo_reader import CFTurboReader
from ..prop import Properties

from .working_point import WorkingPoint
from .volute import VoluteCalculation
from .nozzle_row import NozzleRowCalculation
from .rotor import RotorCalculation
from .diffuser import DiffuserCalculation
from .estimate_massflow import estimateMassflow
from ._report import CalculationReport
from .briefing import briefing      # generate brief report of results
from .multivar_newton import newton_multi
from .choked import ChokedException
from .supercritical_solver import SupercriticalSolver

from logging import *



class CalculationVariation:
    """Since RTKF needs to calculate a series of points to form a field of
    performances, there will be a need to deviate conditions at starting point
    from design points. This variation will be represented by this class, which
    shall be used as input when beginning a calculation.

      The variation is currently done on mass flow and/or total pressure at
    inlet.
    """

    ratioPressure = 1.0
    ratioMassFlow = 1.0
    ratioRotationSpeed = 1.0

    def p(self, ratio):
        self.ratioPressure = ratio
        return self

    def m(self, ratio):
        self.ratioMassFlow = ratio
        return self

    def N(self, ratio):
        self.ratioRotationSpeed = ratio
        return self

    def __str__(self):
        return "VARIATION: pvar=%f mvar=%f Nvar=%f" % (
            self.ratioPressure, self.ratioMassFlow, self.ratioRotationSpeed
        )



class Calculation(CalculationReport):

    JOBS = [
        (VoluteCalculation, "volute"),
        (NozzleRowCalculation, "nozzle_row"),
        (RotorCalculation, "rotor"),
        (DiffuserCalculation, "diffuser"),
    ]

    def __init__(self, cfturbo):
        assert isinstance(cfturbo, CFTurboReader)
        self.config = cfturbo # the project configuration as input

        CalculationReport.__init__(self)

    @property
    def designPoint(self):
        return self.config.designPoint

    def queryProperties(self, *args):
        # shortcut for query properties, as for all components medium is same
        return Properties(self.workingPoint.medium).query(*args)

    def __callJobs(self, variation, autostart=True, p_after={}):
        assert isinstance(variation, CalculationVariation)
        self.workingPoint = WorkingPoint(self, variation)
        for job, name in self.JOBS:
            x = job(self)
            try:
                if autostart: 
                    if not hasattr(x, "SKIP"):
                        x.solve()
                    else:
                        print("SKIPPED CALCULAETION OF COMPONENT %s!!" % name)
                        x.skip()
            except ChokedException as chokedException:
                if chokedException.where in p_after:
                    x.setChokedAfterPressure(p_after[chokedException.where])
                else:
                    setattr(self, name, x)
                    raise chokedException
            setattr(self, name, x)
        return self

    def calculate(self, variation, p_after={}):
        return self.__callJobs(variation, autostart=True, p_after=p_after)

    def estimateMassflow(self, variation):
        self.__callJobs(variation, autostart=False)
        return estimateMassflow(self)

    def brief(self):
        """Generate a brief report of calculation results. Appendable
        to text file."""
        return briefing(self)



    """def recalculate(self, p_after={}):
        assert self.workingPoint
        isDownstream = False
        for job, name in self.JOBS:
            if not isDownstream:
                if name in p_after:
                    comp = getattr(self, name)
                    comp.setChokedAfterPressure(p_after[name])
                    isDownstream = True
            else:
                x = job(self)
                try:
                    x.solve()
                except ChokedException as chokedException:
                    if chokedException.where in p_after:
                        x.setChokedAfterPressure(p_after[chokedException.where])
                    else:
                        setattr(self, name, x)
                        raise chokedException
                setattr(self, name, x) # recalculate downstream components
        return self"""







class PerformanceAnalyzer:

    RELATIVE_ERROR = 1e-3

    def __init__(self, cfturbo):
        self.calculation = Calculation(cfturbo)
        self.designPoint = self.calculation.designPoint


    def assignedMassflow(self, variation):
        """
        The m-variation in `variation` is assigned massflow. If this massflow
        satisfy a condition given by N- and p-variation, then it is accepted.
        Otherwise, m-variation will be reduced to the maximal acceptable value
        to provide a result that does not get choked.
        """

        m_dot_designed = self.calculation.designPoint.m_dot

        m_dot = variation.ratioMassFlow * m_dot_designed
        m_dot_seeked = m_dot
        m_dot_max = m_dot
        m_dot_min = 0

        print("*" * 80, "assigned massflow %f begin" % m_dot)
        m_dot_history = []

        i = 0
        while True:
            i += 1
            if i > 20:
                error("\n".join(m_dot_history))
                raise Exception("mass flow took too long")

            v = CalculationVariation()
            v.p(variation.ratioPressure)
            v.N(variation.ratioRotationSpeed)
            v.m(m_dot / m_dot_designed)

            try:
                result = self.calculation.clear().calculate(v)
                choked = False
            except ChokedException as chokedException:
                choked = True
                m_dot_choke = chokedException.m_dot

            m_dot_history.append("massflow seeked: %f, input: %f, min: %f, max: %f -> %s" % (
                m_dot_seeked, m_dot, m_dot_min, m_dot_max,
                "choked (m_dot_cr=%f)" % m_dot_choke if choked else ""
            ))


            if choked:
                m_dot_max = min(m_dot_choke, m_dot)
                if m_dot_min > 0:
                    m_dot = (m_dot_max + m_dot_min)/2
                else:
                    m_dot = 0.9 * m_dot_max
            else:
                if abs(m_dot_max / m_dot - 1) < self.RELATIVE_ERROR:
                    return {
                        "result": self.calculation,
                        "m_dot": m_dot,
                        "choked": m_dot < m_dot_seeked,
                    }
                else:
                    m_dot_min = m_dot
                    m_dot = (m_dot_min + m_dot_max)/2
    

    def assignedDischargePressure(self, variation):
        """
        Only the N-variation in `variation` will be respected. Massflow
        and pressure will be adjusted, to find a inlet total pressure that
        generates outlet static pressure as specified in CFTurbo Project data.
        """
        p_total_inlet_designed = \
            float(self.calculation.designPoint.p_total_inlet)
        p_static_outlet_designed = \
            float(self.calculation.designPoint.p_static_outlet)
        m_dot_designed = self.calculation.designPoint.m_dot
        
        m_dot_old = 0
        m_dot_in = self.calculation.estimateMassflow(variation)
        p_old = p_total_inlet_designed * variation.ratioPressure
        #     == self.workingPoint.p_total_inlet

        print("Initial guess of m_dot", m_dot_in)

        # searches: (achievable mass flow, outlet static pressure)
        searches = []
        searches.append((0, p_total_inlet_designed*variation.ratioPressure))

        debug_history = []

        i = 0
        while True:
            i += 1
            if i > 20:
                print("\n".join(debug_history))
                raise Exception("assigned discharge pressure took too long.")

            v = CalculationVariation()
            v.p(variation.ratioPressure).N(variation.ratioRotationSpeed)
            v.m(m_dot_in / m_dot_designed)

            print("*" * 90, "assigned discharge pressure begin")
            print("!!!!! assigned mass flow", m_dot_in)
            print("!!!!! assigned input pressure", p_old)

            analysis = self.assignedMassflow(v)

            m_dot = analysis["m_dot"]
            choked = analysis["choked"]
            p_dis = float(analysis["result"].diffuser.p2)
            print("!!!!! achievable mass flow: ", m_dot)

            searches.append((m_dot, p_dis))
            if len(searches) > 20:
                print("\n".join(debug_history))
                exit()

            debug_history.append(str(analysis) + "  p2_dis=%f" % p_dis)

            print("p_dis", p_dis)
            print("*" * 90, "end", i)
            for each in searches: print(each)

            if choked:
                print("confirm choked", p_static_outlet_designed)

                if p_dis > p_static_outlet_designed:
                    warn("Enter supercritical solution routine.")
                    print("\n" * 20)
                    time.sleep(1)


                    supercriticalSolver = SupercriticalSolver(
                        cfturbo=self.calculation.config,
                        variation=v,
                        Calculation=Calculation
                    )
                    p_dis, calc = supercriticalSolver.solve(
                        p_dis_target=p_static_outlet_designed)
                    return {
                        "m_dot": m_dot,
                        "p_dis": p_dis,
                        "result": calc,
                    }

            print("!!!!! given p", p_old)
            print("!!!!! given m_dot", m_dot, " - choked" if choked else "")
            print("!!!!! current discharge pressure", p_dis)
            print("!!!!! required discharge pressure", p_static_outlet_designed)

            relerr = abs(p_dis / p_static_outlet_designed - 1)
            print("Relative error: %E" % relerr)


            if relerr <= self.RELATIVE_ERROR:
                return {
                    "m_dot": m_dot,
                    "p_dis": p_dis,
                    "result": self.calculation,
                }

            x1, y1 = searches[-2]
            x2, y2 = searches[-1]
            k = (y2-y1)/(x2-x1)
            b = y1 - k * x1

            m_dot_in = (p_static_outlet_designed - b)/k
            print("!!!!! next m_dot for search:", m_dot_in)
