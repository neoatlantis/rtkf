#!/usr/bin/env python3

from numpy import *

def briefing(self):
    report = {}

    report["pvar"] = self.workingPoint.variation.ratioPressure
    report["mvar"] = self.workingPoint.variation.ratioMassFlow
    report["nvar"] = self.workingPoint.variation.ratioRotationSpeed

    report["m_dot"] =     self.diffuser.m_dot
    report["p_in_total"] = self.workingPoint.p_total_inlet
    report["T_in_total"] = self.workingPoint.T_total_inlet
    report["p_out"] = self.diffuser.p2
    report["p_out_total"] = self.diffuser.p2_total
    report["n"] = self.workingPoint.n
    report["omega"] = self.workingPoint.omega

    report["p_ratio_total2total"] =\
        float(report["p_in_total"]) / float(report["p_out_total"])

    report["p_ratio_total2static"] =\
        float(report["p_in_total"]) / float(report["p_out"])

    report["b2j_ratio"] = self.rotor.speed_ratio
    report["eta_total2static"] = self.rotor.eta_total2static

    

    return report