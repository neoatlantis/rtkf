#!/usr/bin/env python3

from numpy import sin, cos, arctan
from ._report_writer import ReportWriter
from ..cfturbo_reader._types import YAMLObject

# Formatter
writer = ReportWriter(ReportWriter.modes.XHTML)



class DataHolder:

    def __init__(self, **kvargs):
        self.__dict__["__kvargs"] = kvargs

    def __getattr__(self, key):
        return self.__dict__["__kvargs"]["key"]



class CalculationReport:

    def __init__(self):
        """This class is used to generate a report on calculation progress.
        It tracks each variable assignment and records the value into the
        final report.
          When inheriting this class, call this parent's init method only when
        the calculation is about to begin."""
        self.clear()
        self.set("formatters", writer.formatters, skip=True)

    def override(self, yamlobj):
        assert isinstance(yamlobj, YAMLObject)
        for key in yamlobj.parameters:
            self.set(key, yamlobj.parameters[key])

    def skip(self):
        print("\n\n #### %s #### CALCULATION SKIPPED! ####\n\n" % 
            self.__class__.__name__)

    def clear(self):
        self.__dict__["__report"] = []
        return self

    def comment(self, comment=""):
        self.__dict__["__report"].append(writer.comment(comment))

    def set(self, key, value, comment="", skip=False, formatter=None):
        if isinstance(value, CalculationReport):
            self.comment(comment)
            self.__setattr__(key, value)
            return
        self.__dict__[key] = value
        if not skip and "__report" in self.__dict__:
            self.__dict__["__report"].append(writer.entry(
                key, value, comment=comment, formatter=formatter))

    def __setattr__(self, key, value):
        self.__dict__[key] = value
        if "__report" in self.__dict__:
            if isinstance(value, CalculationReport):
                report = writer.subreport(key, value)
            else:
                report = writer.entry(key, value)
            self.__dict__["__report"].append(report)

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            raise AttributeError("%r object has no attribute %r" %
                (self.__class__.__name__, key))

    def __str__(self):
        return writer.report(
            self.__class__.__name__,
            "\n".join(self.__dict__["__report"]))

    def conclude(self, **args):
        """Records output streaming properties."""
        if "c" in args and "alpha" in args:
            alpha = args["alpha"]
            c = args["c"]
            c_m = c * sin(alpha)
            c_theta = c * cos(alpha)
        elif "c_m" in args and "c_theta" in args:
            c_m = args["c_m"]
            c_theta = args["c_theta"]
            c = sqrt(c_m**2 + c_theta**2)
            alpha = atan(c_m / c_theta)
        else:
            raise Exception(
                self.__class__.__name__ + 
                ": insufficient conclusion on velocity.")

        if "h" in args:
            h = args["h"]
            h_t = h + 0.5 * c**2
        elif "h_t" in args:
            h = args.h_t - 0.5 * c**2
            h_t = args["h_t"]
        else:
            raise Exception(
                self.__class__.__name__ +
                ": insufficient conclusion on enthalpy.")

        if "s" not in args:
            raise Exception(
                self.__class__.__name__ +
                ": insufficient conclusion on entropy.")

        self.__conclusion = DataHolder(
            c = c,
            c_theta = c_theta,
            c_m = c_m,
            alpha = alpha,
            h = h,
            h_t = h_t,
            s = s
        )

    @property
    def conclusion(self):
        return self.__conclusion 

