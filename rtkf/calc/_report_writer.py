#!/usr/bin/env python3

import enum
from math import pi
from ..prop.units import PhysicalQuantity

class ReportWriter:

    class modes(enum.Enum):
        PLAINTEXT = 1
        XHTML = 2

    class formatters(enum.Enum):
        DEFAULT = lambda i: self.__conv(i)
        ANGLE = lambda i: "%.2f\u00B0" % (i * 180 / pi)
        PERCENT = lambda i: "%.2f%%" % (i * 100.0)


    def __conv(self, i):
        if self.mode == self.modes.PLAINTEXT:
            if type(i) == float: return "%E" % i
            return str(i)
        if self.mode == self.modes.XHTML:
            if type(i) == float:
                s = "%E" % i
                a, b = s.split("E")
                return "%s &times;10<sup>%d</sup>" % (a, int(b))
            return str(i)

    def __init__(self, mode):
        assert type(mode) == self.modes
        self.mode = mode

    def comment(self, comment):
        if self.mode == self.modes.PLAINTEXT:
            return comment
        if self.mode == self.modes.XHTML:
            return "<tr><td colspan=3>%s</td></tr>" % comment.strip()

    def entry(self, key, value, comment="", formatter=None):
        if not formatter: formatter = self.__conv
        if self.mode == self.modes.PLAINTEXT:
            return "%25s = %s" % (key, self.__conv(value)) +\
                ("  #%s" % comment if comment else "")
        if self.mode == self.modes.XHTML:
            if isinstance(value, PhysicalQuantity):
                value = "<physical type=\"%s\">%s %s</physical>" % (
                    value.__class__.__name__.lower(),
                    formatter(float(value)),
                    value.defaultUnit
                )
            return "<tr>%s%s%s</tr>" % (
                "<td class=\"name\">%s</td>" % key,
                "<td class=\"value\">%s</td>" % formatter(value),
                "<td class=\"comment\">%s</td>" % comment
            )

    def subreport(self, name, report):
        if self.mode == self.modes.PLAINTEXT:
            suboutput = str(report).split("\n")
            return "\n".join([
                ("%16s + " % key if i == 0 else " " * 16 + " | ") +\
                suboutput[i]
                for i in range(0, len(suboutput))
            ])
        if self.mode == self.modes.XHTML:
            return str(report)

    def report(self, name, content):
        if self.mode == self.modes.PLAINTEXT:
            ret = ("---- Calculation: %s ----\n" % name) + content
            return ret
        if self.mode == self.modes.XHTML:
            return "<table><caption>%s</caption>%s</table>" % (
                name, content
            )
