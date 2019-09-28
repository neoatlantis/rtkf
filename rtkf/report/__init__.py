#!/usr/bin/env python3

from bottle import run, get, static_file
import sys
import os

rootpath = os.path.join(os.path.dirname(sys.argv[0]), "report", "static")


def finalize(report):
    return """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
<head>

<link rel="stylesheet" href="/static/style.css" />

</head>
<body>""" + report + """</body>
</html>
"""


def serve(
    variator, calculator, analyzer,
    action="POINT",
    output=None,
    pVar=1.0, NVar=1.0, mVar=1.0
):

    v = variator().m(mVar).p(pVar).N(NVar)
    if "MASSFLOW" == action:
        print("assigned mass flow")
        result = analyzer.assignedMassflow(v)["result"]
    if "PRESSURE" == action:
        print("assigned outlet static pressure(from project file)")
        result = analyzer.assignedDischargePressure(v)["result"]
    if "POINT" == action:
        result = calculator.calculate(v)

    if output:
        if not output.endswith(".csv"):
            raise Exception("Briefing report must end with .csv")

        brief = result.brief()
        briefKeys = sorted(brief.keys())
        if not os.path.isfile(output):
            open(output, "w+").write(
                ",".join(["\"%s\"" % e for e in briefKeys]) + "\n"
            )
        open(output, "a").write(
            ",".join(["%f" % brief[e] for e in briefKeys]) + "\n"
        )

        exit()

    report = finalize(str(result))
    @get("/")
    def index():
        return report



    @get("/static/<filepath>")
    def static(filepath):
        global rootpath
        return static_file(filepath, root=rootpath)

    run(host="localhost", port=10967)
