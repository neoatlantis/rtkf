#!/usr/bin/env python3

import CoolProp.CoolProp as CP

class Properties:

    def __init__(self, fluid):
        self.fluid = fluid

    def query(self, what, **queries):
        keys = queries.keys()
        assert len(keys) == 2
        return CP.PropsSI(
            what,
            keys[0], queries[keys[0]],
            keys[1], queries[keys[1]],
            self.fluid
        )