#!/usr/bin/env python3

def convertXMLNode(node, defaultType="String"):
    attrib = node.attrib
    t = attrib["Type"] if "Type" in attrib else defaultType

    # If type is another object, do that recursively 
    if t == "Object": return CFTurboObject(node)
    if t == "Vector2": return CFTurboObject(node)
    if t == "Vector3": return CFTurboObject(node)
    # If type is simple, return that
    simpleTypes = {
        "Boolean": lambda: (node.text == "True" or node.text == "1"),
        "Integer": lambda: int(node.text),
        "Float": lambda: float(node.text),
        "String": lambda: node.text,
        "Enum": lambda: CFTurboEnum(name=node.tag, choice=node.text),
    }
    if t in simpleTypes: return simpleTypes[t]()
    # Deal with array
    if t.startswith("Array"): return CFTurboArray(node)

    print("Warning: unknown type - %s" % t)
    return None


class CFTurboEnum:

    def __init__(self, name, choice):
        self.__dict__["choice"] = choice 

    def __setattr__(self, *args):
        raise Exception("Enum type not allowed to change.")

    def __repr__(self):
        return self.__dict__["choice"]

    def __eq__(self, v):
        return self.__dict__["choice"] == v


from scipy import interpolate

class CFTurboArray(list):

    def __init__(self, node, defaultType="Object"):
        childNodes = list(node)

        probeNode = childNodes[0]
        if "Index" not in probeNode.attrib:
            list.__init__(
                self,
                [
                    convertXMLNode(e, defaultType=defaultType)
                    for e in childNodes
                ])
            return

#        startsWithArray = node.attrib["Type"].startswith("Array")
        maxIndex = 0
        for childNode in childNodes:
            if "Index" not in childNode.attrib: 
                print(childNode)
                raise Exception("Not an array.")
            maxIndex = max(maxIndex, int(childNode.attrib["Index"]))

        list.__init__(self, [None] * (maxIndex+1))
        for childNode in childNodes:
            index = int(childNode.attrib["Index"])
            self[index] = convertXMLNode(childNode, defaultType=defaultType)

        # Interpolate None values
        # CFTurbo sometimes provide arrays that have only starting and ending
        # values. Between them values will be linear interpolated.
        doInterpolate = (None in self)
        for e in self:
            if not (type(e) in [float, int] or e is None):
                doInterpolate = False
                break
        if doInterpolate:
            xs, ys = [], []
            for i in range(0, len(self)):
                if self[i] is None: continue
                xs.append(i)
                ys.append(self[i])
            f = interpolate.interp1d(xs, ys)
            for i in range(0, len(self)):
                if self[i] is None:
                    self[i] = float(f(i))
                

    def __repr__(self):
        return "<Array of %d items>" % len(self)



class CFTurboObject:

    def __init__(self, node, checkType=True, defaultType='Object'):
        """Reads all children of a given node and attach them to the object
        as attributes."""
        # Get root node
        #if checkType:
             #assert node.attrib["Type"] in ["Object", "Vector2"]
        # Get all child nodes and convert
        self.__externalAttributes = []
        for childNode in list(node):
            name = childNode.tag
            convertedObj = convertXMLNode(childNode, defaultType=defaultType)

            if (
                "Name" in childNode.attrib or
                "Blade" in childNode.attrib or # dirty!
                "Index" in childNode.attrib
            ): # group as a dict
                if "Name" in childNode.attrib:
                    dictNameIndex = childNode.attrib["Name"]
                elif "Blade" in childNode.attrib:
                    dictNameIndex = childNode.attrib["Blade"]
                else:
                    dictNameIndex = None

                if "Index" in childNode.attrib:
                    dictIndex = childNode.attrib["Index"]
                else:
                    dictIndex = None

                if not hasattr(self, name):
                    setattr(self, name, {})
                self.__externalAttributes.append("%s->%s" % (name, dictNameIndex or dictIndex))
                if dictIndex is not None:
                    getattr(self, name)[dictIndex] = convertedObj
                if dictNameIndex is not None:
                    getattr(self, name)[dictNameIndex] = convertedObj
                    
            else:
                self.__externalAttributes.append(name)
                setattr(self, name, convertedObj)

    def __repr__(self):
        ret = []
        for k in self.__externalAttributes:
            if "->" in k:
                name, dictIndex = k.split("->")
                v = getattr(self, name)[dictIndex]
                k = "%s [\"%s\"]" % (name, dictIndex)
            else:
                v = getattr(self, k)
            if isinstance(v, CFTurboObject):
                suboutput = repr(v).split("\n")
                output = "\n".join([
                    ("%35s + " % k if i == 0 else " " * 35 + " | ") +\
                    suboutput[i]
                    for i in range(0, len(suboutput))
                ])
            else:
                output = "%35s = %s" % (k, repr(v))
            ret.append(output)
        return "\n".join(ret)



class YAMLObject:

    def __init__(self, d):
        self.parameters = d

    def __getattr__(self, key):
        return self.parameters[key]
