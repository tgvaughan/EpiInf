# Basic phylogenetic tree manipulation module

import re

class Node:
    """A node in a phylogenetic tree."""

    def __init__(self):
        self.height = None
        self.time = None
        self.branchLength = None
        self.parent = None
        self.children = []
        self.annotations = {}
        self.label = None

    def isRoot(self):
        return self.parent == None

    def isLeaf(self):
        return len(self.children)==0

    def addChild(self, newChild):
        self.children.append(newChild)
        newChild.parent = self

    def getAllChildren(self):
        """Retrieve a list including this node and all of its decendents."""

        childList = [self]
        for child in self.children:
            childList.extend(child.getAllChildren())

        return childList

    def getLeaves(self):
        """Retrieve the list of leaves descending from this node."""

        if self.isLeaf():
            return [self]
        else:
            leaves = []
            for child in self.children:
                leaves.extend(child.getLeaves())
            return leaves

    def getNewick(self):
        """Retrieve a Newick representation of the tree below this node."""

        newick = ""
        if len(self.children)>0:
            newick += "("
            for i,v in enumerate(self.children):
                if i>0:
                    newick += ","
                newick += v.getNewick()
            newick += ")"

        if self.label != None:
            newick += self.label

        if len(self.annotations)>0:
            newick += '[&'
            isFirst = True
            for k,v in self.annotations.items():
                if isFirst:
                    isFirst = False
                else:
                    newick += ","
                newick += '{}="{}"'.format(k,v)
            newick += ']'

        if self.parent == None:
            if self.origin != None:
                newick += ":{}".format(self.origin - self.height)
            else:
                newick += ":0.0"
        else:
            newick += ":{}".format(self.parent.height - self.height)

        return newick

    def computeTimes(self, offset):
        """Set up time attribute of this node and each node below it.  The offset
        parameter indicates the time of this node's parent."""

        self.time = offset + self.branchLength
        for child in self.children:
            child.computeTimes(self.time)


class Tree:
    """A phylogenetic tree."""

    root = None

    def __init__(self, arg):
        if type(arg) is Node:
            self.root = arg
            return

        # Need to parse string or string from file
        newickString = ""
        if type(arg) is str:
            newickStr = arg

        if type(arg) is file:
            firstLine = arg.readline()

            newickString = ""
            if firstLine.lower().startswith('#nexus'):
                for line in arg.readlines():
                    if line.strip().lower().startswith("tree "):
                        newickString = line[(line.find("=")+1):].strip()
                        break
            else:
                newickString = firstLine.strip()

        self.root = self.loadFromString(newickString)


    def __repr__(self):
        return self.root.getNewick() + ";"


    # Basic tree queries

    def getNodes(self):
        return self.root.getAllChildren()

    def getLeaves(self):
        return self.root.getLeaves()

    # Tree parsing code

    class ParseError(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

    class ParseContext:
        def __init__(self, tokenList, valueList):
            self.tokenList = tokenList
            self.valueList = valueList
            self.idx = 0

        def acceptToken(self, token, manditory=False):
            if self.tokenList[self.idx] == token:
                self.idx += 1
                return True
            else:
                if not manditory:
                    return False
                else:
                    raise ParseError("Error parsing token {} ({})".format(self.tokenList[self.idx], self.valueList[self.idx]))

        def getLastValue(self):
            return self.valueList[self.idx-1]


    def loadFromString(self, string):

        # Lexical analysis

        tokens = [
                ('LPAREN',  '\('),
                ('RPAREN',  '\)'),
                ('COLON',   ':'),
                ('STRING', '"[^"]*"'),
                ('STRING', '\'[^\']*\''),
                ('STRING',   '[a-zA-Z0-9_.-]+'),
                ('OPENA', '\[&'),
                ('EQUALS', '='),
                ('CLOSEA', '\]'),
                ('COMMA',   ','),
                ('SEMI',    ';')
                ]

        idx = 0
        tokenList = []
        valueList = []

        while idx < len(string):

            noMatch = True

            for k in range(len(tokens)):
                match = re.match(tokens[k][1], string[idx:])

                if match != None:
                    tokenList.append(tokens[k][0])
                    idx += len(match.group(0))

                    if tokens[k][0] == 'STRING':
                        value = match.group(0)
                        if len(value)>2:
                            if value[0]=='"' and value[len(value)-1]=='"':
                                value = value[1:(len(value)-1)]
                            else:
                                if value[0]=="'" and value[len(value)-1]=="'":
                                    value = value[1:(len(value)-1)]
                        valueList.append(value)
                    else:
                        valueList.append(None)

                    noMatch = False
                    break

            if noMatch:
                raise Tree.ParseError('Unrecognized character at position ' + str(idx) + ': \'' + string[idx] + '\'')

        # Recursive decent parser

        ctx = Tree.ParseContext(tokenList, valueList)
        root = self.ruleN(None, ctx, 0)
        ctx.acceptToken("SEMI", manditory=True)

        # Compute node heights and times:

        root.computeTimes(0.0)

        maxTime = 0.0
        for node in root.getAllChildren():
            maxTime = max(maxTime, node.time)

        for node in root.getAllChildren():
            node.height = maxTime - node.time

        root.origin = root.height + root.branchLength

        return root

    def ruleN(self, parent, ctx, depth):
        node = Node()
        if parent != None:
            parent.addChild(node)

        self.ruleS(node, ctx, depth)
        self.ruleL(node, ctx, depth)
        self.ruleA(node, ctx, depth)
        self.ruleB(node, ctx, depth)

        #print " "*depth + str(node.label) + ":" + str(node.branchLength)
        return node 

    def ruleS(self, node, ctx, depth):
        if ctx.acceptToken('LPAREN'):
            #print " "*depth + "("
            self.ruleN(node, ctx, depth+1)
            self.ruleQ(node, ctx, depth)
            ctx.acceptToken('RPAREN', manditory=True)
            #print " "*depth + ")"

    def ruleQ(self, node, ctx, depth):
        if ctx.acceptToken('COMMA'):
            self.ruleN(node, ctx, depth+1)
            self.ruleQ(node, ctx, depth)

    def ruleL(self, node, ctx, depth):
        if ctx.acceptToken('STRING'):
            node.label = ctx.getLastValue()

    def ruleA(self, node, ctx, depth):
        if ctx.acceptToken('OPENA'):
            self.ruleC(node, ctx, depth)
            self.ruleD(node, ctx, depth)
            ctx.acceptToken('CLOSEA', manditory=True)

    def ruleC(self, node, ctx, depth):
        ctx.acceptToken('STRING', manditory=True)
        key = ctx.getLastValue()
        ctx.acceptToken('EQUALS', manditory=True)
        ctx.acceptToken('STRING', manditory=True)
        value = ctx.getLastValue()

        node.annotations[key] = value

    def ruleD(self, node, ctx, depth):
        if ctx.acceptToken('COMMA'):
            self.ruleC(node, ctx, depth)
            self.ruleD(node, ctx, depth)

    def ruleB(self, node, ctx, depth):
        if ctx.acceptToken('COLON'):
            ctx.acceptToken('STRING', manditory=True)
            node.branchLength = float(ctx.getLastValue())


    # Plotting code

    def plot(self, **kwargs):

        from matplotlib.pylab import plot, text, gca, xlim, ylim
        import scipy as sp

        if "color" not in kwargs.keys():
            kwargs['color'] = 'black'
        
        leaves = self.getLeaves()
        nodes = self.getNodes()
        pos = sp.zeros(len(nodes))
        
        def computePos(node):
            idx = nodes.index(node)
            if node.isLeaf():
                pos[idx] = leaves.index(node)
            else:
                for child in node.children:
                    pos[idx] += computePos(child)
                pos[idx] /= len(node.children)

            return pos[idx]
        computePos(self.root)

        for i in range(len(nodes)):

            if nodes[i].parent == None:
                plot([nodes[i].height, nodes[i].height+nodes[i].branchLength],
                        [pos[i], pos[i]], **kwargs)
            else:
                pi = nodes.index(nodes[i].parent)
                plot([nodes[i].height, nodes[i].parent.height, nodes[i].parent.height],
                        [pos[i], pos[i], pos[pi]], **kwargs)

            if nodes[i].label != None:
                text(nodes[i].height, pos[i], nodes[i].label, {'ha':'left', 'va':'bottom'}, rotation=45)

        ylim(-0.05*len(leaves), 1.05*len(leaves)-1)
        #xlim(-0.05*self.root.origin, 1.05*self.root.origin)

        ax = gca()
        ax.invert_xaxis()
        ax.yaxis.set_visible(False)
        ax.set_frame_on(False)
        ax.grid()
