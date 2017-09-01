import math
import numpy
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True,
                 "font.size": 20})
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, art3d

# hidden line style
hlw = 1
hec = "grey"
hls = ":"

# node label position offset
noffset = 0.02

# cell params
cParams = {"linewidths":2,
           "edgecolors":"k",
           "alpha":0.0}

# label parameters
laParams = {"ha":"left",
            "va":"top",
            "weight":"bold",
            "backgroundcolor":"w"}
flaParams = {"ha":"center",
            "va":"center",
            "weight":"bold",}


class BaseEntity(object):
    num_entities = 0
    def __init__(self):
        self.id = self.__class__.num_entities
        self.__class__.num_entities += 1

    def __repr__(self):
        return "%s %d" % (self.__class__.__name__, self.id)

    def __hash__(self):
        return self.__repr__()

class Node(BaseEntity):
    def __init__(self, x, y, z):
        super(Node, self).__init__()
        self.x = x
        self.y = y
        self.z = z

    def placeLabel(self, **kwargs):
        plt.gca().text(self.x+noffset, self.y, self.z-noffset,
                       "%d" % self.id,
                       **kwargs)


class Face(BaseEntity):
    _cache = []

    @classmethod
    def _getCache(cls, nodes):
        for face in Face._cache:
            if sorted(face.nodes) == sorted(nodes):
                return face
        return None

    def __new__(cls, nodes, *args, **kwargs):
        existingFace = cls._getCache(nodes)
        if existingFace:
            return existingFace
        face = super(Face, cls).__new__(cls, nodes, *args, **kwargs)
        return face

    def __init__(self, nodes, plane="x", hidden=False):
        if self in self._cache:
            return
        self._cache.append(self)

        super(Face, self).__init__()

        self.nodes = nodes
        self.plane = plane
        self.hidden = hidden

    def place3d(self, **kwargs):
        override = kwargs.get("override", False)
        if override: del kwargs["override"]

        if self.hidden and not override:
            kwargs.update({"linestyles": hls,
                           "linewidths": hlw,
                           "edgecolors": hec,
                           "zorder": -1})
        self._poly = art3d.Poly3DCollection([[(n.x, n.y, n.z)
                                              for n in self.nodes]],
                                            **kwargs)
        # need to do this here because of weird issue with facecolor and alpha
        self._poly.set_facecolor("w")
        plt.gca().add_collection3d(self._poly)

    _center = None
    @property
    def center(self):
        if not self._center:
            self._center = [numpy.mean([n.x for n in self.nodes]),
                            numpy.mean([n.y for n in self.nodes]),
                            numpy.mean([n.z for n in self.nodes])]
        return self._center

    def placeLabel(self, **kwargs):
        kwargs.update({"zdir": self.plane})
        plt.gca().text(self.center[0], self.center[1], self.center[2],
                       "%d" % self.id, **kwargs)

class Cell(BaseEntity):
    def __init__(self, nodelist, is_front, is_right, is_top):
        super(Cell, self).__init__()
        # assumes normalized ordering of
        #  (0,0,0), (1,0,0), (1,1,0), (0,1,0) (0,0,1), (1,0,1), (1,1,1), (0,1,1)
        self.nodes = nodelist
        # this tells us if we are on any of the visible boundaries
        self.is_front = is_front
        self.is_right = is_right
        self.is_top = is_top

        self._construct_faces()

    def _construct_faces(self):
        self.faces = [
            # front - hidden for cells 2,3,6,7
            Face(self.nodes[:2] + self.nodes[5:3:-1], "x",
                 hidden=(not self.is_front)),
            # right - hidden for cells 0,2,4,6
            Face(self.nodes[1:3] + self.nodes[6:4:-1], "y",
                 hidden=(not self.is_right)),
            # back - all hidden
            Face(self.nodes[2:4] + self.nodes[7:5:-1], "x", hidden=True),
            # left - all hidden
            Face([self.nodes[3], self.nodes[0], self.nodes[4], self.nodes[7]],
                 "y",
                 hidden=True),
            # bottom - all hidden
            Face(self.nodes[:4], "x", hidden=True),
            # top - hidden for cells 0,1,2,3
            Face(self.nodes[4:], "x", hidden=(not self.is_top))
        ]

    def place3d(self, **kwargs):
        for f in self.faces:
            f.place3d(**kwargs)

    def draw3d(self, **kwargs):
        for f in self.faces:
            f.place3d(override=True, **kwargs)

    def placeNodes(self, **kwargs):
        xs = [n.x for n in self.nodes]
        ys = [n.y for n in self.nodes]
        zs = [n.z for n in self.nodes]
        plt.gca().scatter(xs, ys, zs, c="k", s=60)
        for i, (x,y,z) in enumerate(zip(xs, ys, zs)):
            plt.gca().text(x-noffset, y, z+noffset, "%d" % i, **kwargs)


class Cube(object):
    _xyfaces = range(0,12)
    _xzfaces = range(12, 24)
    _yzfaces = range(24, 36)
    def __init__(self, nx, ny, nz):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.ncells = nx*nx*nz

        self._build_nodes()

        self._build_cells()

        self.faces = [f for c in self.cells for f in c.faces]

        self._reorder_facenums()

    def _build_nodes(self):
        dx = 1./self.nx; dy = 1./self.ny; dz = 1./self.nz
        self.nodes = []
        for z in numpy.linspace(0.0, 1.0, self.nz+1):
            for y in numpy.linspace(0.0, 1.0, self.ny+1):
                for x in numpy.linspace(0.0, 1.0, self.nx+1):
                    self.nodes.append(Node(x,y,z))

    def _build_cells(self):
        # useful offsets
        yoff = self.nx+1
        zoff = yoff*(self.ny+1)

        # build the cells
        self.cells = []
        for z in range(self.nz):
            for y in range(self.ny):
                for x in range(self.nx):
                    # lower left (front) node of each cell
                    lln = x + y*yoff + z*zoff
                    self.cells.append(Cell(
                        # bottom face - counter-clockwise node orientation
                        # (from +z)
                        self.nodes[lln:lln+2]
                        + self.nodes[lln+yoff+1:lln+yoff-1:-1]
                        # top face  - counter-clockwise node orientation
                        # (from +z)
                        + self.nodes[lln+zoff:lln+zoff+2] +
                        self.nodes[lln+yoff+zoff+1:lln+yoff+zoff-1:-1],
                        # are we on visible boundaries of the domain?
                        is_front=(y == 0),
                        is_right=(x == (self.nx-1)),
                        is_top=(z == (self.nz-1))))

    def _reorder_facenums(self):
        # global face ids are not in the order we added them...
        # doing redundant operations, but oh well
        newIds = [12, 25, 14, 24, 0, 4,   # cell 0
                  13, 26, 15, 25, 1, 5,   # cell 1
                  14, 28, 16, 27, 2, 6,   # cell 2
                  15, 29, 17, 28, 3, 7,   # cell 3
                  18, 31, 20, 30, 4, 8,   # cell 4
                  19, 32, 21, 31, 5, 9,   # cell 5
                  20, 34, 22, 33, 6, 10,  # cell 6
                  21, 35, 23, 34, 7, 11,  # cell 7
        ]
        for ffid, f in zip(newIds, self.faces):
            f.id = ffid

    def draw3d(self, **kwargs):
        for c in self.cells:
            c.place3d(**kwargs)

    def placeNodes(self, **kwargs):
        xs = [n.x for n in self.nodes]
        ys = [n.y for n in self.nodes]
        zs = [n.z for n in self.nodes]
        plt.gca().scatter(xs, ys, zs, **kwargs)

    def placeNodeLabels(self, **kwargs):
        for n in self.nodes:
            n.placeLabel(**kwargs)

    def placeFaceLabels(self, plane="all", **kwargs):
        if plane == "all":
            for f in self.faces:
                f.placeLabel(**kwargs)
        else:
            fids = getattr(self, "_%sfaces" % plane, None)
            for f in filter(lambda f: f.id in fids, self.faces):
                f.placeLabel(**kwargs)

fig = plt.figure()
ax3d = fig.add_subplot(111, projection="3d")

c = Cube(2,2,2)

c.draw3d(**cParams)
c.placeNodes(c="k", s=60)
c.placeNodeLabels(**laParams)

plt.axis("off")
ax3d.margins(0.01)
#plt.show()
plt.savefig("allNodes.svg", bbox_inches="tight")


for plane in ["xy", "xz", "yz"]:
    plt.cla(); plt.clf();

    fig = plt.figure()
    ax3d = fig.add_subplot(111, projection="3d")

    c.draw3d(**cParams)
    c.placeFaceLabels(plane, **flaParams)

    plt.axis("off")
#    plt.show()
    plt.savefig("%s_plane.svg" % plane, bbox_inches="tight")



plt.cla(); plt.clf();

fig = plt.figure()
ax3d = fig.add_subplot(111, projection="3d")

c.cells[0].draw3d(**cParams)
c.cells[0].placeNodes(**laParams)

plt.axis("off")
ax3d.autoscale_view()
ax3d.margins(0.01)
plt.show()
