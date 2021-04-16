# ------------------------------------------------------------------------------
# 1.1 Import modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 14

# scientific format of axis labels
def sci_format():
    fmt = mpl.ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2, 3))
    return (fmt)

# ------------------------------------------------------------------------------
# 1.2 Define analyzer
class Analyze(object):

    def __init__(self,
                 root = ".",  # path to result data
                 x_scale = 1.0,
                 y_scale = 1.0,
                 u_scale = 1.0
                 ):
        self = locals().pop("self")
        for name, val in locals().items():
            print("setting " + name + " to:", val)
            setattr(self, name, val)

        # load all data
        self.load_source()
        self.load_target()
        self.load_error()

    def load_source(self):
        # import points
        self.source_x, self.source_y = np.loadtxt(self.root + '/source.csv', delimiter=',', unpack=True)
        self.source_x *= self.x_scale
        self.source_y *= self.y_scale
        # deduce axis range
        self.xlim_source = [self.source_x.min(), self.source_x.max()]
        self.ylim_source = [self.source_y.min(), self.source_y.max()]
        # import field and unscale
        self.source = np.loadtxt(self.root+'/source.dat', delimiter=',', unpack=True)
        self.source *= self.u_scale

    def load_target(self):
        from scipy.spatial.transform import Rotation as R
        # import points
        self.target_x, self.target_y = np.loadtxt(self.root + '/target.csv', delimiter=',', unpack=True)
        # deduce axis range
        self.xlim_target = [self.target_x.min(), self.target_x.max()]
        self.ylim_target = [self.target_y.min(), self.target_y.max()]
        # fields
        self.remap = np.loadtxt(self.root+'/remap.dat', delimiter=',', unpack=True)

    def load_error(self):
        self.error = np.loadtxt(self.root+'/error.dat', delimiter=',', unpack=True)
        self.exact = np.loadtxt(self.root+'/exact.dat', delimiter=',', unpack=True)

# ------------------------------------------------------------------------------
# optional scaling of coordinates and field
x_scale = 1e-6
y_scale = 1e-4
u_scale = 1e-8

# create an instance
analysis = Analyze(".", x_scale, y_scale, u_scale)

# 2.1 plot wavelets and mesh points
fig, axes = plt.subplots(1,1, figsize=[8,6])
ax=axes
ax.set_title("source & target")
ax.scatter(analysis.source_x, analysis.source_y, s=.1, facecolor='C0', label="source")
ax.scatter(analysis.target_x, analysis.target_y, s=.1, facecolor='C1', label="target")
ax.set_xlim(analysis.xlim_target)
ax.set_ylim(analysis.ylim_target)
ax.set_xlabel(r'$x\prime$')
ax.set_ylabel(r'$y\prime$')
plt.legend(loc="upper right")
print("plotting source and target points")
plt.savefig("points.png")

# ------------------------------------------------------------------------------
# 2.2 plot source, exact, remapped and error map
fig, axes = plt.subplots(2,2, figsize=[16,12])

ax=axes[0,0]
ax.set_title("source field")
sc=ax.scatter(analysis.source_x, analysis.source_y, c=analysis.source, marker='o', s=5, cmap="coolwarm")
ax.set_xlim(analysis.xlim_target)
ax.set_ylim(analysis.ylim_target)
fig.colorbar(sc, ax=ax)
ax.set_xlabel(r'$x\prime$')
ax.set_ylabel(r'$y\prime$')

ax=axes[0,1]
ax.set_title("remapped field")
sc=ax.scatter(analysis.target_x, analysis.target_y, c=analysis.remap, marker='o', s=5, cmap="coolwarm")
ax.set_xlim(analysis.xlim_target)
ax.set_ylim(analysis.ylim_target)
fig.colorbar(sc, ax=ax)
ax.set_xlabel(r'$x\prime$')
ax.set_ylabel(r'$y\prime$')

analysis.load_error()
ax=axes[1,0]
ax.set_title("exact values")
sc=ax.scatter(analysis.target_x, analysis.target_y, c=analysis.exact, marker='o', s=5, cmap="coolwarm")
ax.set_xlim(analysis.xlim_target)
ax.set_ylim(analysis.ylim_target)
fig.colorbar(sc, ax=ax)
ax.set_xlabel(r'$x\prime$')
ax.set_ylabel(r'$y\prime$')

ax=axes[1,1]
ax.set_title("error map")
sc=ax.scatter(analysis.target_x, analysis.target_y, c=analysis.error, marker='.', s=5, cmap="coolwarm")
ax.set_xlim(analysis.xlim_target)
ax.set_ylim(analysis.ylim_target)
fig.colorbar(sc, ax=ax)
ax.set_xlabel(r'$x\prime$')
ax.set_ylabel(r'$y\prime$')
print("plotting source, exact, remapped fields and error map")
plt.savefig("fields.png")
# ------------------------------------------------------------------------------