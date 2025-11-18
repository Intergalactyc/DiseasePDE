import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import pathlib

# Good color maps: viridis, inferno, plasma, magma, coolwarm
CMAP = "coolwarm"

# Load mesh ---------------------------------------------------------------
PATH = (
    pathlib.Path(__file__).parent.parent.parent
    / "data_products" / "datameshes" / "SIR_100.exo"
)

mesh = meshio.read(PATH)

# Geometry ---------------------------------------------------------------
points = mesh.points
triangles = mesh.cells_dict["triangle"]
x = points[:, 0]
y = points[:, 1]

# Cell data: S compartment -----------------------------------------------
# mesh.cell_data[name] returns a list of arrays, one per cell block
triangle_block_index = list(mesh.cells_dict).index("triangle")
S_vals = mesh.cell_data["S"][triangle_block_index]

# LogNorm requires positive values
positive = S_vals > 0
vmin = S_vals[positive].min()
vmax = S_vals.max()
norm = LogNorm(vmin=vmin, vmax=vmax)

# Plot -------------------------------------------------------------------
plttri = tri.Triangulation(x, y, triangles)

fig, ax = plt.subplots(figsize=(6, 6))

tpc = ax.tripcolor(
    plttri,
    facecolors=S_vals,
    shading="flat",
    norm=norm,
    cmap=CMAP
)

cbar = fig.colorbar(tpc, ax=ax)
cbar.set_label("S (log scale)")

# ax.set_aspect("equal")
plt.tight_layout()
plt.show()
