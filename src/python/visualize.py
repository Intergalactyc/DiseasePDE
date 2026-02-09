import meshio
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.tri as tri
from matplotlib.colors import SymLogNorm, Normalize
import argparse
import pathlib
# from glob import glob

# note that this only works currently with cell (not point) data

PRODUCTS = pathlib.Path(__file__).parent.parent.parent / "data_products"
SIMD = PRODUCTS / "simulated"
DATA = PRODUCTS / "datameshes"
TEST = PRODUCTS / "test_results"

class Result:
    def __init__(self, filepaths):
        # assumes all have the same mesh just with different data
        self.timesteps = []
        self.cell_data = {}
        self.maxes = {}
        self.mins = {}
        for i, f in enumerate(filepaths):
            # parse timestep
            t = float(str(f)[:-4].split("/")[-1].split("_")[-1])
            self.timesteps.append(t)
            # load file
            data = meshio.read(f)
            # using the first file loaded, take the actual mesh
            if i == 0:
                self.points = data.points
                self.x = self.points[:, 0]
                self.y = self.points[:, 1]
                self.triangles = data.cells_dict["triangle"]
            # load the cell data for each file
            block_index = list(data.cells_dict).index("triangle")
            self.cell_data[t] = {
                k : data.cell_data[k][block_index]
                for k in data.cell_data.keys()
            }
            for k, v in self.cell_data[t].items():
                if k not in self.mins:
                    self.mins[k] = 1e16
                if k not in self.maxes:
                    self.maxes[k] = 0.
                if (m:=v.min()) < self.mins[k]:# and m > 0:
                    self.mins[k] = m
                if (M:=v.max()) > self.maxes[k]:
                    self.maxes[k] = M
        self.timesteps.sort()

    def plot(self, compartment, index, cmap="coolwarm"):
        # index should be time step index
        time = self.timesteps[index]
        data_t = self.cell_data[time]
        if compartment.lower().endswith("capita"):
            vals = data_t["I"]/sum(data_t[c] for c in data_t.keys())
        elif compartment.lower() == "total":
            vals = sum(data_t[c] for c in data_t.keys())
        else:
            vals = data_t[compartment]
        
        pos = vals > 0
        vmin = vals[pos].min()
        vmax = vals[pos].max()
        norm = SymLogNorm(linthresh=1e-6, vmin=vmin, vmax=vmax)

        plttri = tri.Triangulation(self.x, self.y, self.triangles)

        fig, ax = plt.subplots(figsize=(9, 6))

        tpc = ax.tripcolor(
            plttri,
            facecolors=vals,
            shading="flat",
            norm=norm,
            cmap=cmap
        )

        cbar = fig.colorbar(tpc, ax=ax)
        cbar.set_label(f"{compartment} (log scale)")

        ax.set_title(f"Time: {time} (index {index})")

        plt.tight_layout()
        plt.show(block=False)

    def animate(self, compartment, cmap="coolwarm", log=True):
        def get_vals(i):
            time = self.timesteps[i]
            data_t = self.cell_data[time]
            if compartment.lower().endswith("capita"):
                return time, data_t["I"]/sum(data_t[c] for c in data_t.keys())
            elif compartment.lower() == "total":
                return time, sum(data_t[c] for c in data_t.keys())
            return time, data_t[compartment]

        t0, vals = get_vals(0)

        norm = SymLogNorm(linthresh=1e-6, vmin=0, vmax=1) if log else Normalize(vmin=0, vmax=1)
        # norm = SymLogNorm(linthresh=1e-6, vmin=self.mins[compartment], vmax=self.maxes[compartment])

        plttri = tri.Triangulation(self.x, self.y, self.triangles)

        fig, ax = plt.subplots(figsize=(9, 6))

        tpc = ax.tripcolor(
            plttri,
            facecolors=vals,
            shading="flat",
            norm=norm,
            cmap=cmap
        )
        # FuncAnimation with tpc.set_facecolors()
        cbar = fig.colorbar(tpc, ax=ax)
        cbar.set_label(f"{compartment} (log scale)" if log else compartment)

        label = ax.text(0, 0, f"Time: {t0} (index 0)", ha="left", va="center", fontsize=8, transform=ax.transAxes)

        def update(frame):
            t, vals = get_vals(frame)
            label.set_text(f"Time: {t} (index {frame})")
            tpc.set_array(vals)

        anim = FuncAnimation(fig, update, interval=200, frames=len(self.timesteps))

        plt.tight_layout()
        plt.show()

    def interact(self):
        try:
            while True:
                try:
                    c = input("Compartment: ")
                    i = int(input("Time index: "))
                    self.plot(c, i)
                    print("")
                except Exception as e:
                    if isinstance(e, KeyboardInterrupt):
                        raise
                    print("Invalid input")
        except KeyboardInterrupt:
            print("Exiting.")
            exit()

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--compartment",type=str,default="per capita",help="compartment to plot")
    parser.add_argument("-l","--linear",action="store_true",help="use linear (not logarithmic) scaling for animation")
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parse()

    matches = SIMD.glob("L2SIR_*.exo")
    # matches = TEST.glob("SIR_*.exo")
    # matches = SIMD.glob("ALT_*.exo")
    # matches = SIMD.glob("NEW_*.exo")
    # matches = [str(DATA.joinpath(f"SIR_{t}.exo")) for t in range(100,1000,20)]#glob("SIR_*.exo")
    
    res = Result(matches)
    res.animate(args["compartment"], log=not args["linear"])
