import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import SymLogNorm
import argparse
import pathlib
# from glob import glob

SIMD = pathlib.Path(__file__).parent.parent.parent / "data_products" / "simulated"
DATA = pathlib.Path(__file__).parent.parent.parent / "data_products" / "datameshes"

class Result:
    def __init__(self, filepaths):
        # assumes all have the same mesh just with different data
        self.timesteps = []
        self.cell_data = {}
        for i, f in enumerate(filepaths):
            # parse timestep
            t = float(str(f)[:-4].split("/")[-1].split("_")[-1])
            self.timesteps.append(t)
            # load file
            data = meshio.read(f)
            # for first file loaded, save the actual mesh
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
        plt.show()
        

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--compartment",type=str,default="per capita",help="compartment to plot")
    parser.add_argument("-t","--timestep",type=int,default=0,help="index of timestep to plot")
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = parse()
    matches = SIMD.glob("SIR_*.exo")#[str(DATA.joinpath(f"SIR_{t}.exo")) for t in range(100,1000,20)]#glob("SIR_*.exo")
    res = Result(matches)
    res.plot(args["compartment"], args["timestep"])
    
