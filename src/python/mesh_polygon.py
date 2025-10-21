import os
from lib.paths import DATA_PRODUCTS
from lib.meshes import mesh_polygon
import lib.polygons as polygons
import argparse
import matplotlib.pyplot as plt

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", type = str, default = None)
    parser.add_argument("--area", type = float, default = None)
    parser.add_argument("--show", action = "store_true", default = False)
    parser.add_argument("-o", "--output-name", type = str, default = None)

    args = vars(parser.parse_args())

    return {"name" : args.get("state") or "us", "area" : args.get("area") or (0.001 if args.get("state") else 0.01), "show" : args.get("show")}

def main():
    args = parse()
    name = args["name"]
    area = args["area"]
    out = args["output-name"] or name

    poly_path = os.path.join(DATA_PRODUCTS, "polygons", f"{name}.poly")
    if not os.path.exists(poly_path):
        raise OSError(f"File {poly_path} does not exist. Run `get_polygons.py` first.")
    polygon = polygons.read_polygon_from_poly(poly_path)
    mesh = mesh_polygon(polygon, area)

    if args["show"]:
        mesh.plot()
        plt.show()

    mesh.write(os.path.join(DATA_PRODUCTS, "meshes", f"{out}.exo"))

if __name__ == "__main__":
    main()

