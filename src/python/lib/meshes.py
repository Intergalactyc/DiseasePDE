import shapely
import lib.polygons as polygons
import meshio
import matplotlib.pyplot as plt
import matplotlib.tri as tri

class IOMesh(meshio.Mesh):
    def plot(self, ax = None):
        if ax is None:
            _, ax = plt.subplots()

        points = self.points
        triangles = self.cells_dict['triangle']
        x = [point[0] for point in points]
        y = [point[1] for point in points]

        plttri = tri.Triangulation(x, y, triangles)
        ax.triplot(plttri, linewidth = 0.1, c = 'tab:orange')

def mesh_polygon(polygon: shapely.Polygon|polygons.TriPolygon, area: float, mode: str = "zpq30"):
    if isinstance(polygon, shapely.Polygon):
        polygon = polygons.TriPolygon.from_polygon(polygon)
    if not isinstance(polygon, polygons.TriPolygon):
        raise Exception(f"Invalid polygon {polygon} of type {type(polygon)}")

    t = polygon.triangle()
    t.triangulate(area = area, mode = mode)

    points = [point[0] for point in t.get_points()]
    triangles = [triangle[0] for triangle in t.get_triangles()]

    cells = [("triangle", triangles)]
    mesh = IOMesh(points, cells)

    return mesh
