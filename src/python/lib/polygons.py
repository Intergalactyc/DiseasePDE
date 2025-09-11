from shapely.geometry import Polygon, MultiPolygon, Point, LinearRing
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath
import triangle

# https://www.cs.cmu.edu/~quake/triangle.poly.html

def _representative_point_from_linear_ring(lring: LinearRing):
    # When we look at the holes (polygon.interiors) of a polygon, they are LinearRings rather than Polygons
    # Because of this, we need to convert them to Polygons to find a representative interior point
    xs, ys = lring.xy
    poly = Polygon(shell = zip(xs, ys))
    hx, hy = poly.representative_point().coords[0]
    return hx, hy

def polygon_to_path(polygon):
    # Convert Shapely Polygon to Matplotlib Path
    vertices = list(polygon.exterior.coords)
    codes = [mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(vertices) - 2) + [mpath.Path.CLOSEPOLY]
    
    for hole in polygon.interiors:
        hole_vertices = list(hole.coords)
        vertices.extend(hole_vertices)
        codes.extend([mpath.Path.MOVETO] + [mpath.Path.LINETO] * (len(hole_vertices) - 2) + [mpath.Path.CLOSEPOLY])
    
    return mpath.Path(vertices, codes)

def save_polygon_as_poly(polygon: Polygon, filename: str, silent: bool = False):
    with open(filename, 'w') as f:
        subpolys = [polygon.exterior] + list(polygon.interiors)
        N = sum(len(poly.coords) - 1 for poly in subpolys) # Total number of vertices & equivalently of segments

        # First line: <# of vertices> <dimension (must be 2)> <# of attributes (0 here)> <# of boundary markers (0 here)> 
        f.write(f"{N} 2 0 1\n")

        # Following lines: <vertex #> <x> <y>
        i = 1 # vertex numbers (1-indexed)
        last_i = [1] # Keeps track of number of last vertex 
        for poly in subpolys:
            for x, y in list(poly.coords)[:-1]:
                f.write(f"{i} {x} {y} 1\n")
                i += 1
            last_i.append(i)
            
        # One line: <# of segments> <# of boundary markers (0 here)> 
        f.write(f"{N} 1\n")

        # Following lines: <segment #> <endpoint> <endpoint>
        i = 1 # segment numbers
        for j in range(len(last_i) - 1):
            for k in range(last_i[j], last_i[j+1] - 1):
                f.write(f"{i} {k} {k+1} {int(j==0)}\n")
                i += 1
            f.write(f"{i} {last_i[j+1] - 1} {last_i[j]} {int(j == 0)}\n") # Connect last point of subpoly to its first point
            i += 1

        # One line: <# of holes> 
        f.write(f"{len(polygon.interiors)}\n")

        # Following lines: <hole #> <representative x> <representative y> 
        for i, hole in enumerate(polygon.interiors, 1):
            hx, hy = _representative_point_from_linear_ring(hole)
            f.write(f"{i} {hx} {hy}\n")
    
    if not silent:
        print(f"Saved polygon to {filename}")

    return

def read_polygon_from_poly(filename: str, exterior_first: bool = False) -> Polygon:
    # Does not work for disconnected regions; does not maintain special attributes
    # If exterior_first is True, assumes that polygon containing vertex 1 (or 0) is exterior
    with open(filename, 'r') as f:
        # Function to handle reading next line
        def next_line():
            line = f.readline().strip()
            if line.startswith('#'): # Ignore comments
                return next_line()
            return line.split()

        # Get # of vertices from first line
        vertex_meta = next_line()
        if len(vertex_meta) != 4:
            raise Exception('Invalid vertex meta line')
        n_vertices = int(vertex_meta[0])
        if n_vertices < 3:
            raise Exception('Must be at least 3 vertices')
        if int(vertex_meta[1]) != 2:
            raise Exception('Only 2-dimensional polygons allowed')

        # Grab all of the vertices' coordinates, store in dict indexed by number
        vertices = dict()
        for _ in range(n_vertices):
            vertex = next_line()
            if len(vertex) < 3:
                raise Exception('Invalid vertex encountered')
            i, x, y = vertex[:3]
            vertices[i] = (float(x), float(y))
        
        # Get # of segments
        segs_meta = next_line()
        if len(segs_meta) != 2:
            raise Exception('Invalid segment meta line')
        n_segments = int(segs_meta[0])

        # Get all connections from segments and store as adjacency list
        connections = {i : [] for i in vertices.keys()}
        for _ in range(n_segments):
            segment = next_line()
            if len(segment) < 3:
                raise Exception('Invalid segment encountered')
            v1, v2 = segment[1:3]
            connections[v1].append(v2)
            connections[v2].append(v1)
        
        # Find distinct rings based on connections, and get corresponding coordinates
        rings = [[]]
        ring_i = 0
        exterior_ring = -1
        x, (v1, v2) = connections.popitem()
        while len(connections) > 0: # While there are unseen vertices...
            if exterior_first and (x == 0 or (x == 1 and exterior_ring == -1)):
                exterior_ring = ring_i
            rings[ring_i].append(vertices[x])
            if v1 in connections.keys():
                x, (v1, v2) = v1, connections.pop(v1)
            elif v2 in connections.keys():
                x, (v1, v2) = connections.pop(v2)
            else: # Done with this ring, so move on to another
                rings.append([])
                ring_i += 1
                x, (v1, v2) = connections.popitem()
        # Note that this only works if all rings are distinct (no intersections/s)

        # Get # of holes
        holes_meta = next_line()
        if len(holes_meta) != 1:
            raise Exception('Invalid hole meta line')
        n_holes = int(holes_meta[0])
        if n_holes > len(rings) - 1:
            raise Exception(f'Too few distinct rings found (expected {n_holes+1})')

        if not exterior_first:
            # Determine which rings describe holes
            rings_as_polys = [Polygon(shell = ring) for ring in rings]
            containments = [0 for _ in rings]
            for _ in range(n_holes):
                _, hx, hy = next_line()
                for i, poly in enumerate(rings_as_polys):
                    if poly.contains(Point(float(hx), float(hy))):
                        containments[i] += 1
            try:
                exterior_ring = containments.index(n_holes)
                if containments.count(1) != n_holes:
                    raise Exception('Hole descriptions invalid')
            except:
                raise Exception('No clear exterior detected')

        return Polygon(shell = rings[exterior_ring], holes = rings[:exterior_ring] + rings[(exterior_ring+1):])

def plot_polygon(poly: Polygon, ax=None, *, show_now=None, set_limits = False, color='tab:blue', boundary_color=None, fill=False, fill_color=None, background_fill_color = 'white', width = 1):
    if ax is None:
        if show_now is None:
            show_now = True
        _, ax = plt.subplots()

    if boundary_color is None:
        boundary_color = color
    if fill:
        if fill_color is None:
            fill_color = color
        ax.set_facecolor(background_fill_color)

    path = polygon_to_path(poly)
    patch = patches.PathPatch(path, facecolor=fill_color if fill else 'none', edgecolor=boundary_color, lw=width)
    ax.add_patch(patch)

    if set_limits:
        # Set plot limits
        ax.set_xlim(poly.bounds[0] - 1, poly.bounds[2] + 1)
        ax.set_ylim(poly.bounds[1] - 1, poly.bounds[3] + 1)
        ax.set_aspect('equal')

    if show_now:
        plt.tight_layout()
        plt.show()

def plot_gdf(gdf, ax = None, show_now = True):
    if ax is None:
        fig, ax = plt.subplots()

    for poly in gdf['geometry']:
        if type(poly) is Polygon:
            plot_polygon(poly, ax, show_now = False)
        elif type(poly) is MultiPolygon:
            for p in poly.geoms:
                plot_polygon(p, ax, show_now = False)
    if show_now:
        plt.show()

class TriPolygon:
    def __init__(self, points, segments, holes):
        self._points = points
        self._segments = segments
        self._holes = holes

    def triangle(self) -> triangle.Triangle:
        t = triangle.Triangle()
        t.set_points(self._points)
        t.set_segments(self._segments)
        t.set_holes(self._holes)
        return t

    @classmethod
    def from_polygon(cls, poly: Polygon):
        points, segments, holes = [], [], []

        exterior_points = list(zip(*poly.exterior.xy))[:-1]
        N_exterior = len(exterior_points)
        points += exterior_points
        segments += [(i, i+1) for i in range(N_exterior - 1)] + [(N_exterior-1, 0)]
        
        i_current = N_exterior
        
        for interior in poly.interiors:
            interior_points = list(zip(*interior.xy))[:-1]
            N_interior = len(interior_points)
            points += interior_points
            segments += [(i_current + i, i_current + i + 1) for i in range(N_interior - 1)] + [(i_current + N_interior - 1, i_current)]
            i_current += N_interior
            holes.append(_representative_point_from_linear_ring(interior))

        return cls(points, segments, holes)
