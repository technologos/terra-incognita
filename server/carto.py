import numpy as np
from scipy.spatial import Voronoi

from matplotlib import pyplot as plt
plt.style.use('default')

import logging
logger = logging.getLogger(__name__)

class Map():
    def __init__(self, 
                 pixel_width: int, 
                 pixel_height: int,
                 pixels_per_tile: int):
        self.width = pixel_width
        self.height = pixel_height
        self.pixels_per_tile = pixels_per_tile
        self.number_of_tiles = (self.width // self.pixels_per_tile) * (self.height // self.pixels_per_tile)
        self.all_tiles = []
        self.tiles = []
        self.border_tiles = []
        self.vertices = []
        self.edges = []
        self.plates = []

    def __str__(self):
        return "Map of dimensions {} x {}".format(self.width, self.height)
    

    def setup(self, random_seed: int = 42) -> None:

        logger.info('Generating random points')
        np.random.seed(random_seed)

        points = set()

        # Generate unique points
        while len(points) < self.number_of_tiles:
            new_points = zip(
                np.random.randint(0, self.width, self.number_of_tiles - len(points)),
                np.random.randint(0, self.height, self.number_of_tiles - len(points)),
            )
            points.update(new_points)

        # Unpack the unique points into x and y arrays
        x_values, y_values = zip(*list(points))
        x_values = np.array(x_values)
        y_values = np.array(y_values)

        # create points evenly spaced around the perimeter of a rectangle with dimensions (width+4*pixels_per_tile, height+4*pixels_per_tile)
        logger.info('Generating border points')
        n_width_border_points = (self.width // self.pixels_per_tile) + 4
        n_height_border_points = (self.height // self.pixels_per_tile) + 4
        border_x_values = np.linspace(-2*self.pixels_per_tile, 
                                      self.width + 2*self.pixels_per_tile, 
                                      n_width_border_points)
        border_y_values = np.linspace(-2*self.pixels_per_tile, 
                                      self.height + 2*self.pixels_per_tile, 
                                      n_height_border_points)

        # create a Voronoi diagram with the tile points, plus the bottom, top, left, and right border points
        logger.info('Creating Voronoi diagram')
        vor = Voronoi(np.column_stack((np.concatenate((x_values, 
                                                       border_x_values, 
                                                       border_x_values, 
                                                       np.repeat(-2*self.pixels_per_tile, n_height_border_points), 
                                                       np.repeat(self.width + 2*self.pixels_per_tile, n_height_border_points))), 
                                       np.concatenate((y_values,
                                                       np.repeat(-2*self.pixels_per_tile, n_width_border_points), 
                                                       np.repeat(self.height + 2*self.pixels_per_tile, n_width_border_points),
                                                       border_y_values, 
                                                       border_y_values)))))
        logger.info('Voronoi diagram created')
        
        # create a list of all tiles
        self.all_tiles = [None] * len(vor.points)
        for i in range(len(vor.points)):
            self.all_tiles[i] = Tile(self, 
                                 center_x = vor.points[i][0], 
                                 center_y = vor.points[i][1],
                                 core = (i < self.number_of_tiles))

        # downselect to the tiles that are within the bounds of the map
        self.tiles = self.all_tiles[:self.number_of_tiles]

        # create a list of the border tiles
        self.border_tiles = self.all_tiles[self.number_of_tiles:]

        # create a list of the vertices of the Voronoi diagram
        self.vertices = [None] * len(vor.vertices)
        for i in range(len(vor.vertices)):
            self.vertices[i] = Vertex(self, 
                                      x = vor.vertices[i][0], 
                                      y = vor.vertices[i][1])
            
        # create a list of edges
        self.edges = [None] * len(vor.ridge_vertices)
        for i in range(len(vor.ridge_vertices)):
            edge_tiles = [self.all_tiles[j] for j in vor.ridge_points[i]]
            edge_vertices = [self.vertices[j] for j in vor.ridge_vertices[i] if j != -1]
            self.edges[i] = Edge(self, 
                                 tiles = edge_tiles, 
                                 vertices = edge_vertices)
            # update the tiles and vertices with the new edge
            for tile in edge_tiles:
                tile.add_edge(self.edges[i])
            for vertex in edge_vertices:
                vertex.add_edge(self.edges[i])

        logger.info('Voronoi data converted to objects')

    def simulate_tectonics(self, n_plates: int = 10, grow_method = 'random') -> None:
        # randomly select n_plates tiles to be the starting points of the tectonic plates
        starting_tiles = np.random.choice(self.tiles, n_plates, replace=False)
        remaining_tiles = [tile for tile in self.tiles if tile not in starting_tiles]
        self.plates = [TectonicPlate(self, tile) for tile in starting_tiles]
        remaining_plates = self.plates.copy()

        # grow the plates to cover the entire map
        while len(remaining_tiles) > 0 and len(remaining_plates) > 0:
            # randomly select one of the remaining plates
            plate = np.random.choice(remaining_plates)
            new_tile = plate.grow(grow_method=grow_method)
            if new_tile is not None:
                remaining_tiles.remove(new_tile)
            else:
                remaining_plates.remove(plate)
                logger.info(f"Removed {plate}, {len(remaining_plates)} plates and {len(remaining_tiles)} tiles remaining")


    def _build_plot(self, display_method='save', color_by='tile', filename=None, figwidth=6, **kwargs):
        fig = plt.figure(figsize=(figwidth, figwidth * self.height / self.width))
        ax = fig.add_subplot(111)
        ax.set_xlim(0, self.width)
        ax.set_ylim(0, self.height)

        if color_by == 'plate':
            logger.info('Coloring by plate')
            colormap = plt.cm.get_cmap('tab10', len(self.plates))
            for i, plate in enumerate(self.plates):
                plate_color = colormap(i)
                plate.plot(ax, color = plate_color, **kwargs)

        elif color_by == 'tile':
            logger.info('Coloring by tile')
            colormap = plt.cm.get_cmap('cividis', 20)
            for i, tile in enumerate(self.tiles):
                tile_color = colormap(i)
                tile.plot(ax, color = tile_color, **kwargs)

        if display_method == 'save':
            logger.info(f'Saving to {filename}')
            fig.savefig(filename)
        elif display_method == 'display':
            logger.info('Displaying')
            fig.show()

    def display(self, **kwargs):
        self._build_plot(display_method='display', **kwargs)

    def save(self, filename: str, **kwargs):
        self._build_plot(display_method='save', filename=filename, **kwargs)



class Tile():
    def __init__(self, 
                 map: Map,
                 center_x: int, 
                 center_y: int,
                 core: bool) -> None:
        self.map = map
        self.x = center_x
        self.y = center_y
        self.edges = []
        self.plate = None
        self.core = core

    def add_edge(self, edge):
        if edge not in self.edges:
            self.edges.append(edge)
    
    def remove_edge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)

    def set_plate(self, plate):
        self.plate = plate

    def get_neighbors(self):
        neighbors = [tile for edge in self.edges for tile in edge.tiles if tile != self]
        return neighbors
    
    def get_vertices(self):
        vertex_pairs = [edge.vertices for edge in self.edges]
        # start from the first vertex of the first edge
        vertices = [vertex_pairs[0][0]]
        current_vertex = vertex_pairs[0][0]
        # at every step, traverse to the next vertex in the next edge, then delete that edge from vertex_pairs
        while len(vertex_pairs) > 0:
            for i in range(len(vertex_pairs)):
                if current_vertex in vertex_pairs[i]:
                    next_vertex_pair = vertex_pairs.pop(i)
                    next_vertex = [vertex for vertex in next_vertex_pair if vertex != current_vertex][0]
                    vertices.append(next_vertex)
                    current_vertex = next_vertex
                    break
       
        return vertices
    
    def get_vertex_coords(self):
        vertices = self.get_vertices()
        coords = [(vertex.x, vertex.y) for vertex in vertices]
        return coords
    
    def plot(self, ax, color='b'):
        coords = self.get_vertex_coords()
        ax.fill(*zip(*coords), color=color, alpha=0.5)

    def plot_marker(self, ax, color='b', marker='o'):
        ax.plot(self.x, self.y, f'{color}{marker}')

    def __str__(self):
        return f"Tile at ({self.x}, {self.y})"



class Vertex():
    def __init__(self,
                 map: Map,
                 x: int,
                 y: int) -> None:
        self.map = map
        self.x = x
        self.y = y
        self.edges = []

    def add_edge(self, edge):
        if edge not in self.edges:
            self.edges.append(edge)

    def remove_edge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)

    def get_neighbors(self):
        neighbors = [edge for edge in self.edges for vertex in edge.vertices if vertex != self]
        return neighbors
    
    def plot(self, ax, color='r', marker='o'):
        ax.plot(self.x, self.y, f'{color}{marker}')

    def __str__(self):
        return f"Vertex at ({self.x}, {self.y})"
    


class Edge():
    def __init__(self,
                 map: Map,
                 tiles: list[Tile, Tile],
                 vertices: list[Vertex, Vertex]) -> None:
        self.map = map
        self.tiles = tiles
        self.vertices = vertices
        self.plate_boundary = False

    def plot(self, ax, color='k', linestyle='-'):
        coords = [(vertex.x, vertex.y) for vertex in self.vertices]
        ax.plot(*zip(*coords), f'{color}{linestyle}')



class TectonicPlate():
    def __init__(self, 
                 map: Map,
                 starting_tile: Tile) -> None:
        self.map = map
        self.starting_tile = starting_tile
        self.tiles = [starting_tile]
        starting_tile.set_plate(self)
        self.neighbors = starting_tile.get_neighbors()
        self.centroid_x = starting_tile.x
        self.centroid_y = starting_tile.y

        self.movement_direction = np.random.random() * 2 * np.pi
        self.movement_magnitude = np.random.random() # TODO: scaling?

        self.color = None

        self.boundary_tiles = []
        self.boundary_edges = []
        self.boundary_vertices = []

    def grow(self, grow_method = 'random') -> Tile:
        # identify which neighbors are not already part of a plate
        unclaimed_neighbors = [tile for tile in self.neighbors if tile.plate is None and tile.core]
        if len(unclaimed_neighbors) == 0:
            return None
        
        if grow_method == 'random':
            # randomly select one of these neighbors
            new_tile = np.random.choice(unclaimed_neighbors)
        elif grow_method == 'squared_inv_distances':
            # weight the selection by the inverse of the distance to the starting tile
            squared_inv_distances = [1/((tile.x - self.starting_tile.x)**2 + (tile.y - self.starting_tile.y)**2) for tile in unclaimed_neighbors]
            new_tile = np.random.choice(unclaimed_neighbors, p = squared_inv_distances / np.sum(squared_inv_distances))
        elif grow_method == 'squared_inv_centroid_distances':
            # weight the selection by the inverse of the distance to the centroid of the plate
            squared_inv_distances = [1/((tile.x - self.centroid_x)**2 + (tile.y - self.centroid_y)**2) for tile in unclaimed_neighbors]
            new_tile = np.random.choice(unclaimed_neighbors, p = squared_inv_distances / np.sum(squared_inv_distances))

        # add the new tile to the plate
        self.add_tile(new_tile)
        return new_tile

    def add_tile(self, tile) -> None:
        # add the tile and update the tile's plate
        self.tiles.append(tile)
        tile.set_plate(self)
        # update the plate's centroid incrementally
        self.centroid_x = (self.centroid_x * (len(self.tiles) - 1) + tile.x) / len(self.tiles)
        self.centroid_y = (self.centroid_y * (len(self.tiles) - 1) + tile.y) / len(self.tiles)
        # remove the tile from the list of neighbors
        if tile in self.neighbors:
            self.neighbors.remove(tile)
        # add the tile's neighbors to the list of neighbors if not already there
        new_neighbors = [neighbor for neighbor in tile.get_neighbors() if neighbor.plate is not self and neighbor not in self.neighbors and neighbor.core]
        self.neighbors.extend(new_neighbors)

    def plot(self, ax, color=None, show_starting_tile=False):
        logger.info(f"Plotting plate with {len(self.tiles)} tiles in color {color}")
        for tile in self.tiles:
            tile.plot(ax, color=color)
        if show_starting_tile:
            self.starting_tile.plot_marker(ax, color='k', marker='x')

    def _recalculate_centroid(self):
        self.centroid_x = np.mean([tile.x for tile in self.tiles])
        self.centroid_y = np.mean([tile.y for tile in self.tiles])

    def __str__(self):
        return f"Tectonic Plate composed of {len(self.tiles)} tiles starting at {self.starting_tile}"