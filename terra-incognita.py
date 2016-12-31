#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 20:02:19 2016

@author: technologos

Thanks to mewo2 and amitp for inspiration and tutorials.
Thanks particularly to mewo2 for the framework for the regularization.
"""

import numpy
from numpy import sqrt, square, dot
from numpy.random import random, randint
from numpy.linalg import norm as magnitude
from matplotlib import pyplot
from matplotlib.collections import LineCollection
from scipy.spatial import Voronoi
from random import uniform
from math import degrees, acos


'''
helper functions
'''
def centroid(vertices):
    x = vertices[:, 0]
    y = vertices[:, 1]
    x1 = numpy.concatenate((x[1:], x[:1]))
    y1 = numpy.concatenate((y[1:], y[:1]))
    C = (x * y1) - (y * x1)
    A = sum(C / 2)
    C_x = ((x + x1) * C) / (6 * A)
    C_y = ((y + y1) * C) / (6 * A)
    return (sum(C_x), sum(C_y))

def angle(vector1, vector2):
    return degrees(acos(dot(vector1, vector2) / 
                                  (magnitude(vector1) * magnitude(vector2))))
  
'''
map classes
'''

class MapTile():
    def __init__(self, center):
        self.center = center # coordinates of defining point (~centroid, after relaxation)
        self.vertices = [] # indices of MapVertexes defining the region
        self.vertex_coords = []
        self.neighbors = [] # indices of MapTiles bordering this MapTile
        self.edges = [] #indices of MapEdges bordering this MapTile
        self.biome = []
        self.plate = -1 # index of plate on which this tile sits
        self.properties = []
        self.height = 0 # calculated from edges
        self.boundary = False # is this tile on the outisde of a plate?

class MapEdge():
    def __init__(self, tiles):
        self.tiles = tiles # indices of two adjacent MapTiles
        self.vertices = [] # indices of two endpoints
        self.type = []
        self.properties = []
        self.height = 0 # derived from adjacent MapTiles' velocity vectors
        self.length = 0
        self.boundary = False

class MapVertex():
    def __init__(self, coords):
        self.coords = coords
        self.edges = [] # indices of connected MapEdges
        
class TectonicPlate():
    def __init__(self, initial_tile):
        self.plate_type = []
        self.tiles = [initial_tile]
        self.center = initial_tile.center
        self.velocity = numpy.array([0, 0]) # do plates also rotate?
        self.color = tuple(random(3))
        self.boundaries = [] # MapTiles 
        # do plates also need heights?
        

class Map():    
    def __init__(self, 
                 number_of_tiles = 2**8,
                 smoothing_strength = 2):
        
        # define the fundamental topology
        self.points = random((number_of_tiles, 2))
        
        self.tiles = [] # MapTiles
        self.generate_tiles(smoothing_strength) # center, vertices, vertex_coords
        
        self.edges = [] # MapEdges
        self.generate_adjacencies() # neighbors
        
        self.vertices = [] # MapVertexes
        self.generate_vertices()
        
        # create tectonic plates
        self.plates = [] # TectonicPlates
        self.boundaries = [] # MapEdges
        self.generate_plates()
        
        # move plates
        self.move_plates()
        
        # generate tradewinds
        # generate base heat
        # run simple dynamic weather
        # build biomes
        # select initial settlement locations
        # more...
        
    def generate_tiles(self, smoothing_strength = 2):
        self.regularize_tiles(smoothing_strength)
        #self.points = numpy.append(self.points, [[0,0], [0,1], [1,0], [1,1]], axis = 0)
        self.voronoi = Voronoi(self.points)
        for index, point in enumerate(self.voronoi.points):
            new_tile = MapTile(point)
            new_tile.vertices = [i 
                                 for i in self.voronoi.regions[self.voronoi.point_region[index]] 
                                 if i != -1]
            new_tile.vertex_coords = numpy.asarray([self.voronoi.vertices[i] 
                                                    for i in new_tile.vertices])
            self.tiles.append(new_tile)
        
    def regularize_tiles(self, number_of_iterations):
        for _ in range(number_of_iterations):
            vor = Voronoi(self.points)
            new_points = []
            for index in range(len(vor.points)):
                point = vor.points[index,:]
                region = vor.regions[vor.point_region[index]]
                if -1 in region:
                    new_points.append(point)
                else:
                    region_vertices = numpy.asarray([vor.vertices[i,:] for i in region])
                    region_vertices[region_vertices < 0] = 0
                    region_vertices[region_vertices > 1] = 1
                    new_point = centroid(region_vertices)
                    new_points.append(new_point)
            self.points = numpy.asarray(new_points)
            
    def generate_adjacencies(self):
        for ridge in self.voronoi.ridge_points:
            tiles = [self.tiles[i] for i in ridge]
            new_edge = MapEdge(tiles)
            self.edges.append(new_edge)
            
            for index, tile in enumerate(tiles):
                tiles[index].neighbors.append(tiles[1 - index])
                tiles[index].edges.append(new_edge)
#            for point in ridge:
#                neighbor_point = [x for x in ridge if x != point][0]
#                self.tiles[point].neighbors.append(self.tiles[neighbor_point])
#                self.tiles[point].edges.append(new_edge)
    
    def generate_vertices(self):
        for vertex in self.voronoi.vertices:
            new_vertex = MapVertex(vertex)
            self.vertices.append(new_vertex)
        for index, ridge in enumerate(self.voronoi.ridge_vertices):
            for vertex_index in ridge:
                if vertex_index != -1:
                    self.edges[index].vertices.append(self.vertices[vertex_index])
                    self.vertices[vertex_index].edges.append(self.edges[index])
    
    def generate_plates(self):
        tiles_assigned = 0
        while tiles_assigned < len(self.tiles):
            focus_index = randint(len(self.tiles))
            focus_tile = self.tiles[focus_index]
            if focus_tile.plate == -1:
                new_plate = TectonicPlate(focus_tile)
                
                velocity_vector = numpy.array([uniform(-1, 1), uniform(-1, 1)])
                norm_factor = sqrt(numpy.sum(square(velocity_vector))) # makes magnitude 1
                magnitude = random()
                new_plate.velocity = magnitude * velocity_vector / norm_factor 
                
                focus_tile.plate = new_plate
                tiles_assigned += 1
                self.plates.append(new_plate)
                
            for focus_plate in self.plates:
                tiles_to_add = []
                for tile in focus_plate.tiles:
                    for index, neighbor in enumerate(tile.neighbors):
                        if neighbor.plate == -1:
                            neighbor.plate = tile.plate
                            tiles_assigned += 1
                            tiles_to_add.append(neighbor)
                        elif neighbor.plate is not tile.plate:
                            tile.boundary = True
                            focus_plate.boundaries.append(tile)
                focus_plate.tiles.extend(tiles_to_add)
            
        for edge in self.edges:
            if not edge.boundary:
                edge_plates = [tile.plate for tile in edge.tiles]
                if edge_plates[0] is not edge_plates[1]:
                    edge.boundary = True
                    self.boundaries.append(edge)
                    for tile in edge.tiles:
                        tile.boundary = True
        
    def move_plates(self):
        for edge in self.boundaries:
            fault_vector = edge.vertices[0] - edge.vertices[1]
            plate_vectors = [tile.plate.velocity for tile in edge.tiles]
            fault_angles = [angle(fault_vector, plate_vector) for plate_vector in plate_vectors]
            
            if condition:
                edge.type = 'transform'
            elif condition2:
                edge.type = 'convergent'
            else:
                edge.type = 'divergent'
            
            boundary_directions = []
            velocities = []
            for i in range(2):
                boundary_directions.append(edge.tiles[i].center - edge.tiles[1-i].center)
                velocities.append(edge.tiles[i].plate.velocity - edge.tiles[1-i].plate.velocity)
#                similarities.append(cosine_similarity(edge.tiles[i].plate.velocity,
#                                                      edge.tiles[1-i].plate.velocity))
                # if angle between vectors is <20 of parallel, strike-slip
                # if either is converging, -cosine similarity * sum(magnitudes)
                # if either is diverging, cosine similarity * sum(magnitudes)
                # if one is converging and the other is diverging, 
                #     (back vector - fore vector) * CS
                
    
    def display(self, 
                show_grid = True,
                highlight_tile = [-1], 
                show_centers = False, 
                show_intersections = False,
                show_plates = False,
                show_plate_boundaries = False,
                clean = False,
                xlim = [0.05, .95], 
                ylim = [0.05, .95]):
        
        if clean:
            show_centers = False
            show_intersections = False
            xlim = [.05, .95]
            ylim = [.05, .95]
        

        figure = pyplot.figure()
        axes = figure.gca()
        
        if show_grid:
            line_segments = []
            for edge in self.edges:
                if len(edge.vertices) == 2:
                    line_segments.append([(x, y) for x, y in [vertex.coords 
                                          for vertex in edge.vertices]])
            grid = LineCollection(line_segments,
                                     colors='k',
                                     lw=1.0,
                                     linestyle='solid')
            grid.set_alpha(1.0)
            axes.add_collection(grid)

        
        if show_plate_boundaries:
            line_segments = []
            for edge in self.boundaries:
                if len(edge.vertices) == 2:
                    line_segments.append([(x, y) for x, y in [vertex.coords 
                                          for vertex in edge.vertices]])
            borders = LineCollection(line_segments,
                                     colors='k',
                                     lw=3.0,
                                     linestyle='solid')
            borders.set_alpha(1.0)
            axes.add_collection(borders)
            
        
        if show_plates:
            for plate in self.plates:
                for tile in plate.tiles:
                    pyplot.fill(*zip(*tile.vertex_coords), color = plate.color)
                    
        
        highlight_tile = numpy.asarray(highlight_tile)
            
        if numpy.all(highlight_tile >= 0) and numpy.all(highlight_tile < len(self.tiles)):
            for tile in highlight_tile:
                pyplot.fill(*zip(*self.tiles[tile].vertex_coords), 'y')
        
        pyplot.xlim(xlim[0], xlim[1])
        pyplot.ylim(ylim[0], ylim[1])
        pyplot.show()
        

if __name__ == '__main__':
    pass