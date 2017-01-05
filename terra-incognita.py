#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 20:02:19 2016

@author: technologos

Thanks to mewo2 and amitp for inspiration and tutorials.
Thanks particularly to mewo2 for the framework for the regularization.
"""

import numpy
from numpy import sqrt, dot, sign, mean, pi, cos, sin
from numpy.random import random, randint, normal
from numpy.linalg import norm as magnitude
from matplotlib import pyplot
from matplotlib.collections import LineCollection
from matplotlib import colors as mpl_colors
from matplotlib import cm as colormap
from scipy.spatial import Voronoi
from math import degrees, acos
from functools import wraps
from time import time
from random import choice


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

def timed(function):
    
    @wraps(function)
    def timed_function(*args, **kwargs):
        start_time = time()
        function(*args, **kwargs)
        end_time = time()
        print(round(end_time - start_time, 2), ' seconds')
    
    return timed_function
  
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
        self.elevation = 0 # calculated from edges
        self.boundary = False # is this tile on the outisde of a plate?
        self.slope = 0
        self.roughness = 0
        self.water_depth = 0

class MapEdge():
    def __init__(self, tiles):
        self.tiles = tiles # indices of two adjacent MapTiles
        self.vertices = [] # indices of two endpoints
        self.type = []
        self.properties = []
        self.elevation = 0 # derived from adjacent MapTiles' velocity vectors
        self.length = 0
        self.boundary = False

class MapVertex():
    def __init__(self, coords):
        self.coords = coords
        self.edges = [] # indices of connected MapEdges
        self.elevation = -100
        
class TectonicPlate():
    def __init__(self, initial_tile):
        self.plate_type = []
        self.tiles = [initial_tile]
        self.center = initial_tile.center
        self.velocity = numpy.array([0, 0]) # do plates also rotate?
        self.color = tuple(random(3))
        self.boundaries = [] # MapTiles 
        self.elevation = normal(scale = .25) # TODO: refine
        

class Map():    
    def __init__(self, 
                 number_of_tiles = 4096,
                 smoothing_strength = 2):
        
        # define the fundamental topology
        self.points = random((number_of_tiles, 2))
        self.ocean_elevation = random() - .5 # TODO: refine
        
        self.tiles = [] # MapTiles
        print('Generating tiles...')
        self.generate_tiles(smoothing_strength) # center, vertices, vertex_coords
        
        self.edges = [] # MapEdges
        print('Generating adjacencies...')
        self.generate_adjacencies() # neighbors
        
        self.vertices = [] # MapVertexes
        print('Generating vertices...')
        self.generate_vertices()
        
        # create tectonic plates
        self.plates = [] # TectonicPlates
        self.boundaries = [] # MapEdges
        print('Generating plates...')
        self.generate_plates()
        
        # move plates
        print('Moving plates...')
        self.move_plates()
        #     TODO: island-forming volcanoes?
        
        # update elevations
        print('Calculating elevation...')
        self.calculate_elevation()
        self.fill_oceans()
        
        # generate tradewinds
        latitude_span = randint(0,10) + 1
        min_latitude = randint(-60, 61 - latitude_span)
        self.latitudes = (min_latitude, min_latitude + latitude_span)
        self.calculate_tradewinds()
        
        # TODO: generate base heat
        # TODO: run simple dynamic weather
        #     note: stream length and basin size can be checked against Hack's Law
        # TODO: build biomes
        # TODO: select initial settlement locations
        # more... second-gen settlements?
    
    @timed
    def generate_tiles(self, smoothing_strength = 2):
        self.regularize_tiles(smoothing_strength)
        #self.points = numpy.append(self.points, [[0,0], [0,1], [1,0], [1,1]], axis = 0)
        self.voronoi = Voronoi(self.points)
        for index, point in enumerate(self.voronoi.points):
            new_tile = MapTile(point)
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
    
    @timed
    def generate_adjacencies(self):
        for ridge in self.voronoi.ridge_points:
            tiles = [self.tiles[i] for i in ridge]
            new_edge = MapEdge(tiles)
            self.edges.append(new_edge)
            
            for index, tile in enumerate(tiles):
                tiles[index].neighbors.append(tiles[1 - index])
                tiles[index].edges.append(new_edge)
    
    @timed
    def generate_vertices(self):
        for vertex in self.voronoi.vertices:
            new_vertex = MapVertex(vertex)
            self.vertices.append(new_vertex)
        for index, ridge in enumerate(self.voronoi.ridge_vertices):
            for vertex_index in ridge:
                if vertex_index != -1:
                    self.edges[index].vertices.append(self.vertices[vertex_index])
                    self.vertices[vertex_index].edges.append(self.edges[index])
        for index, tile in enumerate(self.tiles):
            vertex_indices = [i for i in self.voronoi.regions[self.voronoi.point_region[index]]
                              if i != -1]
            tile.vertices = [self.vertices[j] for j in vertex_indices]
            tile.vertex_coords = numpy.asarray([self.voronoi.vertices[k] 
                                                for k in vertex_indices])
                             
    @timed
    def generate_plates(self):
        tiles_assigned = 0
        while tiles_assigned < len(self.tiles):
            focus_index = randint(len(self.tiles))
            focus_tile = self.tiles[focus_index]
            if focus_tile.plate == -1:
                new_plate = TectonicPlate(focus_tile)
                
                direction = random() * 2 * pi
                velocity_vector = numpy.array([cos(direction), sin(direction)])
                magnitude = random()
                new_plate.velocity = magnitude * velocity_vector 
                
                focus_tile.plate = new_plate
                tiles_assigned += 1
                self.plates.append(new_plate)
                
            for focus_plate in self.plates:
                tiles_to_add = []
                for tile in focus_plate.tiles:
                    for index, neighbor in enumerate(tile.neighbors):
                        if neighbor.plate == -1: # TODO: consider prob < 1
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
    @timed    
    def move_plates(self):
        for edge in self.boundaries:
            for index, tile in enumerate(edge.tiles):
                normal_vector = edge.tiles[1-index].center - tile.center
                normal_vector /= magnitude(normal_vector) # sets magnitude to 1
                normal_force = dot(tile.plate.velocity, normal_vector)
                edge.elevation +=  (sign(normal_force) * sqrt(abs(normal_force)) + # TODO: sqrt?
                                    tile.plate.elevation) # TODO: adding plate elevation here?
    
    @timed
    def calculate_elevation(self):
        boundary_vertices = [vertex for boundary in self.boundaries 
                             for vertex in boundary.vertices]
        boundary_vertices = list(set(boundary_vertices))
                             
        for vertex in boundary_vertices:
            vertex.elevation = mean([edge.elevation for edge in vertex.edges 
                                     if edge in self.boundaries])
                             
        current_vertices = boundary_vertices
        completed_vertices = len(boundary_vertices)
        total_vertices = len(self.vertices)
        
#        log_map_size = 1 / log2(len(self.points))
        log_map_size = 8 / sqrt(len(self.points)) # TODO: calibrate
        
        while completed_vertices < total_vertices:
            
            new_vertices = [new_vertex for current_vertex in current_vertices
                            for new_edge in current_vertex.edges
                            for new_vertex in new_edge.vertices
                            if new_vertex.elevation == -100]
            new_vertices = list(set(new_vertices))
            
            old_vertices = []
            to_remove = []
            for new_vertex in new_vertices:
                if random() < .5: # TODO: calibrate this
                    plate_elevation = new_vertex.edges[0].tiles[0].plate.elevation
                    new_vertex.elevation = (((.8 + random() * .2) ** log_map_size) * 
                                            ((mean([vertex.elevation # TODO: coefficient
                                                  for edge in new_vertex.edges
                                                  for vertex in edge.vertices
                                                  if vertex.elevation != -100])) - 
                                                   plate_elevation)) + plate_elevation
                else:
                    old_vertex = choice([old_vertex for old_edge in new_vertex.edges
                                         for old_vertex in old_edge.vertices
                                         if old_vertex.elevation != -100])
                    old_vertices.append(old_vertex)
                    to_remove.append(new_vertex)
            
            for vertex in to_remove:
                new_vertices.remove(vertex)

            completed_vertices += len(new_vertices)
            current_vertices = new_vertices
            current_vertices.extend(old_vertices)

        for tile in self.tiles:
            vertex_elevations = [vertex.elevation - self.ocean_elevation for vertex in tile.vertices]
            tile.elevation = mean(vertex_elevations) - self.ocean_elevation # TODO: refine
            tile.slope = max(vertex_elevations) - min(vertex_elevations)
            # TODO: tile.roughness based on coplanarity
    
    def fill_oceans(self):
        for tile in self.tiles:
            tile.water_depth = max(0, -tile.elevation)
            # TODO: only fill oceans if connected to a map edge?
    
    def calculate_tradewinds(self):
#        if abs(self.latitude) <= 5:
         pass  
    
    @timed        
    def display(self, 
                show_grid = True,
                highlight_tile = [-1], 
                show_centers = False, 
                show_intersections = False,
                show_plates = False,
                show_plate_centers = False,
                show_plate_velocities = False,
                show_plate_boundaries = False,
                show_boundary_elevation = False,
                show_tile_elevation = False,
                show_vertex_elevation = False,
                show_tile_elevation_labels = False,
                show_water = False,
                clean = False,
                plate_test = False,
                elevation_test = False,
                xlim = [0.05, .95], 
                ylim = [0.05, .95]):
        
        if clean:
#            show_centers = False
#            show_intersections = False
            xlim = [.05, .95]
            ylim = [.05, .95]
        
        if plate_test:
            show_plates = True
            show_plate_centers = True
            show_plate_velocities = True
            show_plate_boundaries = True
            show_boundary_elevation = True   
            
        if elevation_test:
            show_plate_boundaries = True
            show_boundary_elevation = True
            show_tile_elevation = True

        figure = pyplot.figure()
        axes = figure.gca()
        
        if show_boundary_elevation or show_tile_elevation or show_vertex_elevation:
            color_norm = mpl_colors.Normalize(vmin = -2, vmax = 2 - self.ocean_elevation)
            color_map = pyplot.get_cmap('gist_earth')
            palette = colormap.ScalarMappable(norm = color_norm, cmap = color_map)
            
        
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
        
        if show_plates:
            for plate in self.plates:
                for tile in plate.tiles:
                    pyplot.fill(*zip(*tile.vertex_coords), color = plate.color)
                if show_plate_centers:
                    pyplot.plot(plate.center[0], plate.center[1], 'ko')
                    if show_plate_velocities:
                        pyplot.arrow(plate.center[0], plate.center[1], 
                                     plate.velocity[0], plate.velocity[1],
                                     label = magnitude(plate.velocity))

        if show_plate_boundaries:
            line_segments = []
            colors = []
            for edge in self.boundaries:
                if len(edge.vertices) == 2:
                    line_segments.append([(x, y) for x, y in [vertex.coords 
                                          for vertex in edge.vertices]])
                    if show_boundary_elevation:
                        colors.append(palette.to_rgba(edge.elevation))
                    else:
                        colors.append('k')
            borders = LineCollection(line_segments,
                                     colors=colors,
                                     lw=3.0,
                                     linestyle='solid')
            borders.set_alpha(1.0)
            axes.add_collection(borders)
        
        if show_tile_elevation or show_tile_elevation_labels:
            for tile in self.tiles:
                if show_tile_elevation:
                    if show_water and tile.elevation <= 0:
                        pyplot.fill(*zip(*tile.vertex_coords), 
                                    color = palette.to_rgba(tile.elevation / 10 - .8))
                        # TODO: fix this to use water_depth
                    else:
                        pyplot.fill(*zip(*tile.vertex_coords), 
                                    color = palette.to_rgba(tile.elevation))
                if show_tile_elevation_labels:
                    pyplot.text(tile.center[0], tile.center[1], round(tile.elevation, 2))
        
        if show_vertex_elevation:
            for vertex in self.vertices:
                pyplot.plot(vertex.coords[0], vertex.coords[1], 
                            color = palette.to_rgba(vertex.elevation),
                            marker = 'o')
                
        
        highlight_tile = numpy.asarray(highlight_tile)
            
        if numpy.all(highlight_tile >= 0) and numpy.all(highlight_tile < len(self.tiles)):
            for tile in highlight_tile:
                pyplot.fill(*zip(*self.tiles[tile].vertex_coords), 'y')
        
        pyplot.xlim(xlim[0], xlim[1])
        pyplot.ylim(ylim[0], ylim[1])
        pyplot.show()
        

if __name__ == '__main__':
    pass