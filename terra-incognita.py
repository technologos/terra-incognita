#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 20:02:19 2016

@author: CMilroy
"""

import numpy
from numpy.random import random
from matplotlib import pyplot
from scipy.spatial import Voronoi, voronoi_plot_2d, Delaunay

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
  
'''
map classes
'''

class MapTile():
    def __init__(self):
        self.center = [] # coordinates of defining point (~centroid, after relaxation)
        self.vertices = [] # indices of MapVertexes defining the region
        self.vertex_coords = []
        self.neighbors = [] # indices of MapTiles bordering this MapTile
        self.edges = [] #indices of MapEdges bordering this MapTile
        self.biome = []
        self.properties = []

class MapEdge():
    def __init__(self):
        self.tiles = [] # indices of two adjacent MapTiles
        self.vertices = [] # indices of two endpoints
        self.type = []
        self.properties = []

class MapVertex():
    def __init__(self):
        self.edges = [] # indices of connected MapEdges
        

class Map():    
    def __init__(self, 
                 number_of_tiles = 2**8,
                 smoothing_strength = 2):
        
        self.points = random((number_of_tiles, 2))
        
        self.tiles = []
        self.generate_tiles(smoothing_strength)
        
        self.edges = []
        self.generate_adjacencies()
        
        self.vertices
        
    def generate_tiles(self, smoothing_strength = 2):
        self.regularize_tiles(smoothing_strength)
        #self.points = numpy.append(self.points, [[0,0], [0,1], [1,0], [1,1]], axis = 0)
        self.voronoi = Voronoi(self.points)
        for index, point in enumerate(self.voronoi.points):
            new_tile = MapTile()
            new_tile.center = point
            new_tile.vertices = self.voronoi.regions[self.voronoi.point_region[index]]
            if -1 in new_tile.vertices:
                continue
            else:
                new_tile.vertex_coords = numpy.asarray([self.voronoi.vertices[i] 
                                                        for i in new_tile.vertices])
                self.tiles.append(new_tile)
        
    def regularize_tiles(self, number_of_iterations):
        #print "Improving points"
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
        self.delaunay = Delaunay(self.points)
        vnv = self.delaunay.vertex_neighbor_vertices
        #num_simplices = len(self.delaunay.simplices)
        self.neighbors = [vnv[1][vnv[0][i]: vnv[0][i+1]] 
                          for i in range(len(vnv[0])-1)]
        for index, tile in enumerate(self.tiles):
            tile.neighbors = list(self.neighbors[index])
        
    
    def display(self, highlight_tile = -1,
                xlim = [0,1], ylim = [0,1]):
        voronoi_plot_2d(self.voronoi)
        if highlight_tile >= 0 and highlight_tile < len(self.tiles):
            pyplot.fill(*zip(*self.tiles[highlight_tile].vertex_coords), 'y')
                
        pyplot.triplot(self.points[:,0], self.points[:,1], self.delaunay.simplices.copy(), 'r-')
        
        pyplot.xlim(xlim[0], xlim[1])
        pyplot.ylim(ylim[0], ylim[1])
        pyplot.show()


# create tectonic plates
# move plates
# generate tradewinds
# generate base heat
# run simple dynamic weather
# build biomes
# select initial settlement locations
# more...

if __name__ == '__main__':
    pass