#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 18:44:27 2016

@author: technologos
"""

    '''
    def clean_voronoi(self):
        # if a vertex is outside the bounding box, replace all references to it with -1
        for vertex, coords in enumerate(self.voronoi.vertices):
            if numpy.any(coords > 1) or numpy.any(coords < 0):
                print('Cleaning out vertex ', vertex)
                for i, region in enumerate(self.voronoi.regions):
                    for j, region_vertex in enumerate(region):
                        if region_vertex == vertex:
                            self.voronoi.regions[i][j] = -1
                            
                for i, ridge_pair in enumerate(self.voronoi.ridge_vertices):
                    for j, ridge_vertex in enumerate(ridge_pair):
                        if ridge_vertex == vertex:
                            self.voronoi.ridge_vertices[i][j] = -1
                 
        # if a ridge is now only between -1 and -1, kill it off
        new_ridge_vertices = []
        new_ridge_points = []
        for i, ridge_pair in enumerate(self.voronoi.ridge_vertices):
            if ridge_pair == [-1, -1]:
                continue
            else:
                new_ridge_vertices.append(ridge_pair)
                new_ridge_points.append(self.voronoi.ridge_points[i])
        self.voronoi.ridge_vertices = new_ridge_vertices
        self.voronoi.ridge_points = new_ridge_points
    '''
        
    '''
    def add_edge_points(self):
        center = self.voronoi.points.mean(axis = 0)
        new_vertices = []
        for points, ridge in zip(self.voronoi.ridge_points, self.voronoi.ridge_vertices):
            ridge = numpy.asarray(ridge)
            if numpy.any(ridge < 0):
                vertex = ridge[ridge >= 0][0] # select the vertex that is not -1
                starting_coords = numpy.asarray(self.voronoi.vertices[vertex])
                tangent = self.voronoi.points[points[1]] - self.voronoi.points[points[0]] # slope of the segment
                tangent /= numpy.linalg.norm(tangent) # puts the vector on the unit circle
                normal = numpy.array([-tangent[1], tangent[0]]) # unit slope of the extended ridge
                midpoint = self.voronoi.points[points].mean(axis = 0) # location of the extended ridge
                direction = numpy.sign(numpy.dot(midpoint - center, normal)) * normal
                # check the relevant intersections with the bounding box
    
                new_vertex = bounding_intersection(starting_coords, direction)
                new_vertices.append(new_vertex)
                
                # append the intersection to the vertex list and update the regions
                self.voronoi.vertices.append(new_vertex)
                vertex_index = len(self.voronoi.vertices)
                
                for point in points:
                    vertex_list = self.voronoi.regions[self.voronoi.point_region[point]]
                    
                
                
        # deal with corners of the bounding box
    '''

def bounding_intersection(point, direction, top_left = [0,0], bottom_right = [1,1]):
    line = numpy.asarray([point, point + direction])
    
    left_bound = numpy.asarray([[0,0], [0,1]])
    right_bound = numpy.asarray([[1,0], [1,1]])
    top_bound = numpy.asarray([[0,1], [1,1]])
    bottom_bound = numpy.asarray([[0,0], [1,0]])
    
    intersections = numpy.asarray((intersection(line, left_bound),
                                  intersection(line, right_bound),
                                  intersection(line, top_bound),
                                  intersection(line, bottom_bound)))
    distances = (intersections - point) / direction
    distances = distances[:,0] # 0 and 1 dimensions should be identical for each
    
    return intersections[distances == min(distances[distances >= 0])]

def intersection(line1, line2):
    x1, y1 = line1[0]
    x2, y2 = line1[1]
    x3, y3 = line2[0]
    x4, y4 = line2[1]
    
    denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if denom != 0:
        A = x1*y2 - y1*x2
        B = x3*y4 - y3*x4
        Px = (A * (x3-x4) - B * (x1-x2)) / denom
        Py = (A * (y3-y4) - B * (y1-y2)) / denom
        return numpy.array((Px, Py))
    else:
        return