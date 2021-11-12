#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:47:40 2021

@author: cosimo
"""
from shapely.geometry import Point
import select_emissions_param as par
def lat_to_index(lat, lon):
    # returns indexes in the netCDF4 matrix relative to input latitude and longitude
    # the round command is used to round lat and long to .05, avoiding .00 (emi_ch4 data are centered at .05).
    # it works like this: 0.05 is added to lat and long, then they are rounded to .1, then 0.05 is subtracted again. 
    # the value 0.05 is subtracted once more in order to obtain a .1 rounded number that can be converted to index multiplying by 10 (the grid points spacing is 0.1°)
    lat_index = int( (round((lat + 0.05)*10)*0.1 - 0.1 + 90)*10 )
    lon_index = int( (round((lon + 0.05)*10)*0.1 - 0.1     )*10 )
    return lat_index, lon_index

def index_to_lat(lat_index, lon_index):
    lat = round( (lat_index - 900) * 0.1 + 0.05, 2)
    lon = round( lon_index * 0.1 + 0.05, 2)
    return lat, lon


def within_internal_boundary(point, boundary): 
    # check if point in inside internal boundary (i.e. at least three out of the 4 angles of the bin are inside the boundary)
    # returns true/false value
    p = [[] for i in range(4)]
    p[0] = Point(point.x +0.5*par.d_lon, point.y +0.5*par.d_lat)
    p[1] = Point(point.x +0.5*par.d_lon, point.y -0.5*par.d_lat)
    p[2] = Point(point.x -0.5*par.d_lon, point.y -0.5*par.d_lat)
    p[3] = Point(point.x -0.5*par.d_lon, point.y +0.5*par.d_lat)
    i=0
    for j in range(4):
        # count number of points inside boundary
        if p[j].within(boundary):
            i=i+1            
    is_within = i>2
    return is_within

def within_external_boundary(point, boundary): 
    # check if point is inside external boundary (i.e. at least one out of the 4 angles of the bin is inside the boundary)
    # returns true/false value    
    p = [[] for i in range(4)]
    p[0] = Point(point.x +0.5*par.d_lon, point.y +0.5*par.d_lat)
    p[1] = Point(point.x +0.5*par.d_lon, point.y -0.5*par.d_lat)
    p[2] = Point(point.x -0.5*par.d_lon, point.y -0.5*par.d_lat)
    p[3] = Point(point.x -0.5*par.d_lon, point.y +0.5*par.d_lat)
    is_within = p[0].within(boundary) or p[1].within(boundary) or p[2].within(boundary) or p[3].within(boundary) or point.within(boundary)
    return is_within

