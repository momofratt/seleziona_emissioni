#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 11:21:26 2021

@author: Cosimo Fratticioli
@contact: c.fratticioli@isac.cnr.it
"""
from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.path import Path
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import itertools
from shapely.geometry import Point
import shapely
from os import sys
import pandas as pd
import timeit 
from sklearn.metrics import r2_score
import read_ISPRA

###############################################################################
###                             FUNCTIONS                                   ###
###############################################################################
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
    p[0] = Point(point.x +0.5*d_lon, point.y +0.5*d_lat)
    p[1] = Point(point.x +0.5*d_lon, point.y -0.5*d_lat)
    p[2] = Point(point.x -0.5*d_lon, point.y -0.5*d_lat)
    p[3] = Point(point.x -0.5*d_lon, point.y +0.5*d_lat)
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
    p[0] = Point(point.x +0.5*d_lon, point.y +0.5*d_lat)
    p[1] = Point(point.x +0.5*d_lon, point.y -0.5*d_lat)
    p[2] = Point(point.x -0.5*d_lon, point.y -0.5*d_lat)
    p[3] = Point(point.x -0.5*d_lon, point.y +0.5*d_lat)
    is_within = p[0].within(boundary) or p[1].within(boundary) or p[2].within(boundary) or p[3].within(boundary) or point.within(boundary)
    return is_within

def eval_emission(spec):
    emission = []
    emission_err = []
    # evaluate emissions for given gas specie
    if spec == 'CH4':  #select year
        start_year = start_year_ch4
        end_year = end_year_ch4
    elif spec == 'CO':
        start_year = start_year_co
        end_year = end_year_co
        
    #open total emission file and write header
    outfile_emi = open('./res_files/total_emission_'+region+'_'+spec+'_'+str(start_year)+'-'+str(end_year)+'.txt', 'w')
    outfile_emi.write('Year tot_emi[kg] tot_emi_error[kg]\n')
    
    for year in range(start_year, end_year+1):
        if spec == 'CH4':
            file_path = './data_CH4/v6.0_CH4_'+str(year)+'_TOTALS.0.1x0.1/'
            file_name = 'v6.0_CH4_'+str(year)+'_TOTALS.0.1x0.1.nc'
        elif spec == 'CO':
            file_path = './data_CO/v50_CO_'+str(year)+'.0.1x0.1/'
            file_name = 'v50_CO_'+str(year)+'.0.1x0.1.nc'
        else:
            print('undefined gas specie')
            sys.exit()
        
        data = Dataset(file_path+file_name, 'r')
        
        # define matrices to store selected data
        sel_data     = np.zeros((nbin,nbin)) 
        sel_data_int = np.zeros((nbin,nbin)) #restricted data selection
        sel_data_ext = np.zeros((nbin,nbin))
       
        # output file with three columns (lat, lon, emi_ch4)
        outfile = open("./res_files/selected_emissions_"+region+"_"+spec+"_"+str(year)+".txt", 'w')
        outfile.write("Lat Lon emi_"+spec+"[kg/m2/s]\n")
        for j in range(0, nbin): # loop over data matrix
            for i in range(0,nbin):
                lat_ind = j_ref - nbin_2 + j
                lon_ind = i_ref - nbin_2 + i
                point_coords = index_to_lat( lat_ind, lon_ind)
                point = Point(point_coords[1], point_coords[0])
               
                if point.within(regione):
                    sel_data[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                    outfile.write(str(point_coords[0]) + " " + str(point_coords[1]) + " " + str(sel_data[j][i]) + "\n")
                
                if within_internal_boundary(point, regione):
                    sel_data_int[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                
                if within_external_boundary(point, regione):
                    sel_data_ext[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                    
        outfile.close()
        
        # ###############################################################################
        ###                      Evaluate total emission                            ###
        ###############################################################################
        tot_emi = 0
        tot_emi_int = 0
        tot_emi_ext = 0
        earth_rad2 = pow(earth_rad*1000,2) # square earth radius [m^2]
        d_lat_rad = d_lat *math.pi/180 # convert grid point separation lat and lon to radiants
        d_lon_rad = d_lon *math.pi/180
        pi = math.pi
       
        for j in range(0, nbin):    # loop over latitudes to evaluate total emission
            temp_lat, _ = index_to_lat(j_ref-nbin_2 + j , 1) # get latitude of jth iteration
            bin_surf = earth_rad2 * d_lat_rad * math.sin( 0.5*pi - abs(temp_lat) *pi/180) * d_lon_rad # [m^2] evaluate the surface of the grid point. dS = r*Dtheta * r*sin(theta)*Dphi
            for i in range(1, nbin-1):      # loop over longitudes
                tot_emi     = tot_emi     + sel_data[j][i]*bin_surf  # add the grid point emission to total emission
                tot_emi_int = tot_emi_int + sel_data_int[j][i]*bin_surf
                tot_emi_ext = tot_emi_ext + sel_data_ext[j][i]*bin_surf
        
        outfile_emi.write(str(year) + " " + str(round(tot_emi*60*60*24*365)) + " " + str( round(max(abs(tot_emi-tot_emi_int), abs(tot_emi-tot_emi_ext))*60*60*24*365) ) +'\n')
        emission.append(tot_emi*60*60*24*365)
        emission_err.append( max(abs(tot_emi-tot_emi_int), abs(tot_emi-tot_emi_ext))*60*60*24*365 )
    
    outfile_emi.close() # close total emission file
    return np.array(emission)/1000, np.array(emission_err)/1000 # convert to Mg


def fit_emis(spec, emi, emi_err, start_year, end_year, fit_order, plot):
    # polynomial fit on EDGAR data for given specie on selected period
  
    # In the case where start_year and end_year are different from the detault ones,  
    # reshape the emi vector to a shorter one in order to fit with the yrs array length
    first_emi_index = start_year - start_year_co 
    last_emi_index  = len(emi) - (end_year_co - end_year)
    
    emi = emi[first_emi_index:last_emi_index]
    emi_err = emi_err[first_emi_index:last_emi_index]
    
    yrs = np.arange(0, end_year-start_year+1, 1)
    fit_res = np.polyfit(yrs, emi, fit_order)
    poly = np.poly1d(fit_res)
    xyears = np.arange(-3, end_year-start_year+6, 1)
    r2 = r2_score(emi, poly(yrs)) # eval r-square
    
    ### Plot
    if plot == True:
        fig2, ax2 = plt.subplots(1,1, figsize = (9,4))
        fig2.suptitle('EDGAR '+spec+' emissions for '+region_full+' region and '+str(fit_order)+' order polynomial fit')
        ax2.grid(True, zorder=0)
        ax2.errorbar(yrs+start_year,emi,emi_err,fmt='.', elinewidth=1, capsize=3)
        ax2.plot(xyears+start_year, poly(xyears), c='C1')
            
        #### define string with fit results as a function of the fit order  ####
        f_string = 'f(x) = '   
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        for i in range(len(fit_res)-1):
            f_string = f_string + 'a$_'+str(i)+'$x$^{'+str(len(fit_res)-1-i)+'}$+'    
        f_string = f_string + 'a$_'+str(i+1)+'$\n'    
        for i in range(len(fit_res)):
            f_string = f_string + ' a$_'+str(i)+'$='+str(round(fit_res[i],2))+'\n'
        f_string = f_string + '$r^2$ = ' + str(round(r2,3))
            
        ax2.text(0.03, 0.05, f_string, transform=ax2.transAxes, bbox=props) 
        ########################################################################
        #### plot predicted points and their values
        missing_points_str = ''
        for i in range(1,2020-end_year+1):
            ax2.scatter(end_year+i, poly(end_year-start_year+i), c='r')
            missing_points_str = missing_points_str + str(end_year+i) + ' = ' + str(round(poly(end_year-start_year+i))) +' t'
            if i < 2020-end_year:
                missing_points_str=missing_points_str+'\n'
        ax2.text(0.8, 0.7, missing_points_str, transform=ax2.transAxes, bbox=props)
        ########################################################################
        
        ax2.set_ylim(0,max(emi+emi_err)*1.05)
        ax2.set_xlabel('years')
        ax2.set_ylabel(spec+' total emission [t]')
        fig2.savefig(res_folder+'EDGAR_'+region+'_'+spec+'_'+str(start_year)+'_'+str(fit_order)+'order_polyfit.pdf', format = 'pdf')
        plt.show()

    return poly, r2


def eval_predicted_values(spec):
    """ 
    Evaluates the predicted values according to past values. The predicted values are obtained by multiple polynomial fits performed on the available data.
    A polinomial fit is performed for different orders of the polynomial and different range of the x axis.
    The predicted values are evaluated selecting only values which time variation is not much different from the variation of the measured ones.
    In particular, the variation respect to the last 5 years of measured data is evaluated and only predicted values which variation 
    respect to the last measured value is less than the previous variation are taken into account.

    Parameters
    spec : str specie('CH4' or 'CO') to be evaluated

    Returns
    None.

    """
    if spec == 'CO':
        emi = emi_CO
        emi_err = emi_err_CO
        end_year = end_year_co
    elif spec == 'CH4':
        emi = emi_CH4
        emi_err = emi_err_CH4
        end_year = end_year_ch4
    predicted_emi_file = open('predicted_'+region+ '_' +spec+'_emi.txt', 'w')
    predicted_emi_file.write('start_year fit_order r2')
    for year in range(end_year, 2021): # add one col for each predicted year
        predicted_emi_file.write(' '+str(year))
    predicted_emi_file.write('\n')
    for start_year in range(1990, end_year-5):
        for fit_order in range(1,5):
            poly, r2 = fit_emis(spec, emi, emi_err, start_year, end_year, fit_order, False)
            
            variation_2020 = abs(  poly(end_year-start_year+5) - emi[end_year-start_year]) # variation of the predicted value respect to the last measurement 
            variation_2010 = abs(emi[end_year-start_year-5] - emi[end_year-start_year]) # variation of the last 5 years of measured data
            
            if (variation_2020 < variation_2010) and (poly(end_year-start_year+5) < poly(end_year-start_year-5)):            
                predicted_emi_file.write(str(start_year)+' '+str(fit_order) + ' ' + str(round(r2,3))+ ' ')
                for year in range(end_year, 2020):
                    predicted_emi_file.write(str(round(poly(year-start_year))) + ' ')
                predicted_emi_file.write(str(round(poly(2020-start_year)))+'\n')
    predicted_emi_file.close()
    


def fit_emis_ISPRA(spec, emi, start_year, end_year, fit_order, plot):
    # polynomial fit on ISPRA data for given specie on selected period
  
    # In the case where start_year and end_year are different from the default ones,  
    # reshape the emi vector to a shorter one in order to fit with the yrs array length
    first_emi_index = int((start_year - 1990)/5)
    last_emi_index  = len(emi)
    
    emi = emi[first_emi_index:last_emi_index]
    
    yrs = np.arange(0, end_year-start_year+1, 5)
    fit_res = np.polyfit(yrs, emi, fit_order)
    poly = np.poly1d(fit_res)
    xyears = np.arange(-5, end_year-start_year+10, 1)
    r2 = r2_score(emi, poly(yrs)) # eval r-square
    
    ### Plot
    if plot == True:
        fig2, ax2 = plt.subplots(1,1, figsize = (9,4))
        fig2.suptitle('ISPRA '+spec+' emissions for '+region_full+' region and '+str(fit_order)+' order polynomial fit')
        ax2.grid(True, zorder=0)
        ax2.scatter(yrs+start_year,emi)
        ax2.plot(xyears+start_year, poly(xyears), c='C1')
            
        #### define string with fit results as a function of the fit order  ####
        f_string = 'f(x) = '   
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
        for i in range(len(fit_res)-1):
            f_string = f_string + 'a$_'+str(i)+'$x$^{'+str(len(fit_res)-1-i)+'}$+'    
        f_string = f_string + 'a$_'+str(i+1)+'$\n'    
        for i in range(len(fit_res)):
            f_string = f_string + ' a$_'+str(i)+'$='+str(round(fit_res[i],2))+'\n'
        f_string = f_string + '$r^2$ = ' + str(round(r2,3))
            
        ax2.text(0.03, 0.05, f_string, transform=ax2.transAxes, bbox=props) 
        ########################################################################
        #### plot predicted points and their values
        missing_points_str = ''
        for i in range(1,2020-end_year+1):
            ax2.scatter(end_year+i, poly(end_year-start_year+i), c='r')
            missing_points_str = missing_points_str + str(end_year+i) + ' = ' + str(round(poly(end_year-start_year+i))) +' t'
            if i < 2020-end_year:
                missing_points_str=missing_points_str+'\n'
        ax2.text(0.8, 0.7, missing_points_str, transform=ax2.transAxes, bbox=props)
        ########################################################################
        
        ax2.set_ylim(0,max(emi)*1.15)
        ax2.set_xlabel('years')
        ax2.set_ylabel(spec+' total emission [t]')
        fig2.savefig(res_folder+region+'_ISPRA_'+spec+'_'+str(start_year)+'_'+str(fit_order)+'order_polyfit.pdf', format = 'pdf')
        plt.show()
    return poly, r2

    
def eval_predicted_values_ISPRA(spec):
    """ 
    Evaluates the predicted values according to past values. The predicted values are obtained by multiple polynomial fits performed on the available data.
    A polinomial fit is performed for different orders of the polynomial and different range of the x axis.
    The predicted values are evaluated selecting only values which time variation is not much different from the variation of the measured ones.
    In particular, the variation respect to the last 5 years of measured data is evaluated and only predicted values which variation 
    respect to the last measured value is less than the previous variation are taken into account.

    Parameters
    spec : str specie('CH4' or 'CO') to be evaluated

    Returns
    None.

    """
    if spec == 'CO':
        emi = emi_CO_ISPRA
    elif spec == 'CH4':
        emi = emi_CH4_ISPRA
        
    end_year = 2015        
    predicted_emi_file = open('predicted_'+region+'_ISPRA_'+spec+'_emi.txt', 'w')
    predicted_emi_file.write('start_year fit_order r2')
    for year in range(end_year, 2021): # add one col for each predicted year
        predicted_emi_file.write(' '+str(year))
    predicted_emi_file.write('\n')
    for start_year in range(1990, 2005, 5):
        for fit_order in range(1,5):
            poly, r2 = fit_emis_ISPRA(spec, emi, start_year, end_year, fit_order, True)
            
            variation_2020 = abs(poly(end_year-start_year+5) - emi[int((end_year-start_year)/5)]) # variation of the predicted value respect to the last measurement 
            variation_2005 = abs( emi[int((end_year-start_year-10)/5)] - emi[int((end_year-start_year)/5)]) # variation of the last 5 years of measured data
            
            if (variation_2020 < variation_2005) and (poly(end_year-start_year+5) < poly(end_year-start_year-5)):            
                predicted_emi_file.write(str(start_year)+' '+str(fit_order) + ' ' + str(round(r2,3))+ ' ')
                for year in range(end_year, 2020):
                    predicted_emi_file.write(str(round(poly(year-start_year))) + ' ')
                predicted_emi_file.write(str(round(poly(2020-start_year)))+'\n')
    predicted_emi_file.close()
   

         
def plot_ISPRA_EDGAR(spec):
    if spec == 'CO':
        yrs = np.arange(start_year_co, end_year_co+1, 1)
        emi = emi_CO
        emi_err = emi_err_CO
        emi_ispra = emi_CO_ISPRA
    elif spec == 'CH4':
        yrs = np.arange(start_year_ch4, end_year_ch4+1, 1)
        emi = emi_CH4
        emi_err = emi_err_CH4
        emi_ispra = emi_CH4_ISPRA
        
        
    fig, ax = plt.subplots(1,1, figsize = (9,5))
    fig.suptitle('ISPRA and EDGAR emissions for ' +region_full+' region')
    ax.errorbar(yrs,emi,emi_err,fmt='.', elinewidth=1, capsize=3, label = 'EDGAR '+spec)
    ax.scatter(anni, emi_ispra, c='C1',label = 'ISPRA '+spec)
    ax.set_xlabel('years')
    ax.set_ylim(0,max(emi+emi_err)*1.05)
    ax.set_ylabel(spec+' total emission [t]')
    ax.grid()
    ax.legend()
    fig.savefig(res_folder+region+'_EDGAR_ISPRA_'+spec+'.pdf', format = 'pdf')
    
def plot_predicted(spec):
    if spec == 'CO':
        yrs = np.arange(start_year_co, end_year_co+1, 1)
        emi = emi_CO
        if IPR:
            emi_ISPRA = emi_CO_ISPRA
        emi_err = emi_err_CO
        start_year = start_year_co
        end_year = end_year_co
    elif spec == 'CH4':
        yrs = np.arange(start_year_ch4, end_year_ch4+1, 1)
        emi = emi_CH4
        if IPR:
            emi_ISPRA = emi_CH4_ISPRA
        emi_err = emi_err_CH4
        start_year = start_year_ch4
        end_year = end_year_ch4
    
    yrs_ISPRA = np.arange(1990,2020,5)
    predicted_values = pd.read_csv('predicted_'+region+'_'+spec+'_emi.txt', sep = ' ')
    predicted_values_ISPRA = pd.read_csv('predicted_ISPRA_'+spec+'_emi.txt', sep = ' ')
    predicted_emi_file = open('predicted_'+region+'_'+spec+'_yearly_emi.txt', 'w')
    predicted_emi_file.write('year emi[t] emi_err[t]\n')
    for i in range(len(emi)):
        predicted_emi_file.write(str(yrs[i]) +' '+ str(round(emi[i])) +' '+ str(round(emi_err[i])) +'\n')
    
    fig, ax = plt.subplots(1,1, figsize = (9,5))
    fig.suptitle('Past and predicted emissions for '+spec)
    ax.errorbar(yrs,emi,emi_err,fmt='.', elinewidth=1, capsize=3, label = 'EDGAR')
    if IPR: 
        ax.scatter(yrs_ISPRA, emi_ISPRA, c='C4', label = 'ISPRA')
    for year in range(end_year+1, 2021):
        sup =  max(np.array(predicted_values[str(year)]) - predicted_values[str(year)].mean())
        inf = -min(np.array(predicted_values[str(year)]) - predicted_values[str(year)].mean())
        last_error = 0
        for i in range(5): #evaluate last error as the mean error of the last 5 years
            last_error = last_error + emi_err[end_year-i-start_year]/5 
        ax.errorbar(year, predicted_values[str(year)].mean(), yerr=max(sup,inf,last_error), c='C1',fmt='.', elinewidth=1, capsize=3)
        if IPR: #plot only if IPR=True
            ax.scatter(year, predicted_values_ISPRA[str(year)].mean(), c='C3')
        
        predicted_emi_file.write(str(year) +' '+ str(round(predicted_values[str(year)].mean())) +' '+ str(round(max(sup,inf,last_error))) +'\n')
    ax.set_xlabel('years')
    ax.set_ylim(0,max(emi+emi_err)*1.05)
    ax.set_ylabel(spec+' total emission [t]')
    ax.grid()
    ax.legend(loc='best')
    fig.savefig(res_folder+'EDGAR_predicted_'+region+'_'+spec+'.pdf', format = 'pdf')

def plot_2D_emis(spec, year): 
    # plot 2d emission and selected area for a given specie in a giver year
    dlat = 0.5 + abs((max(max_lat,min_lat)-lat_ref))
    dlon =   1 + abs((max(max_lon,min_lon)-lon_ref))
    latinf, loninf = lat_to_index(lat_ref-dlat, lon_ref-dlon)
    latsup, lonsup = lat_to_index(lat_ref+dlat, lon_ref+dlon)
    
    if spec == 'CH4':
        file_path = './data_CH4/v6.0_CH4_'+str(year)+'_TOTALS.0.1x0.1/'
        file_name = 'v6.0_CH4_'+str(year)+'_TOTALS.0.1x0.1.nc'
    elif spec == 'CO':
        file_path = './data_CO/v50_CO_'+str(year)+'.0.1x0.1/'
        file_name = 'v50_CO_'+str(year)+'.0.1x0.1.nc'
    else:
        print('undefined gas specie')
        sys.exit()
        
    data = Dataset(file_path+file_name, 'r')
        
    ### Plot selected region ###
    lats = data['lat'][latinf:latsup]
    lons = data['lon'][loninf:lonsup]
    emis = data['emi_'+spec.lower()][latinf:latsup, loninf:lonsup] 
    
    # set up a map
    # ATTENZIONE: usi plt.axes-> non plt.subplots!!!
    # occhio che succede comunque del casino con i grafici precedenti!!!
    ax1 = plt.axes(projection=ccrs.Stereographic(central_latitude=lat_CMN, central_longitude=lon_CMN))
    
    # define the coordinate system that the grid lons and grid lats are on
    data_transf = ccrs.PlateCarree()
    plt.pcolormesh(lons, lats, emis, transform=data_transf, norm=colors.LogNorm())
    
    # Add Colorbar
    plt.colorbar()
    
    # Add Title
    plt.title(spec+' emission [kg m$^{-2}$ s$^{-1}$] for '+str(year))
    
    # Add borders
    border_lats = [b[1] for b in bound_coords]
    border_lons = [b[0] for b in bound_coords]
    ax1.plot(border_lons,border_lats, transform=data_transf,color='m', lw = .5)
    
    # Add cities
    ax1.plot(lon_ref, lat_ref,transform=data_transf, marker='^',color='red', label='reference coordinates', markersize=3)
    #ax.plot(lon_bolo, lat_bolo,transform=data_transf, marker='x',color='orange', label='Bologna', markersize=3)
    
    for j in range(0, nbin): 
            for i in range(0, nbin):  
                  lat_ind = j_ref - nbin_2 + j
                  lon_ind = i_ref - nbin_2 + i
                  point_coords = index_to_lat( lat_ind, lon_ind)
                  point = Point(point_coords[1], point_coords[0])
                
                  if point.within(regione):
                      ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='o', color='y', markersize=2)
                  
                  if within_internal_boundary(point, regione):
                      ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=2)
                  
                  if within_external_boundary(point, regione):
                      ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=1)
                   
                      
    # Add coastlines and rivers
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.RIVERS)
    
    gl = ax1.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    plt.savefig(spec+'_emission_'+region+'_'+str(year)+'_ER.pdf', format='pdf')
    plt.show()


###############################################################################
###                         Read file and parameters                        ###
###############################################################################

#################################################
################ USER PARAMETERS ################
res_folder = './results/'

# seleziona province. NB: le province selezionate DEVONO essere adiacenti. e.g. [Ferrara, Bologna, Modena] è una buona selezione, [Ferrara, Bologna, Parma] invece NO

# ############## EMILIA ROMAGNA ###############
region = 'ER' # suffix for ouput files
region_full = 'Emilia Romagna'
provinces = ['Modena', 'Reggio Emilia', 'Parma', 'Piacenza', 'Bologna', 'Ferrara', 'Ravenna', 'Forlì-Cesena', 'Rimini']
IPR=True
###############################################
# ############## TOSCANA ###############
# region = 'TOS' # suffix for ouput files
# region_full = 'Toscana'
# provinces = ['Grosseto', 'Livorno', 'Pisa', 'Massa-Carrara', 'Lucca', 'Pistoia', 'Prato', 'Firenze', 'Arezzo', 'Siena']
# IPR=True
########################################
start_year_ch4 = 1990
end_year_ch4   = 2018
start_year_co  = 1990
end_year_co    = 2015
lon_CMN, lat_CMN   = 10.70111, 44.19433 # CMN coordinates 
lon_bolo, lat_bolo = 11.3435 , 44.4947  # Bologna
#################################################
#################################################

#################################################
################ FIXED PARAMETERS ###############
########### DO NOT MODIFY THIS SECTION ##########
d_lon = 0.1 # separation of grid points in deg
d_lat = 0.1
earth_rad = 6371 # [km]
j_CMN, i_CMN = lat_to_index(lat_CMN, lon_CMN)
lat_CMN, lon_CMN = index_to_lat(j_CMN, i_CMN) # obtain approximated lat_CMN, lon_CMN cendered on grid points 
#################################################
#################################################

###############################################################################
###                         Select emission region                          ###
###############################################################################

states_geoframe = gpd.GeoDataFrame.from_file('admin1.geojson')        # read geojson file

prov_geoframe = gpd.GeoDataFrame()  # select provinces

for prov in provinces:
    new_prov = states_geoframe[states_geoframe['name']==prov] # read line from the GeoDataFrame that correspond to the selected province
    if type(new_prov['geometry'].to_list()[0]) == shapely.geometry.multipolygon.MultiPolygon: # if prov is a multipolygon (e.g. islands) uses only the largest polygon 
        prov_poly = new_prov['geometry'].to_list()    
        max_prov_poly = max(prov_poly[0], key=lambda a: a.area)
        new_prov['geometry'] = max_prov_poly
    
    prov_geoframe = prov_geoframe.append(new_prov, ignore_index=True)

# union of provinces and evaluation of external boundaries
geoms = prov_geoframe['geometry'].tolist()
intersection_iter = gpd.GeoDataFrame(gpd.GeoSeries([poly[0].union(poly[1]) for poly in  itertools.combinations(geoms, 2)]), columns=['geometry'])
regione = intersection_iter.unary_union
bound_coords = list(regione.exterior.coords)

bound_coords_index = [lat_to_index(coord[1], coord[0]) for coord in bound_coords] # convert boundary coordinates to indexes
bound_coords_index = list(dict.fromkeys(bound_coords_index))                      # remove duplicate elements in the list

max_lat = max(b[1] for b in bound_coords)
min_lat = min(b[1] for b in bound_coords)
max_lon = max(b[0] for b in bound_coords)
min_lon = min(b[0] for b in bound_coords)

max_coords_ind = lat_to_index(max_lat, max_lon)
min_coords_ind = lat_to_index(min_lat, min_lon)

nbin = max(abs(max_coords_ind[0] - min_coords_ind[0]), abs(max_coords_ind[1] - min_coords_ind[1])) #maximum extension in number of bins between N-S and E-W directions
nbin = nbin+10 # add bins to avoid problems at the borders

if nbin % 2 != 0:
    nbin=nbin+1    
nbin_2=int(0.5*nbin) # half width of nbin

j_ref, i_ref = lat_to_index(0.5*(max_lat -min_lat)+min_lat, 0.5*(max_lon-min_lon)+min_lon)  # reference coordinates to center the sel_ data matrix
lat_ref, lon_ref = index_to_lat(j_ref, i_ref)

###############################################################################
###   Select emissions inside emission region and write output on files     ###
###############################################################################
# define arrays to store selected data with respective latitudes and longitudes
# evaluate arrays for emission and relative error 
emi_CH4, emi_err_CH4 = eval_emission('CH4')
emi_CO, emi_err_CO = eval_emission('CO')

if IPR:
    # get ISPRA emissions and write on file
    emi_CO_ISPRA, emi_CH4_ISPRA, anni = read_ISPRA.select_ISPRA(region)
    
    file_emi_ch4_ISPRA = open('./res_files/total_emissions_'+region+'_CH4_ISPRA.txt', 'w')
    file_emi_ch4_ISPRA.write('year tot_emi[t]\n')
    file_emi_co_ISPRA = open('./res_files/total_emissions_'+region+'_CO_ISPRA.txt', 'w')
    file_emi_co_ISPRA.write('Year emi_CO[t]\n')
    for i in range(len(emi_CH4_ISPRA)):
        file_emi_co_ISPRA.write(str(anni[i]) +' '+ str(round(emi_CO_ISPRA[i])) +'\n')
        file_emi_ch4_ISPRA.write(str(anni[i]) +' '+ str(round(emi_CH4_ISPRA[i])) +'\n')
    
    file_emi_co_ISPRA.close()
    file_emi_ch4_ISPRA.close()

###############################################################################
###                               FIT                                       ###
###############################################################################

fit_emis('CO', emi_CO, emi_err_CO, start_year_co, end_year_co, 2, True)
fit_emis('CH4', emi_CH4, emi_err_CH4, start_year_ch4, end_year_ch4, 4, True)
fit_emis('CO', emi_CO, emi_err_CO, 2000, 2015, 3, True)


eval_predicted_values('CO')
eval_predicted_values('CH4')

if IPR:
    fit_emis_ISPRA('CO', emi_CO_ISPRA,  1990, 2015, 2, True)
    eval_predicted_values_ISPRA('CO')
    eval_predicted_values_ISPRA('CH4')


###############################################################################
###                               PLOT                                      ###
###############################################################################
if IPR:
    plot_ISPRA_EDGAR('CO')
    plot_ISPRA_EDGAR('CH4')    

plot_predicted('CO')
plot_predicted('CH4')
plot_2D_emis('CH4', 2015)

################################################################
####                    TO BE CORRECTED                   ######
################################################################

# output file with a lat-long matrix which elements are the emission values of CH4 at the respective coordinates
# outfile2D = open("selected_emissions_2D_ER.txt", 'w')

# for i in range(0, nbin-1): # write header (longitudes)
#     _, temp_long = index_to_lat(j_ref-nbin_2 , i_ref-nbin_2+i)
#     outfile2D.write(str(temp_long) + " ")

# i=i+1 # write last longitude and "\n"
# _, temp_long = index_to_lat(j_ref-nbin_2 , i_ref-nbin_2+i)
# outfile2D.write(str(temp_long) + "\n")

# for j in range(0, nbin):       # loop over latitudes 
#     i=0
#     temp_lat, _ = index_to_lat(j_ref-nbin_2 + j , i_ref-nbin_2)
#     outfile2D.write(str(temp_lat) + " ")
#     for i in range(0, nbin-1):        # loop over longitudes
#         outfile2D.write(str(sel_data[j][i]) + " ")
#     i=nbin-1
#     outfile2D.write(str(sel_data[j][i]) + "\n")
      
# outfile2D.close()