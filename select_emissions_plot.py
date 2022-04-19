#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:07:11 2021

@author: cosimo
"""
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib.colors as colors
from shapely.geometry import Point
import numpy as np
import select_emissions_param as par
import select_emissions_netcdf as net
import sys


def eval_max_emi_error(emi, emi_err):
    new_error = np.zeros(len(emi_err))
    if len(emi) == len(emi_err):
        for i in range(len(emi_err)):
            new_error[i] = max(emi_err[i], 0.1*emi[i])
    else:
        print('error')
        sys.exit()
    return new_error

def plot_ISPRA_EDGAR(emi, emi_err, spec, anni=0, emi_ispra=0):
    if spec == 'CO':
        yrs = np.arange(par.start_year_co, par.end_year_co+1, 1)
    elif spec == 'CH4':
        yrs = np.arange(par.start_year_ch4, par.end_year_ch4+1, 1)

    emi_err = eval_max_emi_error(emi, emi_err) 
    fig, ax = plt.subplots(1,1, figsize = (9,5))
    fig.suptitle('ISPRA and EDGAR emissions for ' +par.region_full+' region')
    ax.errorbar(yrs,emi,emi_err,fmt='.', elinewidth=1, capsize=3, label = 'EDGAR '+spec)
    if (emi_ispra!=0) & (anni!=0):
        ax.scatter(anni, emi_ispra, c='C1',label = 'ISPRA '+spec)
        ax.set_ylim(0,max(max(emi+emi_err),max(emi_ispra))*1.05)
    else:
        ax.set_ylim(0,max(emi+emi_err)*1.05)
        
    ax.set_xlabel('years')
    ax.set_ylabel(spec+' total emission [t]')
    ax.grid()
    ax.legend()
    fig.savefig(par.res_folder+par.region+'_EDGAR_ISPRA_'+spec+'.pdf', format = 'pdf')
    plt.close(fig)

def plot_EDGAR_monthly(emi, emi_err, emi_monthly, spec):
    if spec == 'CO':
        yrs = np.arange(par.start_year_co, par.end_year_co+1, 1)
        end_year = par.end_year_co
    elif spec == 'CH4':
        yrs = np.arange(par.start_year_ch4, par.end_year_ch4+1, 1)
        end_year = par.end_year_ch4


    predicted_emi_frame = pd.read_csv(par.predicted_res_folder+'predicted_'+par.region+'_'+spec+'_yearly_emi.txt', sep =' ')    
    predicted_emi_frame = predicted_emi_frame[predicted_emi_frame['year']>end_year]

    emi_err_12      = eval_max_emi_error(emi/12, emi_err/12)
    emi_err_pred_12 = eval_max_emi_error((predicted_emi_frame['emi[t]']/12).to_numpy(),(predicted_emi_frame['emi_err[t]']/12).to_numpy())

    fig, ax = plt.subplots(1,1, figsize = (9,5))
    fig.suptitle('EDGAR monthly emissions for ' +par.region_full+' region')
    ax.errorbar(yrs,emi/12,emi_err_12,fmt='.', elinewidth=1, capsize=3, label = 'EDGAR '+spec+ ' monthly average')
    ax.errorbar(predicted_emi_frame['year'],predicted_emi_frame['emi[t]']/12, emi_err_pred_12, fmt='.', elinewidth=1, capsize=3, label = 'EDGAR '+spec+ ' predicted monthly average', c='purple')
    #ax.scatter(anni, emi_ispra, c='C1',label = 'ISPRA '+spec)

    months = np.arange(0,1,1/12)
    anno = 2000
    for year_data in emi_monthly[0:len(emi_monthly)-1]:
        ax.plot(anno+months, year_data, '-o', c='C1', marker='^', ms=3, lw=1)
        anno = anno+1    

    ax.plot(anno+months, emi_monthly[len(emi_monthly)-1], '-o', c='C1', marker='^', ms=3, lw=1, label = 'EDGAR '+spec+' monthly profile')


    ax.set_xlabel('years')
    ax.set_xlim(2000,)
    ax.set_xticks(np.arange(1998, 2022,2))
    ax.set_ylim(bottom=0)
    ax.set_ylabel(spec+' total emission [t]')
    ax.grid()
    ax.legend()
    fig.savefig(par.res_folder+par.region+'_EDGAR_monthly_'+spec+'.pdf', format = 'pdf')
    plt.close(fig)

def plot_predicted(emi, emi_err, spec, emi_ISPRA=0,):
    if spec == 'CO':
        yrs = np.arange(par.start_year_co, par.end_year_co+1, 1)
        start_year = par.start_year_co
        end_year = par.end_year_co
    elif spec == 'CH4':
        yrs = np.arange(par.start_year_ch4, par.end_year_ch4+1, 1)
        start_year = par.start_year_ch4
        end_year = par.end_year_ch4
    emi_err = eval_max_emi_error(emi, emi_err)
    
    predicted_values = pd.read_csv(par.predicted_res_folder+'predicted_'+par.region+'_'+spec+'_emi.txt', sep = ' ')
    predicted_emi_file = open(par.predicted_res_folder+'predicted_'+par.region+'_'+spec+'_yearly_emi.txt', 'w')
    predicted_emi_file.write('year emi[t] emi_err[t]\n')
    for i in range(len(emi)):
        predicted_emi_file.write(str(yrs[i]) +' '+ str(round(emi[i])) +' '+ str(round(emi_err[i])) +'\n')

    fig, ax = plt.subplots(1,1, figsize = (9,5))
    fig.suptitle('Past and predicted emissions for '+spec)
    ax.errorbar(yrs,emi,emi_err,fmt='.', elinewidth=1, capsize=3, label = 'EDGAR')
    if par.IPR: 
        yrs_ISPRA = np.arange(1990,2020,5)
        ax.scatter(yrs_ISPRA, emi_ISPRA, c='C4', label = 'ISPRA')
    for year in range(end_year+1, par.end_year+1):
        sup =  max(np.array(predicted_values[str(year)]) - predicted_values[str(year)].mean())
        inf = -min(np.array(predicted_values[str(year)]) - predicted_values[str(year)].mean())
        last_error = 0
        for i in range(5): #evaluate last error as the mean error of the last 5 years
            last_error = last_error + emi_err[end_year-i-start_year]/5 
        ax.errorbar(year, predicted_values[str(year)].mean(), yerr=max(sup,inf,last_error), c='C1',fmt='.', elinewidth=1, capsize=3)
        if par.IPR: #plot only if IPR=True
            predicted_values_ISPRA = pd.read_csv(par.predicted_res_folder+'predicted_'+par.region+'_ISPRA_'+spec+'_emi.txt', sep = ' ')
            ax.scatter(year, predicted_values_ISPRA[str(year)].mean(), c='C3')

        predicted_emi_file.write(str(year) +' '+ str(round(predicted_values[str(year)].mean())) +' '+ str(round(max(sup,inf,last_error))) +'\n')
    ax.set_xlabel('years')
    if emi_ISPRA!=0:
        ax.set_ylim(0,max(max(emi+emi_err),max(emi_ISPRA))*1.05)
    else:
        ax.set_ylim(0,max(emi+emi_err)*1.05)
        
    ax.set_ylabel(spec+' total emission [t]')
    ax.grid()
    ax.legend(loc='best')
    fig.savefig(par.res_folder+'EDGAR_predicted_'+par.region+'_'+spec+'.pdf', format = 'pdf')
    plt.close(fig)

def plot_2D_emis(spec, year,plot_included_bins): 
    # plot 2d emission and selected area for a given specie in a giver year
    dlat = 0.5 + abs((max(par.max_lat,par.min_lat)-par.lat_ref))
    dlon =   1 + abs((max(par.max_lon,par.min_lon)-par.lon_ref))

    latinf, loninf = net.lat_to_index(par.lat_ref-dlat, par.lon_ref-dlon)
    latsup, lonsup = net.lat_to_index(par.lat_ref+dlat, par.lon_ref+dlon)

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
    ax1 = plt.axes(projection=ccrs.Stereographic(central_latitude=par.lat_CMN, central_longitude=par.lon_CMN))

    # define the coordinate system that the grid lons and grid lats are on
    data_transf = ccrs.PlateCarree()
    plt.pcolormesh(lons, lats, emis, transform=data_transf, norm=colors.LogNorm())

    # Add Colorbar
    plt.colorbar()

    # Add Title
    plt.title(spec+' emission [kg m$^{-2}$ s$^{-1}$] for '+str(year))

    # Add borders
    border_lats = [b[1] for b in par.bound_coords]
    border_lons = [b[0] for b in par.bound_coords]
    ax1.plot(border_lons,border_lats, transform=data_transf,color='m', lw = .5)

    # Add cities
    ax1.plot(par.lon_ref, par.lat_ref,transform=data_transf, marker='x',color='red', label='reference coordinates', markersize=3)
    ax1.plot(par.lon_stat, par.lat_stat,transform=data_transf, marker='^',color='red', label=par.stat_nm, markersize=3.5)

    if plot_included_bins == True:
        for j in range(0, par.nbin): 
                for i in range(0, par.nbin):  
                      lat_ind = par.j_ref - par.nbin_2 + j
                      lon_ind = par.i_ref - par.nbin_2 + i
                      point_coords = net.index_to_lat( lat_ind, lon_ind)
                      point = Point(point_coords[1], point_coords[0])

                      if point.within(par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='o', color='y', markersize=2)

                      if net.within_internal_boundary(point, par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=2)

                      if net.within_external_boundary(point, par.regione):
                          ax1.plot(point_coords[1], point_coords[0], transform=data_transf, marker='x', color='g', markersize=1)


    # Add coastlines and rivers
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.RIVERS)

    gl = ax1.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    plt.savefig(spec+'_emission_'+par.region+'_'+str(year)+'_ER.pdf', format='pdf')
    plt.close('all')
