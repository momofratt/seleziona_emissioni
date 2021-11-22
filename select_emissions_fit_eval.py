#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:36:55 2021

@author: cosimo
"""
from netCDF4 import Dataset
import sys
import numpy as np 
import math 
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from shapely.geometry import Point
import select_emissions_param as par
import select_emissions_netcdf as net

def eval_emission(spec):
    emission = []
    emission_err = []
    # evaluate emissions for given gas specie
    if spec == 'CH4':  #select year
        start_year = par.start_year_ch4
        end_year = par.end_year_ch4
    elif spec == 'CO':
        start_year = par.start_year_co
        end_year = par.end_year_co
        
    #open total emission file and write header
    outfile_emi = open('./res_files/total_emission_'+par.region+'_'+spec+'_'+str(start_year)+'-'+str(end_year)+'.txt', 'w')
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
        sel_data     = np.zeros((par.nbin,par.nbin)) 
        sel_data_int = np.zeros((par.nbin,par.nbin)) #restricted data selection
        sel_data_ext = np.zeros((par.nbin,par.nbin))
       
        # output file with three columns (lat, lon, emi_ch4)
        outfile = open("./res_files/selected_emissions_"+par.region+"_"+spec+"_"+str(year)+".txt", 'w')
        outfile.write("Lat Lon emi_"+spec+"[kg/m2/s]\n")
        for j in range(0, par.nbin): # loop over data matrix
            for i in range(0,par.nbin):
                lat_ind = par.j_ref - par.nbin_2 + j
                lon_ind = par.i_ref - par.nbin_2 + i
                point_coords = net.index_to_lat( lat_ind, lon_ind)
                point = Point(point_coords[1], point_coords[0])
               
                if point.within(par.regione):
                    sel_data[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                    outfile.write(str(point_coords[0]) + " " + str(point_coords[1]) + " " + str(sel_data[j][i]) + "\n")
                
                if net.within_internal_boundary(point, par.regione):
                    sel_data_int[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                
                if net.within_external_boundary(point, par.regione):
                    sel_data_ext[j][i] = (( data['emi_'+spec.lower()][ lat_ind ][ lon_ind ] ))
                    
        outfile.close()
        
        # ###############################################################################
        ###                      Evaluate total emission                            ###
        ###############################################################################
        tot_emi = 0
        tot_emi_int = 0
        tot_emi_ext = 0
        earth_rad2 = pow(par.earth_rad*1000,2) # square earth radius [m^2]
        d_lat_rad = par.d_lat *math.pi/180 # convert grid point separation lat and lon to radiants
        d_lon_rad = par.d_lon *math.pi/180
        pi = math.pi
       
        for j in range(0, par.nbin):    # loop over latitudes to evaluate total emission
            temp_lat, _ = net.index_to_lat(par.j_ref-par.nbin_2 + j , 1) # get latitude of jth iteration
            bin_surf = earth_rad2 * d_lat_rad * math.sin( 0.5*pi - abs(temp_lat) *pi/180) * d_lon_rad # [m^2] evaluate the surface of the grid point. dS = r*Dtheta * r*sin(theta)*Dphi
            for i in range(1, par.nbin-1):      # loop over longitudes
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
    first_emi_index = start_year - par.start_year_co 
    last_emi_index  = len(emi) - (par.end_year_co - end_year)
    
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
        fig2.suptitle('EDGAR '+spec+' emissions for '+par.region_full+' region and '+str(fit_order)+' order polynomial fit')
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
        fig2.savefig(par.res_folder+'EDGAR_'+par.region+'_'+spec+'_'+str(start_year)+'_'+str(fit_order)+'order_polyfit.pdf', format = 'pdf')
        plt.show()

    return poly, r2


def eval_predicted_values(emi, emi_err, spec):
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
        end_year = par.end_year_co
    elif spec == 'CH4':
        end_year = par.end_year_ch4
        
    predicted_emi_file = open(par.predicted_res_folder+'predicted_'+par.region+ '_' +spec+'_emi.txt', 'w')
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
        fig2.suptitle('ISPRA '+spec+' emissions for '+par.region_full+' region and '+str(fit_order)+' order polynomial fit')
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
        fig2.savefig(par.res_folder+par.region+'_ISPRA_'+spec+'_'+str(start_year)+'_'+str(fit_order)+'order_polyfit.pdf', format = 'pdf')
        plt.show()
    return poly, r2

    
def eval_predicted_values_ISPRA(emi, region, spec):
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

    end_year = 2015        
    predicted_emi_file = open(par.predicted_res_folder+'predicted_'+region+'_ISPRA_'+spec+'_emi.txt', 'w')
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
   
