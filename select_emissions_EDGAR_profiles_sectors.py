#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:22:41 2021

@author: cosimo
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import select_emissions_param as par

def get_region_df(profile_region):
    """
    Parameters
    ----------
    profile_region : int
        yearly profile mapping number (see Crippa et al. 2019)
    Returns
    -------
    profile_frame : DataFrame
    """
    profile_frame = pd.read_csv("./sector_monthly_gridmaps/EDGAR_temporal_profiles_r1.csv", sep = ' ')
    profile_frame = profile_frame[profile_frame['Region/country']==str(profile_region)]
    return profile_frame

def get_sector_profile(year, sector, profile_frame):
    """
    get the profile for a given emission sector during a selected year

    Parameters
    ----------
    year : int
        selected year for the analysis.
    sector : str
        emission sector according to IPCC 2006 (see Crippa et al. 2019).
    profile_frame: DataFrame
        frame to perform the computation of the yearly profile.
    Returns
    -------
    out_array : array
        array containing the monthly coefficients to compute yearly variation
    """    
    no_year_sectors = profile_frame[profile_frame['Year']==0]['IPCC_2006_source_category'] # select sector with no year variability
    if sector in list(no_year_sectors): # set year == 0 if no yearly variability is present for the given sector
        year = 0
    if year > 2017: # use last available profile if year>2017
        year=2017
    
    temp_array = profile_frame[(profile_frame['IPCC_2006_source_category']==str(sector)) & 
                               (profile_frame['Year']==year)][['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']]
    temp_array=temp_array.to_numpy()
    out_array = np.zeros(12)
    for line in temp_array:
        out_array = out_array + line
    norm = np.sum(out_array)
    out_array = out_array/norm
    return out_array
    
def get_monthly_profile(yearly_emi, sectors, profile_region):
    profile_frame = pd.read_csv("./sector_monthly_gridmaps/EDGAR_temporal_profiles_r1.csv", sep = ' ')
    profile_frame = profile_frame[profile_frame['Region/country']==str(profile_region)]['IPCC_2006_source_category'] # select frame rows corresponding to the region of interest
    if sectors == 'ALL': # create a list that contains all sector names
        all_sectors = profile_frame
        all_sectors = list( dict.fromkeys(all_sectors) )
    else: 
        all_sectors = sectors
        
    region_frame = get_region_df(profile_region) # select dataframe for the given region
    all_time_profile = []
    for yr_emi in yearly_emi:
        yearly_profile = np.zeros(12)
        for sec in all_sectors:
            yearly_profile = yearly_profile + get_sector_profile(yr_emi[0], sec, region_frame)
        norm = np.sum(yearly_profile)
        yearly_profile = yearly_profile/norm
        # plt.scatter(np.arange(1,13,1), yearly_profile, s =2, label = str(yr_emi[0]))
        # plt.legend()
        yearly_profile = yearly_profile * yr_emi[1]
        all_time_profile.append(yearly_profile)
        
    return all_time_profile

def write_monthly_emissions(region, spec, monthly_profile):
  
    month_str=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August','September','October','November','December']
    out_file = open(par.predicted_res_folder+'predicted_'+region+'_'+spec+'_monthly_emi.txt','w')
    out_file.write('year month emi[t]\n')
    for year in range(21):
        for month in range(12):
            out_file.write(str(2000+year) +' '+month_str[month]+' '+str(monthly_profile[year][month])+'\n')
    out_file.close()
