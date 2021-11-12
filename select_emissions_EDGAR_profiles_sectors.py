#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 11:22:41 2021

@author: cosimo
"""
import pandas as pd
import numpy as np

def get_sector_profile(year, sector, profile_region):
    profile_frame = pd.read_csv("./sector_monthly_gridmaps/EDGAR_temporal_profiles_r1.csv", sep = ' ')
    profile_frame = profile_frame[profile_frame['Region/country']==str(profile_region)]
    
    no_year_sectors = profile_frame[profile_frame['Year']==0]['IPCC_2006_source_category'] # select sector with no year variability
    if sector in list(no_year_sectors): # set year == 0 if no yearly variability is present for the given sector
        year = 0
    
    temp_array = profile_frame[(profile_frame['IPCC_2006_source_category']==str(sector)) & 
                               (profile_frame['Year']==year)][['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']]
    temp_array=temp_array.to_numpy()
    out_array = np.zeros(12)
    for line in temp_array:
        out_array = line
    yearly_sum = np.sum(out_array)
    out_array = out_array/yearly_sum

    return out_array
    
def get_monthly_profile(yearly_emi, sectors, profile_region):
    profile_frame = pd.read_csv("./sector_monthly_gridmaps/EDGAR_temporal_profiles_r1.csv", sep = ' ')
    profile_frame = profile_frame[profile_frame['Region/country']==str(profile_region)]['IPCC_2006_source_category'] # select frame rows corresponding to the region of interest
    if sectors == 'ALL': # create a list that contains all sector names
        all_sectors = profile_frame
        all_sectors = list( dict.fromkeys(all_sectors) )
    else: 
        all_sectors = sectors
    all_time_profile = []
    for yr_emi in yearly_emi:
        yearly_profile = np.zeros(12)
        for sec in all_sectors:
            yearly_profile = yearly_profile + yr_emi[1] * get_sector_profile(yr_emi[0], sec, profile_region)
        all_time_profile.append(yearly_profile)
        print(yr_emi[0])
    return all_time_profile
    