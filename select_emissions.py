#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 11:21:26 2021

@author: Cosimo Fratticioli
@contact: c.fratticioli@isac.cnr.it
"""

import read_ISPRA
import select_emissions_fit_eval as fitev
import select_emissions_param as par
import select_emissions_plot as plot
import select_emissions_EDGAR_profiles_sectors as eps
import pandas as pd
###############################################################################
###   Select emissions inside emission region and write output on files     ###
###############################################################################
# define arrays to store selected data with respective latitudes and longitudes
# evaluate arrays for emission and relative error 
emi_CH4, emi_err_CH4 = fitev.eval_emission('CH4')
emi_CO, emi_err_CO = fitev.eval_emission('CO')

emi_err_CH4 = plot.eval_max_emi_error(emi_CH4, emi_err_CH4) # set the error equal to the max between the 10% of the emission (EDGAR default????) and the error estimated from boundary selection
emi_err_CO  = plot.eval_max_emi_error(emi_CO, emi_err_CO)

if par.IPR:
    print('read ISPRA')
    # get ISPRA emissions and write on file
    emi_CO_ISPRA, emi_CH4_ISPRA, anni = read_ISPRA.select_ISPRA(par.prov_num)
    
    file_emi_ch4_ISPRA = open('./res_files/total_emissions_'+par.region+'_CH4_ISPRA.txt', 'w')
    file_emi_ch4_ISPRA.write('year tot_emi[t]\n')
    file_emi_co_ISPRA = open('./res_files/total_emissions_'+par.region+'_CO_ISPRA.txt', 'w')
    file_emi_co_ISPRA.write('Year emi_CO[t]\n')
    for i in range(len(emi_CH4_ISPRA)):
        file_emi_co_ISPRA.write(str(anni[i]) +' '+ str(round(emi_CO_ISPRA[i])) +'\n')
        file_emi_ch4_ISPRA.write(str(anni[i]) +' '+ str(round(emi_CH4_ISPRA[i])) +'\n')
    
    file_emi_co_ISPRA.close()
    file_emi_ch4_ISPRA.close()

###############################################################################
###                               FIT                                       ###
###############################################################################
print('fit emission')
fitev.fit_emis('CO', emi_CO, emi_err_CO, par.start_year_co, par.end_year_co, 2, True)
fitev.fit_emis('CH4', emi_CH4, emi_err_CH4, par.start_year_ch4, par.end_year_ch4, 4, True)
fitev.fit_emis('CO', emi_CO, emi_err_CO, 2000, 2015, 3, True)

print('eval predicted')
fitev.eval_predicted_values(emi_CO, emi_err_CO, 'CO')
fitev.eval_predicted_values(emi_CH4,emi_err_CH4, 'CH4')

print('eval predicted ISPRA')
if par.IPR:
    fitev.fit_emis_ISPRA('CO', emi_CO_ISPRA,  1990, 2015, 2, True)
    fitev.eval_predicted_values_ISPRA(emi_CO_ISPRA, par.region, 'CO')
    fitev.eval_predicted_values_ISPRA(emi_CH4_ISPRA, par.region, 'CH4')


###############################################################################
###                               PLOT                                      ###
###############################################################################

if par.IPR:
    print('plot ISPRA EDGAR')
    plot.plot_ISPRA_EDGAR(emi_CO, emi_err_CO, emi_CO_ISPRA, 'CO', anni)
    plot.plot_ISPRA_EDGAR(emi_CH4, emi_err_CH4, emi_CH4_ISPRA, 'CH4', anni)    

print('plot predicted')
plot.plot_predicted(emi_CO, emi_err_CO, emi_CO_ISPRA,'CO')
plot.plot_predicted(emi_CH4, emi_err_CH4, emi_CH4_ISPRA,'CH4')

print('Plot 2D')
plot.plot_2D_emis('CH4', 2015, False)


###############################################################################
###                           MONTHLY PROFILES                              ###
###############################################################################
print('plot monthly')
yrs_profile = [x for x in range(2000,2018)]

yearly_emi_CO_frame = pd.read_csv(par.predicted_res_folder+'predicted_'+par.region+'_CO_yearly_emi.txt', sep=' ',usecols=['year','emi[t]'])
yearly_emi_CO_frame = yearly_emi_CO_frame[yearly_emi_CO_frame['year']>1999]
yearly_emi_CO = list(zip(yearly_emi_CO_frame.year, yearly_emi_CO_frame['emi[t]'] ))

monthly_profile_CO = eps.get_monthly_profile( yearly_emi_CO ,'ALL', profile_region=par.profile_region) 
eps.write_monthly_emissions(par.region, 'CO', monthly_profile_CO)
plot.plot_EDGAR_monthly(emi_CO, emi_err_CO, monthly_profile_CO, 'CO')


yearly_emi_CH4_frame = pd.read_csv(par.predicted_res_folder+'predicted_'+par.region+'_CH4_yearly_emi.txt', sep=' ',usecols=['year','emi[t]'])
yearly_emi_CH4_frame = yearly_emi_CH4_frame[yearly_emi_CH4_frame['year']>1999]
yearly_emi_CH4 = list(zip(yearly_emi_CH4_frame.year, yearly_emi_CH4_frame['emi[t]'] ))

monthly_profile_CH4 = eps.get_monthly_profile( yearly_emi_CH4 ,'ALL', profile_region=par.profile_region) 
eps.write_monthly_emissions(par.region, 'CH4', monthly_profile_CH4)
plot.plot_EDGAR_monthly(emi_CH4, emi_err_CH4, monthly_profile_CH4, 'CH4')


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