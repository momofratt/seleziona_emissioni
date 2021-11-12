#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:41:43 2021

@author: Cosimo Fratticioli 
@contact: c.fratticioli@isac.cnr.it
"""
# la selezione può essere effettuata anche utilizzando le province ed i rispettivi codici:
# province = [('Piacenza',33), ('Parma',34),('Reggio Emilia',35),('Modena',36),('Bologna',37),('Ferrara',38),('Ravenna',39),('Forli',40),('Rimini',99)]
import pandas as pd    

def select_ISPRA(region):
    # returns ISPRA total emissions for CO and CH4 in Emilia Romagna 
    # i dati devono essere suddivisi in files diferenti. Ogni contiene i dati relativi ad un macrosettore
    # I vari macrosettori sono riportati nei dati disaggregati ISPRA ( http://www.sinanet.isprambiente.it/it/sia-ispra/inventaria/disaggregazione-dellinventario-nazionale-2015/view )
    if region == 'ER':
        region_nm = 'Emilia Romagna'
        region_id = 8
    elif region == 'TOS':
        region_nm = 'Toscana'
        region_id = 9
    data_path = './data_ISPRA/'
    macro_settore = [1,2,3,4,5,6,7,8,9,10,11]
    regioni = [(region_nm, region_id)]  # ad ogni regione è associato un numero identificativo- Per selezionare altre regioni bisogna indicare il rispettivo numero
    anni = [1990, 1995, 2000, 2005, 2010, 2015]
    
    def select_regional_data():
        # seleziona dati riferiti alla regione di interesse
        sel_data_frame = pd.DataFrame()
        for reg in regioni:
            for sett in macro_settore: 
                if sett<10: # definisce il prefisso del file realativo al macrosettore
                    prefix = 'sect_0'
                else:
                    prefix = 'sect_'
                frame = pd.read_csv(data_path+prefix+str(sett)+'.csv', sep = ' ')    # legge frame del macrosettore    
                frame = frame[frame['COD_REGI'] == reg[1]]  # seleziona dati regionali per il macrosettore
                sel_data_frame = sel_data_frame.append(frame) # appende i dati del macrosettore al frame totale
        return sel_data_frame
    
    def select_specie(in_frame, specie):
        # seleziona dati per la specie selezionata
        out_frame = pd.DataFrame()
        out_frame = out_frame.append(in_frame[in_frame['COD_POL']==specie])
        return out_frame
        
    regional_frame = select_regional_data()
    
    regional_ch4_frame = select_specie(regional_frame, '004')
    regional_co_frame = select_specie(regional_frame, '005')
    
    co_annuale  = []
    ch4_annuale = []
    for yr in anni:
        co_annuale.append( regional_co_frame[str(yr)].sum())
        ch4_annuale.append(regional_ch4_frame[str(yr)].sum())

    return co_annuale, ch4_annuale, anni