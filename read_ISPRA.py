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
import matplotlib.pyplot as plt

def select_specie(in_frame, specie):
    # seleziona dati per la specie selezionata
    out_frame = pd.DataFrame()
    out_frame = out_frame.append(in_frame[in_frame['COD_POL']==specie])
    return out_frame

def select_regional_data(regioni, macro_settore):
    # seleziona dati riferiti alla regione di interesse
    sel_data_frame = pd.DataFrame()
    for reg in regioni:
        if type(macro_settore)==int:
            macro_settore = [macro_settore]
        for sett in macro_settore: 
            if sett<10: # definisce il prefisso del file realativo al macrosettore
                prefix = 'sect_0'
            else:
                prefix = 'sect_'
            frame = pd.read_csv( './data_ISPRA/'+prefix+str(sett)+'.csv', sep = ' ')    # legge frame del macrosettore    
            frame = frame[frame['COD_REGI'] == reg[1]]  # seleziona dati regionali per il macrosettore
            sel_data_frame = sel_data_frame.append(frame) # appende i dati del macrosettore al frame totale
    return sel_data_frame
    

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
    
    macro_settore = [1,2,3,4,5,6,7,8,9,10,11]
    regioni = [(region_nm, region_id)]  # ad ogni regione è associato un numero identificativo- Per selezionare altre regioni bisogna indicare il rispettivo numero
    anni = [1990, 1995, 2000, 2005, 2010, 2015]
        
    regional_frame = select_regional_data(regioni, macro_settore)
    
    regional_ch4_frame = select_specie(regional_frame, '004')
    regional_co_frame = select_specie(regional_frame, '005')
    
    co_annuale  = []
    ch4_annuale = []
    for yr in anni:
        co_annuale.append( regional_co_frame[str(yr)].sum())
        ch4_annuale.append(regional_ch4_frame[str(yr)].sum())

    return co_annuale, ch4_annuale, anni

def plot_regional_sources():
    
    macro_settore = [1,2,3,4,5,6,7,8,9,10,11]
    ER  = [('Emilia Romagna', 8)]
    TOS = [('Toscana',9)]
    anno = '2015'
    
    settori_TOS_CO = []
    settori_ER_CO = []
    settori_TOS_CH4 = []
    settori_ER_CH4 = []
    
    for sett in macro_settore:  
        frame_ER = select_regional_data(ER, sett)
        frame_TOS = select_regional_data(TOS, sett)
    
        ch4_frame_ER = select_specie(frame_ER, '004')
        co_frame_ER = select_specie(frame_ER, '005')
        ch4_frame_TOS = select_specie(frame_TOS, '004')
        co_frame_TOS = select_specie(frame_TOS, '005') 
        
        settori_ER_CH4.append(  ch4_frame_ER[anno].sum() )
        settori_ER_CO.append(    co_frame_ER[anno].sum())
        settori_TOS_CH4.append(ch4_frame_TOS[anno].sum())
        settori_TOS_CO.append(  co_frame_TOS[anno].sum())
    
    fig, ax = plt.subplots(2,2, figsize = (8,8))
    fig.suptitle('Sector contribution to CO and CH$_4$ emissions for TOS and ER during '+anno+'\nfrom ISPRA database', x=0.63)
    
    ax[0][0].pie(settori_ER_CH4)
    ax[0][0].set_title('ER CH$_4$')
    
    ax[0][1].pie(settori_ER_CO)
    ax[0][1].set_title('ER CO')
    
    ax[1][0].pie(settori_TOS_CH4)
    ax[1][0].set_title('TOS CH$_4$')

    labels =  ['Combustione nell industria e impianti energetici',
               'Impianti di combustione non industriale',
               'Processi produttivi (combustione nell industria manufatturiera)',
               'Processi produttivi (combustione senza contatto)',
               'Estrazione e distribuzione di combustibili fossili ed energia geotermica',
               'Uso di solventi ed altri prodotti',
               'Trasporti stradali',
               'Altre sorgenti mobili e macchinari mobili (trasporti fuori strada)',
               'Trattamento dei rifiuti e discariche',
               'Agricoltura',
               'Altre emissioni ed assorbimenti']
    
    ax[1][1].pie(settori_TOS_CO, labels=labels, labeldistance=None)
    ax[1][1].set_title('TOS CO')
    plt.subplots_adjust(bottom=0.3, right=0.94, left= 0.3, top=0.9)
    fig.legend(loc='lower right')

    plt.savefig('./results/pie_chart_sectors_ER_TOS.pdf', format ='pdf', bbox_inches="tight")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    