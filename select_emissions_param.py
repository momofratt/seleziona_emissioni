import geopandas as gpd
import itertools
import shapely
import select_emissions_netcdf as net
import pandas as pd
#################################################
################ USER PARAMETERS ################
res_folder = './results/'
predicted_res_folder = './res_predicted_emi/'

lon_CMN, lat_CMN   = 10.70111, 44.19433 # CMN coordinates 
lon_bolo, lat_bolo = 11.3435 , 44.4947  # Bologna

# seleziona province. NB: le province selezionate DEVONO essere adiacenti. e.g. [Ferrara, Bologna, Modena] è una buona selezione, [Ferrara, Bologna, Parma] invece NO

# the variable description is reported in the Puy de Dome section

# ############## EMILIA ROMAGNA ###############
# region = 'ER' # suffix for ouput files
# region_full = 'Emilia Romagna'
# prov_numname = [(33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia'"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
#lat_stat, lon_stat = lat_CMN, lon_CMN
#country='Italy'
#stat_nm = 'Mt. Cimone'
# IPR=True
###############################################
# ############## TOSCANA ###############
# region = 'TOS' # suffix for ouput files
# region_full = 'Toscana'
# provinces = ['Grosseto', 'Livorno', 'Pisa', 'Massa-Carrara', 'Lucca', 'Pistoia', 'Prato', 'Firenze', 'Arezzo', 'Siena']
#lat_stat, lon_stat = lat_CMN, lon_CMN
#stat_nm = 'Mt. Cimone'
#country='Italy'
# IPR=True
###############################################
# ############## PO' VALLEY ###############
region = 'PO'
region_full = "Po' Valley"
prov_numname = [(3, 'Novara'),(12, 'Varese'),(13, 'Como'),(15, 'Milano'),(16, 'Bergamo'),(17, 'Brescia'),(18, 'Pavia'),(19, 'Cremona'),(20, 'Mantova'),(97, 'Lecco'),(98, 'Lodi'), (103, 'Monza e Brianza'),(23, 'Verona'),(24, 'Vicenza'),(26, 'Treviso'),(27, 'Venezia'),(28, 'Padova'),(29, 'Rovigo'),(33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
lat_stat, lon_stat = lat_CMN, lon_CMN
country='Italy'
stat_nm = 'Mt. Cimone'
IPR = True
###############################################
# ############## PO' VALLEY extended ###############
# region = 'PO-ext'
# region_full = "Po' Valley"
# prov_numname = [(97,'Lecco'),(14,'Sondrio'),(25, 'Belluno'),(22,'Trento'),
#                 (21,'Bozen'),(25,'Belluno'),(93,'Pordenone'), (30,'Udine'),
#                 (3, 'Novara'),(12, 'Varese'),(13, 'Como'),(15, 'Milano'),(16, 'Bergamo'),(17, 'Brescia'),
#                 (18, 'Pavia'),(19, 'Cremona'),(20, 'Mantova'),(97, 'Lecco'),(98, 'Lodi'), (103, 'Monza e Brianza'),
#                 (23, 'Verona'),(24, 'Vicenza'),(26, 'Treviso'),(27, 'Venezia'),(28, 'Padova'),(29, 'Rovigo'),
#                 (33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),
#                 (39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
# lat_stat, lon_stat = lat_CMN, lon_CMN
# country='Italy'
# stat_nm = 'Mt. Cimone'
# IPR = True
###############################################
###############################################
# ############## PO' VALLEY EST ###############
#altre province: (1, 'Turin'),(4, 'Cuneo')(2, 'Vercelli'),(96, 'Biella'),(5, 'Asti'),(25, 'Belluno'),(6, 'Alessandria'),
# region = 'PO-east'
# region_full = "Po' Valley EAST"
# prov_numname = [(17, 'Brescia'),(20, 'Mantova'),(23, 'Verona'),(24, 'Vicenza'),(26, 'Treviso'),(27, 'Venezia'),(28, 'Padova'),(29, 'Rovigo'),(35, "Reggio Emilia"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(25, 'Belluno'),(22,'Trento'),(21,'Bozen'),(25,'Belluno'),(93,'Pordenone'), (30,'Udine')]
# lat_stat, lon_stat = lat_CMN, lon_CMN
# stat_nm = 'Mt. Cimone'
# IPR = True
###############################################
###############################################
# ############## BAVARIA ###############
# altre province: (1, 'Turin'),(4, 'Cuneo')(2, 'Vercelli'),(96, 'Biella'),(5, 'Asti'),(25, 'Belluno'),(6, 'Alessandria'),
# region = 'GER'
# region_full = "Bayern-Baden-Württemberg"
# #prov_numname = [(3, 'Novara'),(12, 'Varese'),(13, 'Como'),(15, 'Milano'),(16, 'Bergamo'),(17, 'Brescia'),(18, 'Pavia'),(19, 'Cremona'),(20, 'Mantova'),(97, 'Lecco'),(98, 'Lodi'), (103, 'Monza e Brianza'),(23, 'Verona'),(24, 'Vicenza'),(26, 'Treviso'),(27, 'Venezia'),(28, 'Padova'),(29, 'Rovigo'),(33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
# prov_numname=[(0,'Bayern'),(0,'Baden-Württemberg')]
# lat_stat, lon_stat = 47.8011, 11.0246 
# stat_nm = 'Hohenpeissenberg'
# country = 'Germany' # country name in the admin1 geojson file so avoi selecting regions with same name from different countries.
# IPR = False
# ###############################################
# ###############################################
# # ############## Puy de dome ###############

#region = 'PUY' # name for the datafiles
#region_full = "Rhone-Auvergne" # name for the plots
#prov_numname=[(0,'Puy-de-Dôme'),(0,'Creuse'),(0,'Corrèze'),(0,'Cher'),(0,'Rhône'),(0,'Lot'),(0,'Lozère'),(0,'Aveyron'),(0,'Loire'),(0,'Haute-Loire'),(0,'Haute-Vienne'),(0,'Jura'),(0,'Saône-et-Loire'),(0,'Cantal'),(0,'Allier')]
#lat_stat, lon_stat = 45.7719, 2.9658 # latitude of the station. It is used to plot a triangle on that point
#stat_nm = 'Puy-de-Dôme' # name for the plots
#country = 'France' # country name in the admin1 geojson file so avoi selecting regions with same name from different countries.
#IPR = False

# ###############################################
# ###############################################
# ## Observatoire perenne de l'environnement ###

# region = 'OPE' # name for the datafiles
# region_full = "Grand-Est" # name for the plots
# prov_numname=[(0,'Meuse'),(0,'Bas-Rhin'),(0,'Vosges'),(0,'Marne'),(0,'Meurhe-et-Moselle'),(0,'Moselle'),(0,'Haute-Marne'),(0,'Ardennes'),(0,'Aube')]
# lat_stat, lon_stat = 45.7719, 2.9658 # latitude of the station. It is used to plot a triangle on that point
# stat_nm = 'Observatoire Perenne de l\'Environnement' # name for the plots
# country = 'France' # country name in the admin1 geojson file so avoi selecting regions with same name from different countries.
# IPR = False
###############################################

profile_region = 17 # region profile number according to EDGAR dataset (see Crippa et al. 2020)

start_year_ch4 = 1990 # first year of the CH4 EDGAR netcdf database
end_year_ch4   = 2018 # last year of the CH4 EDGAR netcdf database
start_year_co  = 1990 # first year of the CO EDGAR netcdf database
end_year_co    = 2018 # last year of the CO EDGAR netcdf database
end_year       = 2021 # last year for emission prediction
end_year_IPR   = 2019 # last year of ISPRA data
#################################################
#################################################

#################################################
################ FIXED PARAMETERS ###############
d_lon = 0.1 # separation of grid points in deg
d_lat = 0.1
earth_rad = 6371 # [km]
j_CMN, i_CMN = net.lat_to_index(lat_CMN, lon_CMN)
lat_CMN, lon_CMN = net.index_to_lat(j_CMN, i_CMN) # obtain approximated lat_CMN, lon_CMN cendered on grid points 

prov_num  = [prov[0] for prov in prov_numname] # province number according to ISPRA dataset
provinces = [prov[1] for prov in prov_numname] # province name according to admin1 geoframe
########################################

###############################################################################
###                         Select emission region                          ###
###############################################################################
# defines other parameters and variables according to the provided user parameters

states_geoframe = gpd.GeoDataFrame.from_file('admin1.geojson')        # read geojson file

prov_geoframe = gpd.GeoDataFrame()  # select provinces

for prov in provinces:
    new_prov = states_geoframe[(states_geoframe['name']==prov) & (states_geoframe['country']==country)] # read line from the GeoDataFrame that correspond to the selected province
    if type(new_prov['geometry'].to_list()[0]) == shapely.geometry.multipolygon.MultiPolygon: # if prov is a multipolygon (e.g. islands) uses only the largest polygon 
        prov_poly = new_prov['geometry'].to_list()    
        max_prov_poly = max(prov_poly[0], key=lambda a: a.area)
        new_prov['geometry'] = max_prov_poly
    
    # prov_geoframe = prov_geoframe.append(new_prov, ignore_index=True)
    prov_geoframe = pd.concat([prov_geoframe, new_prov],ignore_index=True)


# union of provinces and evaluation of external boundaries
geoms = prov_geoframe['geometry'].tolist()
intersection_iter = gpd.GeoDataFrame(gpd.GeoSeries([poly[0].union(poly[1]) for poly in  itertools.combinations(geoms, 2)]), columns=['geometry'])
regione = intersection_iter.unary_union
bound_coords = list(regione.exterior.coords)

bound_coords_index = [net.lat_to_index(coord[1], coord[0]) for coord in bound_coords] # convert boundary coordinates to indexes
bound_coords_index = list(dict.fromkeys(bound_coords_index))                      # remove duplicate elements in the list

max_lat = max(b[1] for b in bound_coords)
min_lat = min(b[1] for b in bound_coords)
max_lon = max(b[0] for b in bound_coords)
min_lon = min(b[0] for b in bound_coords)

max_coords_ind = net.lat_to_index(max_lat, max_lon)
min_coords_ind = net.lat_to_index(min_lat, min_lon)

nbin = max(abs(max_coords_ind[0] - min_coords_ind[0]), abs(max_coords_ind[1] - min_coords_ind[1])) #maximum extension in number of bins between N-S and E-W directions
nbin = nbin+10 # add bins to avoid problems at the borders

if nbin % 2 != 0:
    nbin=nbin+1    
nbin_2=int(0.5*nbin) # half width of nbin

j_ref, i_ref = net.lat_to_index(0.5*(max_lat -min_lat)+min_lat, 0.5*(max_lon-min_lon)+min_lon)  # reference coordinates to center the sel_ data matrix
lat_ref, lon_ref = net.index_to_lat(j_ref, i_ref)
