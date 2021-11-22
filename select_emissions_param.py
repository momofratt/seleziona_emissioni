import geopandas as gpd
import itertools
import shapely
import select_emissions_netcdf as net
#################################################
################ USER PARAMETERS ################
res_folder = './results/'
predicted_res_folder = './res_predicted_emi/'
# seleziona province. NB: le province selezionate DEVONO essere adiacenti. e.g. [Ferrara, Bologna, Modena] è una buona selezione, [Ferrara, Bologna, Parma] invece NO

# ############## EMILIA ROMAGNA ###############
# region = 'ER' # suffix for ouput files
# region_full = 'Emilia Romagna'
# prov_numname = [(33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia'"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
# IPR=True
###############################################
# ############## TOSCANA ###############
# region = 'TOS' # suffix for ouput files
# region_full = 'Toscana'
# provinces = ['Grosseto', 'Livorno', 'Pisa', 'Massa-Carrara', 'Lucca', 'Pistoia', 'Prato', 'Firenze', 'Arezzo', 'Siena']
# IPR=True
###############################################
# ############## PO' VALLEY ###############
region = 'PO'
region_full = "Po' Valley"
prov_numname = [(1, 'Turin'),(2, 'Vercelli'),(3, 'Novara'),(4, 'Cuneo'),(5, 'Asti'),(6, 'Alessandria'),(96, 'Biella'),(103, 'Verbano-Cusio-Ossola'),(12, 'Varese'),(13, 'Como'),(14, 'Sondrio'),(15, 'Milano'),(16, 'Bergamo'),(17, 'Brescia'),(18, 'Pavia'),(19, 'Cremona'),(20, 'Mantova'),(97, 'Lecco'),(98, 'Lodi'), (103, 'Monza e Brianza'),(23, 'Verona'),(24, 'Vicenza'),(25, 'Belluno'),(26, 'Treviso'),(27, 'Venezia'),(28, 'Padova'),(29, 'Rovigo'),(33, 'Piacenza'),(34, 'Parma'),(35, "Reggio Emilia"),(36, 'Modena'),(37, 'Bologna'),(38, 'Ferrara'),(39, 'Ravenna'),(40, "Forlì-Cesena"),(99, 'Rimini')]
IPR = True

prov_num  = [prov[0] for prov in prov_numname] # province number according to ISPRA dataset
provinces = [prov[1] for prov in prov_numname] # province name according to admin1 geoframe
########################################
profile_region = 17 # region profile number according to EDGAR dataset (see Crippa et al. 2020)

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
j_CMN, i_CMN = net.lat_to_index(lat_CMN, lon_CMN)
lat_CMN, lon_CMN = net.index_to_lat(j_CMN, i_CMN) # obtain approximated lat_CMN, lon_CMN cendered on grid points 


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
