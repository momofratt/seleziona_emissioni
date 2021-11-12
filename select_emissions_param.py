import geopandas as gpd
import itertools
import shapely
import select_emissions_param as par
import select_emissions_netcdf as net
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

for prov in par.provinces:
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
