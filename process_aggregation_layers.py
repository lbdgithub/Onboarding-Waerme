import geopandas as gpd
import os
import fiona
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon


def process_districts(path_in, path_res, target_epsg):
    # read input data
    districts_gdf = gpd.read_file(path_in + 'ortsteile.shp')
    districts_gdf = districts_gdf.to_crs(target_epsg)
    # print(districts_gdf.info())

    # filter aggregation_layers
    districts_filtered_gdf = districts_gdf.loc[
        (districts_gdf['GEMEINDE'] == 'Treuen, Stadt') | (districts_gdf['GEMEINDE'] == 'Neuensalz')]
    districts_reduced_gdf = districts_filtered_gdf.loc[:, ['SCHLüSSEL', 'ORTSTEIL', 'geometry']]
    districts_reduced_gdf = districts_reduced_gdf.rename(columns={'SCHLüSSEL': 'district_key', 'ORTSTEIL': 'name'})

    # create municipalitys
    municipality_gdf = districts_filtered_gdf[['GEMEINDE', 'geometry']]
    municipality_grouped_gdf = municipality_gdf.dissolve(by='GEMEINDE')
    municipality_grouped_gdf.geometry = municipality_grouped_gdf.geometry.explode(index_parts=True)[0:].values
    municipality_grouped_gdf['name'] = municipality_grouped_gdf.index
    municipality_grouped_gdf.insert(0, 'municipality_key', [14523270, 14523430], True)
    municipality_grouped_gdf.insert(0, 'area_ha', [3349.0, 4374.0], True)
    municipality_grouped_gdf.insert(0, 'population', [2045, 7704], True)

    # save results
    districts_reduced_gdf.to_file(path_res + '__INFRA_Onboarding_districts.geojson', driver='GeoJSON')
    municipality_grouped_gdf.to_file(path_res + '__INFRA_Onboarding_municipalitys.geojson', driver='GeoJSON')


def process_building_blocks(path_in, path_res, target_epsg):
    # read input data
    building_blocks_neuensalz_gdf = gpd.read_file(path_in + 'Neuensalz_UEbergabe.gpkg', layer='1a_Neuensalz_mittlGasV')
    building_blocks_neuensalz_reduced_gdf = building_blocks_neuensalz_gdf.loc[:, ['geometry']]
    building_blocks_neuensalz_reduced_gdf = building_blocks_neuensalz_reduced_gdf.to_crs(target_epsg)
    building_blocks_treuen_gdf = gpd.read_file(path_in + 'Treuen_UEbergabe.gpkg', layer='1a_Treuen_mittlGasV')
    building_blocks_treuen_reduced_gdf = building_blocks_treuen_gdf.loc[:, ['geometry']]
    building_blocks_treuen_reduced_gdf = building_blocks_treuen_reduced_gdf.to_crs(target_epsg)

    dataframe_list = [building_blocks_neuensalz_reduced_gdf, building_blocks_treuen_reduced_gdf]
    building_blocks_gdf = gpd.GeoDataFrame(pd.concat(dataframe_list, ignore_index=True))
    building_blocks_gdf.to_file(path_res + '__INFRA_Onboarding_building_blocks.geojson', driver='GeoJSON')


def create_list(path, extension):
    lst = os.listdir(path)
    lst = list(filter(lambda k: extension in k, lst))
    lst = [ f for f in os.listdir(path) if f[(len(f) - len(extension)):len(f)].find(extension)>=0 ]
    #print(lst)
    return lst


def get_layernames(path):
    layernames = []
    for layername in fiona.listlayers(path):
        with fiona.open(path, layer=layername) as src:
            layernames.append([layername, len(src)])


def process_alkis_data(path_in, path_res, target_epsg):
    xml_files = create_list(path_in, '.xml')
    dataframe_list = []

    for file in xml_files[::1]:
        # open file
        path_input_file = path_in + file
        current_file_gdf = gpd.read_file(path_input_file, layer='AX_Flurstueck')
        # extract relevant information
        current_file_reduced_gdf = current_file_gdf.loc[:, ['flurstueckskennzeichen', 'geometry']]
        current_file_reduced_gdf = current_file_reduced_gdf.rename(columns={'flurstueckskennzeichen': 'estate_key'})
        # append
        dataframe_list.append(current_file_reduced_gdf)

    # merge and save
    fluerstuecke_gdf = gpd.GeoDataFrame( pd.concat( dataframe_list, ignore_index=True) ).set_crs('epsg:25833')
    fluerstuecke_gdf = fluerstuecke_gdf.to_crs(target_epsg)

    fluerstuecke_gdf.to_file(path_res + '__INFRA_Onboarding_flurstuecke.geojson', driver='GeoJSON')





if __name__ == "__main__":
    # define input variables
    path_input = 'aggregation_layers/'
    path_results = path_input + '_processed_data/'
    target_epsg_infra_onboarding = 'epsg:4326'

    process_districts_bool = False
    process_building_blocks_bool = False
    process_alkis_bool = True

    if process_districts_bool:
        print("----- START PROCESSING DISTRICTS -----")
        process_districts(path_input, path_results, target_epsg_infra_onboarding)
        print('')

    if process_building_blocks_bool:
        print("----- START PROCESSING BUILDING BLOCKS -----")
        process_building_blocks(path_input, path_results, target_epsg_infra_onboarding)
        print('')

    if process_alkis_bool:
        print("----- START PROCESSING ALKIS DATA -----")
        process_alkis_data('ALKIS/', 'ALKIS/_processed_data/', target_epsg_infra_onboarding)
        print('')



