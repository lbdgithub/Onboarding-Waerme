import os
import geopandas as gpd
import pandas as pd
from pathlib import Path
from shapely import wkt, Point, LineString, MultiPolygon, MultiLineString
from shapely.geometry import Polygon
import re
import xml.etree.ElementTree as ET
import math
import numpy as np





def remove_z_coordinate(wkt_string):
    # Use regular expression to find and remove the Z coordinate and the "Z" character
    pattern = re.compile(r'Z|\b(\d+\.\d+)\s(\d+\.\d+)\s(\d+\.\d+)\b')
    result = re.sub(pattern, r'\1 \2', wkt_string)
    return result





def extract_ground_surfaces(citygml_version,xml_data, path_results, epsg, layer):
    # Parse the XML data
    tree = ET.parse(xml_data)
    root = tree.getroot()

    # Define a GeoDataFrame    
    gdf_data = pd.DataFrame(columns=('alkis_id','nr','function', 'height_m', 'roof_type', 'surface', 'surface_id','community_key','community_name', 'postcode', 'street', 'house_number', 'geometry'))
    
    nr = 1
    
    # Iterate through Building elements
    for building_elem in root.findall('.//{{http://www.opengis.net/citygml/building/{version}}}Building'.format(version=citygml_version)):
    
        # Extract attributes on building level
        gml_id = building_elem.get('{http://www.opengis.net/gml}id')

		#extract building function for nested attributes
        function = None
        for function_elem in building_elem.findall('./{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="function")):
            function = function_elem.text

		#extract building height for nested attributes
        height = None
        for measured_height_elem in building_elem.findall('./{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="measuredHeight")):
            height = float(measured_height_elem.text)

		#extract roof type for nested attributes
        roofType = None
        for roofType_elem in building_elem.findall('./{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="roofType")):

            roofType = roofType_elem.text

        #extract gemeindeschluessel in stringAttribute element
        #define the namespaces
        namespaces = {
    		'gen': 'http://www.opengis.net/citygml/generics/{version}'.format(version=citygml_version),
    		'core': 'http://www.opengis.net/citygml/{version}'.format(version=citygml_version),
			}
        
        gemeindeschluessel = None
        for gemeindeschluessel_elem in building_elem.findall('.//gen:stringAttribute[@name="Gemeindeschluessel"]/gen:value', namespaces=namespaces):
            gemeindeschluessel = gemeindeschluessel_elem.text

        #extract city, street, housenumber and postal code
        #define the namespaces
        namespaces = {
    		'xAL': 'urn:oasis:names:tc:ciq:xsdschema:xAL:2.0',
    		'core': 'http://www.opengis.net/citygml/{version}'.format(version=citygml_version),
    		'bldg': 'http://www.opengis.net/citygml/building/{version}'.format(version=citygml_version)
			}

        #extract city
        city = None
        for city_elem in building_elem.findall('.//bldg:address/core:Address/core:xalAddress/xAL:AddressDetails/xAL:Country/xAL:Locality/xAL:LocalityName', namespaces=namespaces):
            city = city_elem.text

        #extract street
        street_name = None
        for street_name_elem in building_elem.findall('.//bldg:address/core:Address/core:xalAddress/xAL:AddressDetails/xAL:Country/xAL:Locality/xAL:Thoroughfare/xAL:ThoroughfareName', namespaces=namespaces):
            street_name = street_name_elem.text

        #extract housenumber
       	housenumber = None
        for housenumber_elem in building_elem.findall('.//bldg:address/core:Address/core:xalAddress/xAL:AddressDetails/xAL:Country/xAL:Locality/xAL:Thoroughfare/xAL:ThoroughfareNumber', namespaces=namespaces):
            housenumber = housenumber_elem.text

        #extract postal code
        postcode = None
        for postcode_elem in building_elem.findall('.//bldg:address/core:Address/core:xalAddress/xAL:AddressDetails/xAL:Country/xAL:Locality/xAL:PostalCode/xAL:PostalCodeNumber', namespaces=namespaces):
            postcode = postcode_elem.text

		#if gml consists of building parts rooftype will be none; in this case go to else and iterate the building parts
        if roofType != None:

        	#iterate through surfaces within building: ground vs. roof (layer)
        	for surface_elem in building_elem.findall('.//{{http://www.opengis.net/citygml/building/{version}}}{layer}'.format(version=citygml_version, layer=layer)):

        		surface_geometry = []
        
        		#iterate through surface parts
        		for pos_list_elem in surface_elem.findall('.//{http://www.opengis.net/gml}posList'):

        			#extract coordinates from posList                
        			coordinates = [float(coord) for coord in pos_list_elem.text.split()]
        			coordinates_str = str(pos_list_elem.text)
                
        			#extract the surface id
        			surface_id = surface_elem.get('{http://www.opengis.net/gml}id')

        			#write extracted data to dataframe row
        			gdf_data.loc[nr] = [str(gml_id),nr,str(function),str(height), str(roofType), str(layer), str(surface_id), str(gemeindeschluessel), str(city), str(postcode), str(street_name), str(housenumber) ,str(coordinates_str)]
        			nr = nr + 1
        			#print(str(layer) + "__" + str(nr))
                		
        else:
        
        	for building_part_elem in building_elem.findall('.//{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="BuildingPart")):
        		
        		#extract measured height from building parts
        		for measured_height_elem in building_part_elem.findall('.//{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="measuredHeight")):
        			height = float(measured_height_elem.text)

        		#extract roof type from building parts
        		for roofType_elem in building_part_elem.findall('.//{{http://www.opengis.net/citygml/building/{version}}}{element}'.format(version=citygml_version, element="roofType")):

        			roofType = roofType_elem.text

				#iterate the surfaces within building parts
        		for surface_elem in building_part_elem.findall('.//{{http://www.opengis.net/citygml/building/{version}}}{layer}'.format(version=citygml_version, layer=layer)):

        			surface_geometry = []
        
        			#iterate through surface parts
        			for pos_list_elem in surface_elem.findall('.//{http://www.opengis.net/gml}posList'):

        				#extract coordinates from posList                
        				coordinates = [float(coord) for coord in pos_list_elem.text.split()]
        				coordinates_str = str(pos_list_elem.text)
                
        				#extract the surface id
        				surface_id = surface_elem.get('{http://www.opengis.net/gml}id')

        				#write extracted data to dataframe row
        				gdf_data.loc[nr] = [str(gml_id),nr,str(function),str(height), str(roofType), str(layer), str(surface_id), str(gemeindeschluessel), str(city), str(postcode), str(street_name), str(housenumber) ,str(coordinates_str)]
        				nr = nr + 1
        				#print(str(layer) + "__" + str(nr))
        
    #prepare geometry string for wkt-format
    gdf_data['geometry'] = gdf_data['geometry'].apply(lambda x: ' '.join([val + ',' if i % 3 == 2 else val for i, val in enumerate(x.split())])).str[:-1]

    """
	#manipulate first and last coordinate set (equalize them) to close the linear ring for polygon transformation
    gdf_data['value_split_first'] = gdf_data['geometry'].str.rsplit(', ').str[0]
    #print(gdf_data['value_split_first'])
    gdf_data['all_without_last_coordinates'] = gdf_data['geometry'].str.rsplit(', ', n=1).str[0]
    gdf_data['geometry'] = gdf_data['all_without_last_coordinates'].astype(str)+ ", " + gdf_data['value_split_first'].astype(str) 
    #print(gdf_data.geometry)
    
    #delete temporary columns
    gdf_data = gdf_data.drop(['value_split_first'], axis=1)
    gdf_data = gdf_data.drop(['all_without_last_coordinates'], axis=1)
    """
	
	#transform coordinates to polygon from linear ring
    gdf_data.geometry = "POLYGON Z (( "+gdf_data.geometry+" ))"
    
    #remove z-coordinate from geometry-string
    gdf_data.geometry =  gdf_data.geometry.apply(lambda x: remove_z_coordinate(x))
    gdf_data['geometry'] = gdf_data['geometry'].apply(wkt.loads)
    gdf = gpd.GeoDataFrame(gdf_data, geometry='geometry').reset_index(drop=True)

    gdf.crs = epsg
    
    #calculate floor area
    if layer == "GroundSurface":
    	gdf["floor_area_m2"] =  gdf["geometry"].area
    else:
    	pass

    #add community name and number if data is None
    gdf = add_community_name_and_number_if_none(gdf,gdf_radiation_ger,epsg)
    
    #export results to geojson
    #gdf['geometry'] = gdf['geometry'].buffer(0)
    if not gdf.empty:
    	gdf.to_file(path_results, driver='GeoJSON')
    return gdf





def change_file_extension_from_xml_to_gml(path_to_directory):
	for filename in os.listdir(path_to_directory):
		infilename = os.path.join(path_to_directory,filename)
		if not os.path.isfile(infilename): continue
		oldbase = os.path.splitext(filename)
		newname = infilename.replace('.xml', '.gml')
		output = os.rename(infilename, newname)





def create_list(path,extension):
	lst = os.listdir(path)
	lst = list(filter(lambda k: extension in k, lst))
	lst = [ f for f in os.listdir(path) if f[(len(f) - len(extension)):len(f)].find(extension)>=0 ]
	#print(lst)
	return lst





def concat_geojson(path,lst):
	nr = 1
	for i in lst:
		if nr == 1:
			gdf = gpd.read_file(path+i)
			#print(gdf)
		else:
			gdf_tmp = gpd.read_file(path+i)
			gdf = pd.concat([gdf, gdf_tmp], axis=0, ignore_index=True)
		#print(gdf)
		nr = nr + 1
	return gdf





def dropping_geom_duplicates(gdf):
	print(str(len(gdf.index)) + "(length BEFORE drop duplicates)")
	gdf = gdf.drop_duplicates(subset="geometry", keep="first")
	print(str(len(gdf.index)) + "(length AFTER drop duplicates)")
	return gdf





def add_community_name_and_number_if_none(gdf,gdf_radiation_ger,epsg):
	#change projection to gml data protection
	right_df = gdf_radiation_ger.to_crs(epsg)
	#perform spatial join between gml data and radiation table ger with community names and numbers
	gdf_sjoin = gpd.sjoin(gdf, right_df[['geometry','community_key_ger']], how='left', predicate='intersects')
	#keep only the first matching element for each row in left_df
	gdf_sjoin = gdf_sjoin.reset_index(drop=True)
	gdf_sjoin = gdf_sjoin.drop_duplicates(subset="geometry", keep='first')
	#replace if key is None
	gdf_sjoin.community_key = np.where(gdf_sjoin.community_key == 'None', gdf_sjoin.community_key_ger, gdf_sjoin.community_key )
	gdf_sjoin = gdf_sjoin.drop(['community_key_ger'], axis=1)
	gdf_sjoin = gdf_sjoin.drop(['index_right'], axis=1)
	#add name if name is none or isdigit
	gdf_sjoin = pd.merge(gdf_sjoin,right_df[['community_key_ger','community_name_ger']],left_on="community_key", right_on="community_key_ger", how="left")
	gdf_sjoin.community_name = np.where(gdf_sjoin.community_name == 'None', gdf_sjoin.community_name_ger, gdf_sjoin.community_name )
	gdf_sjoin.community_name = np.where(gdf_sjoin.community_name.str.isdigit(), gdf_sjoin.community_name_ger, gdf_sjoin.community_name )
	gdf_sjoin = gdf_sjoin.drop(['community_name_ger'], axis=1)
	gdf_sjoin = gdf_sjoin.drop(['community_key_ger'], axis=1)
	return gdf_sjoin
	




def calculate_roof_orientation(roof_gdf, ground_gdf,epsg,gdf_radiation_ger):
    #calculate  centroid of roof geometry
    #roof_gdf.crs = epsg
    #ground_gdf.crs = epsg
    #roof_gdf = roof_gdf.to_crs('epsg:4326')
    gdf_centroid = roof_gdf.copy()
    
    gdf_centroid.geometry = gdf_centroid.geometry.centroid
    gdf_centroid = gdf_centroid.to_crs('epsg:4326')
    gdf_centroid['x'] = gdf_centroid.geometry.x
    gdf_centroid['y'] = gdf_centroid.geometry.y

    #calculate centroid of ground surface
    
    ground_gdf.geometry = ground_gdf.geometry.centroid
    ground_gdf = ground_gdf.to_crs('epsg:4326')
    ground_gdf['x_reference'] = ground_gdf.geometry.x
    ground_gdf['y_reference'] = ground_gdf.geometry.y
    
    #merge roof and ground
    gdf_centroid = pd.merge(gdf_centroid, ground_gdf[['alkis_id','x_reference','y_reference']], on="alkis_id", how="left")

    #calculate angular orientation between centroid of roof to centroid of ground surface
    gdf_centroid['angle'] = np.arctan2(np.radians(gdf_centroid["y_reference"]) - np.radians(gdf_centroid["y"]),
                                        np.radians(gdf_centroid["x_reference"]) - np.radians(gdf_centroid["x"]))
	#calcualte bearing from angle
    gdf_centroid.loc[gdf_centroid.angle >= 0, 'bearing'] = np.degrees(gdf_centroid['angle'] + ((0 * math.pi)))
    gdf_centroid.loc[gdf_centroid.angle < 0, 'bearing'] = np.degrees(gdf_centroid['angle'] + ((2 * math.pi)))
    gdf_centroid['bearing'] = -(gdf_centroid['bearing'].astype(float) + 90 -360) #+ (math.pi/2)
    gdf_centroid['bearing'] = np.where( gdf_centroid['bearing']<0, gdf_centroid['bearing'] + 360,gdf_centroid['bearing'])

    #add cardinal direction
    #https://dev.qweather.com/en/docs/resource/wind-info/
    conditions = [
        (gdf_centroid['bearing'] >= 0) & (gdf_centroid['bearing'] < 22.5),
        (gdf_centroid['bearing'] >= 22.5) & (gdf_centroid['bearing'] < 45),
        (gdf_centroid['bearing'] >= 45) & (gdf_centroid['bearing'] < 67.5),
        (gdf_centroid['bearing'] >= 67.5) & (gdf_centroid['bearing'] < 90),
        (gdf_centroid['bearing'] >= 90) & (gdf_centroid['bearing'] < 112.5),
        (gdf_centroid['bearing'] >= 112.5) & (gdf_centroid['bearing'] < 135),
        (gdf_centroid['bearing'] >= 135) & (gdf_centroid['bearing'] < 157.5),
        (gdf_centroid['bearing'] >= 157.5) & (gdf_centroid['bearing'] < 180),
        (gdf_centroid['bearing'] >= 180) & (gdf_centroid['bearing'] < 202.5),
        (gdf_centroid['bearing'] >= 202.5) & (gdf_centroid['bearing'] < 225),
        (gdf_centroid['bearing'] >= 225) & (gdf_centroid['bearing'] < 247.5),
        (gdf_centroid['bearing'] >= 247.5) & (gdf_centroid['bearing'] < 270),
        (gdf_centroid['bearing'] >= 270) & (gdf_centroid['bearing'] < 292.5),
        (gdf_centroid['bearing'] >= 292.5) & (gdf_centroid['bearing'] < 315),
        (gdf_centroid['bearing'] >= 315) & (gdf_centroid['bearing'] < 337.5),      
        (gdf_centroid['bearing'] >= 337.5) & (gdf_centroid['bearing'] <= 360)
    ]

    directions = ["N", "NO", "NO", "O", "O", "SO", "SO", "S", "S", "SW", "SW", "W", "W", "NW", "NW", "N"]

    gdf_centroid['direction'] = np.select(conditions, directions)
    gdf_centroid['direction'] = np.where(gdf_centroid['roof_type'] == "1000", "flat", gdf_centroid['direction'])
    
    #add global yearly radiation to dataset
    gdf_centroid = gpd.sjoin(gdf_centroid,gdf_radiation_ger[['geometry','yearly_radiation_mean']], how='left', predicate='intersects')

    #merge the results back with the original roof_gdf    
    result_gdf = pd.merge(roof_gdf, gdf_centroid[['nr', 'angle', 'bearing', 'direction', 'yearly_radiation_mean']], on="nr", how="left")
    
    #finalize gdf
    result_gdf = result_gdf.drop(['angle'], axis=1)
    result_gdf = dropping_geom_duplicates(result_gdf)
    result_gdf = result_gdf.to_crs(epsg)
    return result_gdf





def calculate_roof_area(surface_area_2d, slope_angle_degrees):
    #convert slope angle from degrees to radians
    slope_angle_radians = math.radians(slope_angle_degrees)
    #calculate roof area 
    roof_area = surface_area_2d / math.cos(slope_angle_radians)
    return roof_area





def calc_solar_factor(slope, azimuth, df_solar_factor):
	#add solar-factor based on direction/azimuth and slop
    df = df_solar_factor[df_solar_factor.azimuth == azimuth].reset_index(drop=True)
    slope_round = slope#.astype(float).round(0)
    val = df[f'slope_{slope_round}'][0]
    return val





def pv_operating_hours(radiation,solar_factor,df_solar_technical_input):
	#get technical input data
	pv_efficiency = df_solar_technical_input[df_solar_technical_input['variable'] == "pv_efficiency"].reset_index(drop=True).value_input.iloc[0]
	inverter_efficiency = df_solar_technical_input[df_solar_technical_input['variable'] == "inverter_efficiency"].reset_index(drop=True).value_input.iloc[0]
	performance_ratio = df_solar_technical_input[df_solar_technical_input['variable'] == "performance_ratio"].reset_index(drop=True).value_input.iloc[0]
	deduction_factor = df_solar_technical_input[df_solar_technical_input['variable'] == "deduction_factor"].reset_index(drop=True).value_input.iloc[0]
	spec_pv_size = 1 / pv_efficiency
	#calculate op hours
	pv_operating_hours = radiation * spec_pv_size * solar_factor * pv_efficiency * inverter_efficiency * performance_ratio	
	return pv_operating_hours





def estimate_pv_power(gdf,df_solar_technical_input):
	#calculate pv power
	deduction_factor = df_solar_technical_input[df_solar_technical_input['variable'] == "deduction_factor"].reset_index(drop=True).value_input.iloc[0]
	pv_efficiency = df_solar_technical_input[df_solar_technical_input['variable'] == "pv_efficiency"].reset_index(drop=True).value_input.iloc[0]
	spec_pv_size = 1 / pv_efficiency
	#estimate pv_power
	gdf['roof_area_qualified_m2'] = gdf['roof_area_m2'] * gdf['roof_qualification']
	gdf['pv_power'] = gdf['roof_area_qualified_m2'] * (1-deduction_factor) / spec_pv_size
	gdf['pv_power'] = gdf['pv_power'].round(0)
	gdf['pv_production'] = gdf['pv_operating_hours'] * gdf['pv_power']
	return gdf





def estimate_solar_thermal(gdf,df_solar_technical_input):
	#share of area... is a indicative assumption (technical correct would be to identify the household size...)
	#see also p. 44: Leitfaden KWP BW: https://um.baden-wuerttemberg.de/fileadmin/redaktion/m-um/intern/Dateien/Dokumente/2_Presse_und_Service/Publikationen/Energie/Leitfaden-Kommunale-Waermeplanung-barrierefrei.pdf
	share_of_area_solar_heat_collectorsize = df_solar_technical_input[df_solar_technical_input['variable'] == "share_of_area_solar_heat_collectorsize"].reset_index(drop=True).value_input.iloc[0]
	spec_solar_heat_power = df_solar_technical_input[df_solar_technical_input['variable'] == "spec_solar_heat_power"].reset_index(drop=True).value_input.iloc[0]
	gdf['solar_heat_area_m2'] = gdf['roof_area_qualified_m2'] * share_of_area_solar_heat_collectorsize
	gdf['solar_heat_production'] = gdf['solar_heat_area_m2'] * spec_solar_heat_power * gdf['pv_operating_hours'] / 1000 #in kWh/a
	gdf['solar_heat_power'] = gdf['solar_heat_production'] / gdf['pv_operating_hours']
	return gdf





def roof_qualification(gdf, df_roof_qualification):
	#identify roof directions from results by direction for exclusion in aggregations on buildings (in future: shall single roofs be excluded simply by direction (e.g. N) or by threshold of operating hours)
	gdf = pd.merge(gdf,df_roof_qualification, on="direction", how="left" )
	return gdf





def aggregate_pv_solar_roof_results_on_alkisid(gdf):
	#perform aggregation
	gdf_agg = gdf.copy()
	gdf_agg = gdf.groupby(['alkis_id']).agg({'roof_area_ground_m2':'sum', \
										 'roof_area_m2':'sum', \
										 'roof_area_qualified_m2':'sum', \
										 'pv_power':'sum', \
										 'pv_production':'sum', \
										 'solar_heat_area_m2':'sum', \
										 'solar_heat_production':'sum', \
										 'solar_heat_power':'sum'}).reset_index()	
	gdf_agg['pv_operating_hours'] = gdf_agg['pv_production'] / gdf_agg['pv_power']
	
	#rename column and add "total"
	gdf_agg.rename({'roof_area_ground_m2': 'total_roof_area_ground_m2', \
					'roof_area_m2': 'total_roof_area_m2',
					'roof_area_qualified_m2': 'total_roof_area_qualified_m2',
					'pv_power': 'total_pv_power', \
					'pv_production': 'total_pv_production', \
					'pv_operating_hours': 'total_pv_operating_hours', \
					'solar_heat_area_m2': 'total_solar_heat_area_m2', \
					'solar_heat_production': 'total_solar_heat_production', \
					'solar_heat_power': 'total_solar_heat_power'}, axis=1, inplace=True)

	#merge to original dataframe
	cols = ['alkis_id','total_roof_area_ground_m2','total_roof_area_m2','total_roof_area_qualified_m2','total_pv_power','total_pv_production','total_pv_operating_hours','total_solar_heat_area_m2','total_solar_heat_production','total_solar_heat_power']
	gdf_agg = gdf_agg[cols]
	gdf = pd.merge(gdf, gdf_agg, on="alkis_id", how="left")
	gdf.total_pv_operating_hours = np.where(gdf.total_pv_production == 0, gdf.pv_operating_hours, gdf.total_pv_operating_hours)
	return gdf





def calculate_capex(gdf):
	#calcualte capex for pv and solar
	spec_pv_invest = df_solar_technical_input[df_solar_technical_input['variable'] == "spec_pv_invest"].reset_index(drop=True).value_input.iloc[0]
	spec_solar_invest = df_solar_technical_input[df_solar_technical_input['variable'] == "spec_solar_invest"].reset_index(drop=True).value_input.iloc[0]
	gdf['total_pv_invest_euro'] = spec_pv_invest * gdf['total_pv_power']
	gdf['total_solar_invest_euro'] =spec_solar_invest * gdf['total_solar_heat_area_m2']
	return gdf





def roof_solar_calculation(gdf_roof, df_rooftype_slope, df_solar_factor):
	#merge dataframes by roof_type
	df_rooftype_slope.roof_type = df_rooftype_slope.roof_type .astype(str)
	gdf_roof = pd.merge(gdf_roof,df_rooftype_slope,on="roof_type", how="left")

	#set exceptions
	gdf_roof.slope = np.where(gdf_roof.slope.isna(),30,gdf_roof.slope)
	gdf_roof.roof_name = np.where(gdf_roof.roof_name.isna(),"Sonstiges",gdf_roof.roof_name)
	
	#calculate solar_factor
	gdf_roof['solar_factor'] = gdf_roof.apply(lambda x: calc_solar_factor(x.slope, x.direction, df_solar_factor), axis=1)
	
	#calculate roof area
	gdf_roof['roof_area_ground_m2'] = gdf_roof.geometry.area
	gdf_roof['roof_area_m2'] = gdf_roof.apply(lambda x: calculate_roof_area(x.roof_area_ground_m2, x.slope), axis=1)
	return gdf_roof





def technical_dimensioning_pv_solar(gdf_roof,df_solar_technical_input, df_roof_qualification):
	#ACHTUNG: ab hier müsste in Zukunft in INFRA gerechnet werden, weil hier editierbare Inputs mit betrachtet werden -> welche Anlage mit welchen technischen Prämissen etc... --> Prämissen/Szenariorahmen
	#add roof qualification by direction (further improvement: by threshold of operating hours)
	gdf_roof = roof_qualification(gdf_roof, df_roof_qualification)
	#calculate pv operating hours
	gdf_roof['pv_operating_hours'] = gdf_roof.apply(lambda x: pv_operating_hours(x.yearly_radiation_mean,x.solar_factor,df_solar_technical_input), axis=1)
	#estimate the possible pv power
	gdf_roof = estimate_pv_power(gdf_roof,df_solar_technical_input)
	#estimate solar thermal power and production
	gdf_roof = estimate_solar_thermal(gdf_roof,df_solar_technical_input)
	#aggregate data on buildings/groundsurface by alkis_id
	gdf_roof = aggregate_pv_solar_roof_results_on_alkisid(gdf_roof)
	#calculate capex
	gdf_roof = calculate_capex(gdf_roof)
	return gdf_roof
	
	
    


def change_datatypes_for_dataframe(gdf):	
	#datatypes according to infra onboarding documentation
	gdf.height_m = gdf.height_m.astype(float)
	gdf.function = gdf.function.astype(str)
	gdf.roof_type = gdf.roof_type.astype(int)
	gdf.community_key = gdf.community_key.astype(str)
	gdf.postcode = gdf.postcode.astype(str)
	return gdf

	


	
def finalize_roof_solar_data(gdf):
	cols = ['bearing','yearly_radiation_mean','roof_area_ground_m2','roof_area_m2','pv_operating_hours',
			'roof_area_qualified_m2','pv_power', 'pv_production',
       		'solar_heat_area_m2', 'solar_heat_production', 'solar_heat_power',
       		'total_roof_area_ground_m2', 'total_roof_area_m2',
       		'total_roof_area_qualified_m2', 'total_pv_power', 'total_pv_production',
       		'total_pv_operating_hours', 'total_solar_heat_area_m2',
       		'total_solar_heat_production', 'total_solar_heat_power',
       		'total_pv_invest_euro', 'total_solar_invest_euro']
	gdf[cols] = gdf[cols].round(4)
	return gdf




		
def create_infra_onboarding_file_lod1(gdf, avg_floorheight_m, construction_year):
	#set the standard for the infra onboarding file for lod1 data
	#calculate floor count
	gdf['floor_count'] = gdf['height_m'].astype(float) / avg_floorheight_m
	gdf['floor_count'] = np.floor(gdf['floor_count'])
	gdf['floor_count'] = np.maximum(1,gdf['floor_count'])
	#create additional mandatory field
	gdf['district_name'] = None
	gdf['construction_year'] = construction_year
	
	#filter only buildings
	gdf['function_prefix'] = gdf['function'].str.split('_').str[0]
	gdf['function_prefix'] = gdf['function_prefix'].astype(str)
	gdf['function'] = gdf['function'].str.split('_').str[1]
	gdf = gdf[gdf.function_prefix == "31001"]
	
	#select columns
	cols = ['alkis_id',
			'geometry',
			'height_m',
			'floor_count',
			'floor_area_m2',
			'function',
			'roof_type',
			'construction_year',
			'community_key',
			'community_name',
			'postcode',
			'street',
			'house_number',
			'district_name']
	gdf = gdf[cols]

	#error handlung
	gdf['community_key'] = gdf['community_key'].replace('None', 0) 

	#set datatypes
	gdf['geometry'] = gdf['geometry'].apply(lambda x: MultiPolygon([x]) if x.geom_type == 'Polygon' else x)
	gdf = gdf.astype({'alkis_id': 'str', \
					'geometry': 'str', \
					'height_m': 'float', \
					'floor_count': 'int', \
					'floor_area_m2': 'float', \
					'function': 'int', \
					'roof_type': 'int', \
					'construction_year': 'int', \
					'community_key': 'int', \
					'community_name': 'str', \
					'postcode': 'str', \
					'street': 'str', \
					'house_number': 'str', \
					'district_name': 'str'})

	# ids have to be unique for processing in INFRA. Therefore, alkis_ids occurring multiple times are modified
	multiple_alkis_ids_df = gdf.groupby('alkis_id').filter(lambda x: len(x) > 1).drop_duplicates(subset='alkis_id')
	print(multiple_alkis_ids_df.info())
	for i, row_i in multiple_alkis_ids_df.iterrows():
		current_alkis_id = row_i['alkis_id']
		count = 0
		for j, row_j in gdf.loc[(gdf['alkis_id'] == current_alkis_id)].iterrows():
			count += 1
			gdf.loc[(gdf['geometry'] == row_j['geometry']), 'alkis_id'] = current_alkis_id + '_' + str(count).zfill(2)

	#export as csv
	gdf.to_csv(path_results + "__INFRA_Onboarding_BuildingsLOD1.csv", index=False)





def round_coordinates(geom, precision=6):
    if geom.geom_type == 'Point':
        return Point(round(geom.x, precision), round(geom.y, precision))
    elif geom.geom_type == 'LineString':
        return LineString([(round(x, precision), round(y, precision)) for x, y in geom.coords])
    elif geom.geom_type == 'Polygon':
        exterior = [(round(x, precision), round(y, precision)) for x, y in geom.exterior.coords]
        interiors = [[(round(x, precision), round(y, precision)) for x, y in interior.coords] for interior in geom.interiors]
        return Polygon(exterior, interiors)
    elif geom.geom_type == 'MultiPolygon':
        polygons = [round_coordinates(p, precision) for p in geom.geoms]
        return MultiPolygon(polygons)
    elif geom.geom_type == 'MultiLineString':
        return MultiLineString([round_coordinates(line, precision) for line in geom])
    else:
        # Handle other geometry types if needed
        pass
        #return geom





def create_infra_onboarding_file_lod2(gdf_roof):
	#set the standard for the infra onboarding file for lod2 data - Potenziale
	#filter relevant attributes
	gdf = gdf_roof.copy()

	#filter only buildings
	gdf['function_prefix'] = gdf['function'].str.split('_').str[0]
	gdf['function_prefix'] = gdf['function_prefix'].astype(str)
	gdf['function'] = gdf['function'].str.split('_').str[1]
	gdf = gdf[gdf.function_prefix == "31001"]

	cols = ['geometry','alkis_id','roof_type','direction','yearly_radiation_mean','total_roof_area_m2','total_roof_area_qualified_m2', \
			'total_pv_power', 'total_pv_production', \
			'total_pv_operating_hours', \
			'total_solar_heat_area_m2','total_solar_heat_production', 'total_solar_heat_power', \
			'total_pv_invest_euro','total_solar_invest_euro']
	gdf = gdf[cols]
	
	#rename accordingly to standard of upload file
	gdf.rename({'alkis_id': 'gml_id',
					'yearly_radiation_mean': 'globalradiation_kwh_per_m2',  \
					'total_roof_area_m2': 'area_roof_m2',
					'total_roof_area_qualified_m2': 'area_roof_available_m2',
					'total_pv_operating_hours': 'operating_hours_h_per_year', \
					'total_pv_power': 'pv_power_kw', \
					'total_pv_production': 'pv_production_kwh', \
					'total_solar_heat_area_m2': 'solar_heat_area_m2', \
					'total_solar_heat_power': 'solar_heat_power_kw', \
					'total_solar_heat_production': 'solar_heat_production_kwh', \
					'total_pv_invest_euro': 'pv_invest_euro', \
					'total_solar_invest_euro': 'solar_invest_euro'}, axis=1, inplace=True)
	
	#make valid
	gdf['geometry'] = gdf['geometry'].buffer(0)
	
	#change geom precision to 6 digits
	gdf = gdf.dropna(subset=['geometry'])
	gdf = gdf[~gdf['geometry'].is_empty & gdf['geometry'].notna()]
	gdf['geometry'] = gdf['geometry'].apply(round_coordinates)

	#export
	cols_pv = ['geometry','gml_id','direction','globalradiation_kwh_per_m2','roof_type','area_roof_m2','area_roof_available_m2','operating_hours_h_per_year','pv_power_kw','pv_production_kwh','pv_invest_euro']
	cols_solar = ['geometry','gml_id','direction','globalradiation_kwh_per_m2','roof_type','area_roof_m2','area_roof_available_m2','solar_heat_area_m2','operating_hours_h_per_year','solar_heat_power_kw','solar_heat_production_kwh','solar_invest_euro']
	gdf[cols_pv].to_file(path_results + '__INFRA_Onboarding_PV_Roof.geojson', driver='GeoJSON')
	gdf[cols_solar].to_file(path_results + '__INFRA_Onboarding_Solar_Roof.geojson', driver='GeoJSON')





if __name__ == "__main__":
	#ToDo: ggf. nochmal das neueste offizielle Shape zu Gemeinden runterladen und die Datengrundlage des radiation-files überarbeiten anpassen (/input/radiation_postcode_community.geojson) --> Stimmen die Namen noch, gab es Gebietsreformen etc...
	#alkis_id muss unique sein... polygone müssen zu multipolygonen werden... Problem: z.B. 2 Gebäudegrundflächen mit derselben id, haben unetrschiedliche dachtypen und gebäudehöhen...




	#TODO: toal op hours passen tlw. nicht zu der stromerzeugung: problem in zeile 438-437????


	#----------------------------------------
	#----------------------------------------
	#NOTES:
	#CityGML version and EPSG must be set according to the input GML dataset! 
	#Example: 
	#Baden-Württemberg EPSG:25832, Version 1.0
	#Berlin EPSG:25833, Version 2.0
	#Brandenburg EPSG:25833, Version 2.0
	#Hamburg EPSG:25832, Version 1.0
	#Niedersachsen EPSG:25832, Version 1.0
	#Sachsen EPSG:25833, Version 1.0
	#Thüringen EPSG:25832, Version 1.0
	#
	#OPTIONAL:
	#-filter_community_names_bool: if buildings should be filtered by community name this variable must be set to true
	# and the relevant communitys must be listed in community_names_list
	#-replace_unknown_building_function_bool: if too many buildings have an unknown function this value can be set to
	# true, and they will be replaced with target_building_function
	#----------------------------------------
	#----------------------------------------




	#set the gml version
	citygml_version = '1.0'


	
	#set the epsg of the gml data source
	epsg = 'epsg:25833'



	#set the community names to filter by
	filter_community_names_bool = True
	community_names_list = ['Neuensalz', 'Treuen']



	#set parameters for replacing unknown building functions
	replace_unknown_building_function_bool = True
	target_building_function = '31001_1000'
	min_floor_area_m2 = 40



	#define input variables
	path_base = os.getcwd() + '/'
	path_input = path_base + 'lod2/'
	path_results = path_input + '_processed_data/'
	path_input_solar_calculation = path_base + 'input/'
	target_epsg_infra_onboarding = 'epsg:4326'



	#set variables for infra onboarding files
	avg_floorheight_m = 3.8
	construction_year = 2000 #only if data not given. additional function needs to be integrated if external data is available and can be joined to gml lod data



	#load solar calculation input data
	gdf_radiation_ger = gpd.read_file(path_input_solar_calculation + "radiation_postcode_community.geojson")
	df_rooftype_slope = pd.read_csv(path_input_solar_calculation + "roof_type_slope.csv")
	df_solar_factor = pd.read_csv(path_input_solar_calculation + "solar_factor.csv")
	df_solar_technical_input = pd.read_csv(path_input_solar_calculation + "input_variables.csv")
	df_roof_qualification = pd.read_csv(path_input_solar_calculation + "roof_qualification.csv")
	
	
	
	#change file extension from xml to gml
	change_file_extension_from_xml_to_gml(path_input)
	
	
	
	#start parser and convert gml-files to geojson
	print("")
	print("----- START PARSING GML -----")
	print("")
	surface_layer_lst = ['GroundSurface','RoofSurface']
	gml_files = create_list(path_input,'.gml')
	
	

	#iterate gml-data in input folder
	file_nr = 1
	for file in gml_files[::1]:
		print("")
		print("********  " + str(file_nr) + "  ********")
		print("********  " + str(file) + "  ********")
		print("")
		file_nr = file_nr + 1
		path_in = path_input + file
		
		#GroundSurface
		surface_layer = "GroundSurface"
		path_out = path_results + file
		path_out = path_out.split('.')[0] + '__'+str(surface_layer)+'.geojson'
		gdf_ground = extract_ground_surfaces(citygml_version,path_in,path_out,epsg,surface_layer)
		if gdf_ground.empty:
			print("GML data is empty!")
			continue

		#RoofSurfaces
		surface_layer = "RoofSurface"
		path_out = path_results + file
		path_out = path_out.split('.')[0] + '__'+str(surface_layer)+'.geojson'
		gdf_roof = extract_ground_surfaces(citygml_version,path_in,path_out,epsg,surface_layer)

		#RoofSurfaces: add roof_directions to roof gdf
		gdf_roof = calculate_roof_orientation(gdf_roof, gdf_ground,epsg,gdf_radiation_ger)
		
		#RoofSurfaces: solar calculation
		gdf_roof = roof_solar_calculation(gdf_roof, df_rooftype_slope, df_solar_factor)
		
		#RoofSurfaces: technical dimensioning
		gdf_roof = technical_dimensioning_pv_solar(gdf_roof, df_solar_technical_input, df_roof_qualification)

		#change datatypes according to infra onboarding documentation
		gdf_ground = change_datatypes_for_dataframe(gdf_ground)
		gdf_roof = change_datatypes_for_dataframe(gdf_roof)
		gdf_roof = finalize_roof_solar_data(gdf_roof)

		#export results as geojson
		path_out = path_results + file
		path_out = path_out.split('.')[0] + '__'+str(surface_layer)+'.geojson'
		gdf_roof.to_file(path_out, driver='GeoJSON')



	#merge parsed data to one geojson-file for each surface
	#iterate
	print("")
	print("----- START MERGING GEOJSON -----")
	print("")
	for surface_layer in surface_layer_lst:

		#merge as geojson processed gml files
		parser_result_files = create_list(path_results,'__'+str(surface_layer)+'.geojson')
		gdf = concat_geojson(path_results,parser_result_files)
		
		#fix geometries
		gdf.geometry = gdf.geometry.buffer(0)
		
		#set epsg for infra onboardings files
		gdf = gdf.to_crs(target_epsg_infra_onboarding)

		#filter buildings that aren't in relevant communitys
		if filter_community_names_bool:
			print(str(len(gdf.index)) + "(length BEFORE filter community names)")
			gdf = gdf.loc[gdf['community_name'].isin(community_names_list)]
			print(str(len(gdf.index)) + "(length AFTER filter community names)")

		#check percentage of buildings with unknown function
		percentage_unknown_function = gdf.loc[(gdf['function'] == '31001_9998'), 'function'].shape[0]/gdf.shape[0] * 100
		if (surface_layer == "GroundSurface") and (percentage_unknown_function >= 30):
			if replace_unknown_building_function_bool:
				gdf.loc[(gdf['function'] == '31001_9998') & (gdf['floor_area_m2'] < min_floor_area_m2), 'function'] = target_building_function
				print(percentage_unknown_function, '% of buildings had an unknown function, which was replaced with the function', target_building_function)
			else:
				print('WARNING:', percentage_unknown_function, '% of buildings have an unknown function.')
		
		#export merged data from gml
		gdf.to_file(path_results + '__' + str(surface_layer) + '__.geojson', driver='GeoJSON')

		#create infra onboarding file
		if surface_layer == "GroundSurface":
			create_infra_onboarding_file_lod1(gdf, avg_floorheight_m, construction_year)
		else:
			create_infra_onboarding_file_lod2(gdf)

	
		
		
	
	
	
