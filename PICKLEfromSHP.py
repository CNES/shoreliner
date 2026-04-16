import geopandas as gpd
import pickle
import numpy as np
import sys
sys.path.insert(0, r"/work/scratch/data/graffim/shorelines/Codes/shoreliner_main/Codes")
from functions import Tools

def create_transect_dict(shp_file, output_pickle_file, prefix="PF",proj=False,projepsg=4326):
    """
    Convert a shapefile with transects to a dictionary and save it as a pickle file.

    Parameters:
    shp_file (str): Path to the input .shp file containing transects.
    output_pickle_file (str): Path to the output .p file for the dictionary.
    prefix (str): Prefix for the dictionary keys (default is 'PF').
    """
    # Load the shapefile using GeoPandas
    gdf = gpd.read_file(shp_file)
    
    # Initialize an empty dictionary to store transects
    transect_dict = {}
    
    # Iterate through each row in the GeoDataFrame
    for idx, row in gdf.iterrows():
        # Extract the geometry (assumes LineString)
        geometry = row.geometry
        if geometry.geom_type == "LineString" and len(geometry.coords) >= 2:
            # Get the first and last coordinates of the LineString as (lon, lat)
            point1 = geometry.coords[0]
            point2 = geometry.coords[-1]
            
            # Create a key for the transect (e.g., PF1, PF2, ...)
            transect_key = f"{prefix}{idx + 1}"
            
            # Store the two points in the dictionary under the 'transect' key
            transect_dict[transect_key] = {
                "transect": [(point1[0], point1[1]), (point2[0], point2[1])]
            }
    if proj==True:
        for i in transect_dict:
            transect_dict[i]['transect_proj'] = transect_dict[i]['transect']
        for i in transect_dict:
            transect_dict[i]['transect'] = Tools.convert_epsg(np.array(transect_dict[i]['transect_proj']),projepsg,4326)[:,:2][:,::-1]   
        
    # Save the dictionary to a pickle file
    print()
    with open(output_pickle_file, 'wb') as p_file:
        print(p_file)
        pickle.dump(transect_dict, p_file)
    
    print(f"Transect dictionary saved to {output_pickle_file}")

# Example usage
shp_file = r"C:\Users\pc\Downloads\transects.shp"  # Replace with the actual path to your .shp file
output_pickle_file = r"C:\Users\pc\Downloads\AugustusTransects.p"       # Replace with your desired output .p file name
create_transect_dict(shp_file, output_pickle_file,proj=True,projepsg = 32613)