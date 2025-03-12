from shapely.geometry import LineString, Point
import geopandas as gpd
import os 

def transectToSHPfile(transects, current_path, crs=3857):
    """
    Create a .shp file of the transects to visualize it on QGIS
    """
    trsc_lines = [LineString(transects[tr]['transect_proj']) for tr in transects]
    gdftr = gpd.GeoDataFrame(geometry=trsc_lines,crs=crs)
    gdftr.to_file(os.path.join(current_path,'TRSC.shp'))

def waterlineToSHPfile(wl_tmp, tag_idx, seg_method, image_name, current_path, crs=3857):
    """
    Create a .shp file of the waterline to visualize it on QGIS
    """
    if not('WL_'+tag_idx in os.listdir(current_path)):
        os.mkdir('WL_'+tag_idx)
    if not(seg_method in os.listdir(os.path.join(current_path, 'WL_'+tag_idx))):
        os.mkdir(os.path.join(current_path,'WL_'+tag_idx,seg_method))
    save_path = os.path.join(current_path,'WL_'+tag_idx,seg_method)
    points = [Point(p) for p in wl_tmp]
    gdf = gpd.GeoDataFrame(geometry=points,crs=crs)
    gdf.to_file(os.path.join(save_path,image_name[:-4]+'_'+seg_method+'_WL.shp'))

