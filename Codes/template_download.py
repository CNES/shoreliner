#---------------------------------------------------------------------------------------------
#--------------- Index maps generation and extraction from GEE -------------------------------
#--------------- Marcan Graffin (largely inspired by Kilian Vos' CoastSat tool) -------------
#---------------------------------------------------------------------------------------------
import ee
import os
import pickle
import sys
sys.path.insert(0, r"../../../Codes/functions")
import yaml
from yaml.loader import SafeLoader
import GEE, Tools


global epsg_target,polygon_geom,S2,L5,L7,L8,L9
file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

if not(inputs['Waterline'] or inputs['SandBar'] or inputs['Vegetataion']):
    sys.exit('No Feature to extract, please modify the config.yaml file')



useS2,useL9,useL8,useL7,useL5=GEE.satMissions(inputs['Missions'])

ee.Initialize(project=inputs['EE_ID'])

start,end = GEE.dates(inputs['Dates'])



filepath = './'
polygon=pickle.load(open('poly.json','rb'))
transects=pickle.load(open('transects.p','rb'))
try:
    lonEPSG = transects[list(transects.keys())[len(transects)//2]]['transect'][0][0]
    latEPSG = transects[list(transects.keys())[len(transects)//2]]['transect'][0][1]
except:
    lonEPSG=(inputs['ROI_array'][0][0]+inputs['ROI_array'][2][0])//2
    latEPSG=(inputs['ROI_array'][0][1]+inputs['ROI_array'][2][1])//2
epsg_target = Tools.convert_wgs_to_utm(lonEPSG, latEPSG)


GEE.folderCreation(inputs)

if inputs['Waterline']:
    wl_f = GEE.index_f(inputs['WaterlineIndex'])

if inputs['SandBar']:
    sb_f = GEE.index_f(inputs['SandBarIndex'])

if inputs['Vegetation']:
    vg_f = GEE.index_f(inputs['VegetationIndex'])

if not('RGB' in os.listdir()) and inputs['RGB']:
    os.mkdir('RGB')

polygon_geom=ee.Geometry.Polygon(polygon)
polygon_coll = ee.FeatureCollection(ee.Feature(polygon_geom))

if inputs['Waterline']:
    wl_data=ee.ImageCollection([])
if inputs['SandBar']:
    sb_data=ee.ImageCollection([])
if inputs['Vegetation']:
    vg_data=ee.ImageCollection([])
if inputs['RGB']:
    rgb_data=ee.ImageCollection([])
    

def add_quality_score(image):
    valid_pixels = image.select('red').reduceRegion(
       reducer=ee.Reducer.count(),
       geometry=polygon_geom,
       scale=image.get('res'),maxPixels=1e13
       ).get('red')
    return image.set('valid_pixels', valid_pixels)

def best_image_per_dateS2(img):
    date = img.get('simple_date')
    same_day = S2.filter(ee.Filter.eq('simple_date', date)).sort('valid_pixels', False)
    return same_day.first()

def best_image_per_dateL9(img):
    date = img.get('simple_date')
    same_day = L9.filter(ee.Filter.eq('simple_date', date)).sort('valid_pixels', False)
    return same_day.first()

def best_image_per_dateL8(img):
    date = img.get('simple_date')
    same_day = L8.filter(ee.Filter.eq('simple_date', date)).sort('valid_pixels', False)
    return same_day.first()

def best_image_per_dateL7(img):
    date = img.get('simple_date')
    same_day = L7.filter(ee.Filter.eq('simple_date', date)).sort('valid_pixels', False)
    return same_day.first()

def best_image_per_dateL5(img):
    date = img.get('simple_date')
    same_day = L5.filter(ee.Filter.eq('simple_date', date)).sort('valid_pixels', False)
    return same_day.first()

def ensureProjS2(image):
    return image.reproject(
        **{
      'crs':'EPSG:'+epsg_target,
      'scale': 10}).copyProperties(image,['res','date','satellite'])

def ensureProjLandsat(image):
    return image.reproject(
        **{
      'crs':'EPSG:'+epsg_target,
      'scale': 15}).copyProperties(image,['res','date','satellite'])
#Image Collection generation and SCoWI computation 
#%%S2 ini

if useS2:
    S2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
    S2 = S2.filterDate(start,end).filter(ee.Filter.bounds(polygon_coll)).filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",inputs['MaxCloudCover']))
    if inputs['SpecSeasonalPeriod']:
        S2 = S2.filter(ee.Filter.calendarRange(inputs['SpecDates'][0],inputs['SpecDates'][1],'month'))
    S2 = S2.map(GEE.normalizeBandNamesS2).map(add_quality_score).map(GEE.add_simple_date)
    if inputs['FilterNonFullImages']:
        S2 = GEE.filterNonFullImages(S2,inputs)
    S2_distinct=S2.distinct(['simple_date'])
    S2 = ee.ImageCollection(S2_distinct.map(best_image_per_dateS2))
    S2 = GEE.resampleL5S2(S2,inputs)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(S2.map(wl_f).map(ensureProjS2))
    if inputs['SandBar']:
        sb_data=sb_data.merge(S2.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(S2.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(S2.map(GEE.RGB).map(ensureProjS2))
    
    
    
    
#%%L8 ini

if useL8:
    L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
    L8 = L8.filterDate(start,end).filter(ee.Filter.bounds(polygon_coll)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    if inputs['SpecSeasonalPeriod']:
        L8 = L8.filter(ee.Filter.calendarRange(inputs['SpecDates'][0],inputs['SpecDates'][1],'month'))
    L8 = L8.map(GEE.normalizeBandNamesL8).map(add_quality_score).map(GEE.add_simple_date)
    L8_distinct=L8.distinct(['simple_date'])
    L8 = ee.ImageCollection(L8_distinct.map(best_image_per_dateL8))
    L8 = GEE.resampleL7L8L9(L8, inputs)

    if inputs['Waterline']:
        wl_data=wl_data.merge(L8.map(wl_f).map(ensureProjLandsat))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L8.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L8.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L8.map(GEE.RGB).map(ensureProjLandsat))
#%%L5 ini

if useL5:
    L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_TOA")
    L5 = L5.filterDate(start,end).filter(ee.Filter.bounds(polygon_coll)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    if inputs['SpecSeasonalPeriod']:
        L5 = L5.filter(ee.Filter.calendarRange(inputs['SpecDates'][0],inputs['SpecDates'][1],'month'))
    L5 = L5.map(GEE.normalizeBandNamesL5).map(add_quality_score).map(GEE.add_simple_date)
    L5_distinct=L5.distinct(['simple_date'])
    L5 = ee.ImageCollection(L5_distinct.map(best_image_per_dateL5))
    L5 = GEE.resampleL5S2(L5,inputs)
    
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L5.map(wl_f).map(ensureProjLandsat))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L5.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L5.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L5.map(GEE.RGB).map(ensureProjLandsat))
#%%L7 ini

if useL7:
    L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA")
    L7 = L7.filterDate(start,end).filter(ee.Filter.bounds(polygon_coll)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    if inputs['SpecSeasonalPeriod']:
        L7 = L7.filter(ee.Filter.calendarRange(inputs['SpecDates'][0],inputs['SpecDates'][1],'month'))
    L7 = L7.map(GEE.normalizeBandNamesL7).map(add_quality_score).map(GEE.add_simple_date)
    L7_distinct=L7.distinct(['simple_date'])
    L7 = ee.ImageCollection(L7_distinct.map(best_image_per_dateL7))
    L7 = GEE.resampleL7L8L9(L7, inputs)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L7.map(wl_f).map(ensureProjLandsat))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L7.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L7.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L7.map(GEE.RGB).map(ensureProjLandsat))
#%%L9 ini

if useL9:
    L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA")
    L9 = L9.filterDate(start,end).filter(ee.Filter.bounds(polygon_coll)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    if inputs['SpecSeasonalPeriod']:
        L9 = L9.filter(ee.Filter.calendarRange(inputs['SpecDates'][0],inputs['SpecDates'][1],'month'))
    L9 = L9.map(GEE.normalizeBandNamesL9).map(add_quality_score).map(GEE.add_simple_date)
    L9_distinct=L9.distinct(['simple_date'])
    L9 = ee.ImageCollection(L9_distinct.map(best_image_per_dateL9))
    L9 = GEE.resampleL7L8L9(L9, inputs)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L9.map(wl_f).map(ensureProjLandsat))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L9.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L9.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L9.map(GEE.RGB).map(ensureProjLandsat))

#%% IMAGE SAVING ON THE LOCAL COMPUTER/SERVER
if inputs['Waterline']:
    sizeDataset = wl_data.size()
elif inputs['SandBar']:
    sizeDataset = sb_data.size()
elif inputs['Vegetation']:
    sizeDataset = sb_data.size()
    



#wl_data = wl_data.map(ensureProjLandsat)

if inputs['Waterline']:
    wl_list=wl_data.toList(sizeDataset)
    filepath_wl = os.path.join(os.getcwd(),inputs['WaterlineIndex'])
if inputs['SandBar']:
    sb_list = sb_data.toList(sizeDataset)
    filepath_sb = os.path.join(os.getcwd(),inputs['SandBarIndex'])
if inputs['Vegetation']:
    vg_list=vg_data.toList(sizeDataset)
    filepath_vg = os.path.join(os.getcwd(),inputs['VegetationIndex'])
if inputs['RGB']:
    rgb_list=rgb_data.toList(sizeDataset)
    filepath_rgb = os.path.join(os.getcwd(),'RGB')
length=sizeDataset.getInfo()

#CREATE A .CSV FILE TO CHECK THE COLLECTION OF IMAGES
#records = Tools.generate_attributes(wl_data.getInfo()['features'])
#Tools.upsert_satellite_csv(records)


for i in range(length):
    
    if inputs['Waterline']:
        img = ee.Image(wl_list.get(i))
        
        try:
            GEE.saveImage(img, i, length, filepath_wl, inputs['WaterlineIndex'], polygon_geom)
        except:
            continue
    
    if inputs['SandBar']:
        img = ee.Image(sb_list.get(i))
        GEE.saveImage(img, i, length, filepath_sb, inputs['SandBarIndex'], polygon_geom)
    
    if inputs['Vegetation']:
        img = ee.Image(vg_list.get(i))
        GEE.saveImage(img, i, length, filepath_sb, inputs['VegetationIndex'], polygon_geom)
        
    if inputs['RGB']:
        img = ee.Image(rgb_list.get(i))
        GEE.saveImage(img, i, length, filepath_rgb, ['red','green','blue'], polygon_geom)
        
#UPDATE THE .CSV FILE TO INFORM ON THE PROGRESS OF DOWNLOADED IMAGES
#records = Tools.generate_attributes(wl_data.getInfo()['features'])
#Tools.upsert_satellite_csv(records)
        
        





