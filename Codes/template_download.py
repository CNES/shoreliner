#---------------------------------------------------------------------------------------------
#--------------- Index maps generation and extraction from GEE -------------------------------
#--------------- Marcan Graffin (largely inspired by Kilian Vos' CoastSat tool) -------------
#---------------------------------------------------------------------------------------------
import ee
import time
from urllib.request import urlretrieve
import os
import shutil
import zipfile
from osgeo import gdal, osr
import pickle
import numpy as np
from shapely import geometry
import sys
sys.path.insert(0, r"../../../Codes/functions")
import yaml
from yaml.loader import SafeLoader
import Tools, GEE


global inputs
file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)

if not(inputs['Waterline'] or inputs['SandBar'] or inputs['Vegetataion']):
    sys.exit('No Feature to extract, please modify the config.yaml file')


useS2,useL9,useL8,useL7,useL5=Tools.satMissions(inputs['Missions'])

ee.Initialize(project=inputs['EE_ID'])

start,end = Tools.dates(inputs['Dates'])



filepath = './'
polygon=pickle.load(open('poly.json','rb'))
transects=pickle.load(open('transects.p','rb'))

GEE.folderCreation(inputs)

if inputs['Waterline']:
    wl_f = GEE.index_f(inputs['WaterlineIndex'])

if inputs['SandBar']:
    sb_f = GEE.index_f(inputs['SandBarIndex'])

if inputs['Vegetation']:
    vg_f = GEE.index_f(inputs['VegetationIndex'])

if not('RGB' in os.listdir()) and inputs['RGB']:
    os.mkdir('RGB')

global polygon_geom
polygon_geom=ee.Geometry.Polygon(polygon)
polygon_tmp = ee.FeatureCollection(ee.Feature(polygon_geom))

if inputs['Waterline']:
    wl_data=ee.ImageCollection([])
if inputs['SandBar']:
    sb_data=ee.ImageCollection([])
if inputs['Vegetation']:
    vg_data=ee.ImageCollection([])
if inputs['RGB']:
    rgb_data=ee.ImageCollection([])


#Image Collection generation and SCoWI computation 
#%%S2 ini

if useS2:
    S2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
    S2 = S2.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE",inputs['MaxCloudCover']))
    S2 = S2.map(GEE.normalizeBandNamesS2)
    S2 = GEE.resampleL5S2(S2,inputs)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(S2.map(wl_f))
    if inputs['SandBar']:
        sb_data=sb_data.merge(S2.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(S2.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(S2.map(GEE.RGB))
#%%L8 ini

if useL8:
    L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
    L8 = L8.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    L8 = L8.map(GEE.normalizeBandNamesL8)
    L8 = GEE.resampleL7L8L9(L8, inputs)

    if inputs['Waterline']:
        wl_data=wl_data.merge(L8.map(wl_f))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L8.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L8.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L8.map(GEE.RGB))
#%%L5 ini

if useL5:
    L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_TOA")
    L5 = L5.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    L5 = L5.map(GEE.normalizeBandNamesL5)
    L5 = GEE.resampleL5S2(L5,inputs)
    
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L5.map(wl_f))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L5.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L5.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L5.map(GEE.RGB))
#%%L7 ini

if useL7:
    L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_TOA")
    L7 = L7.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    L7 = L7.map(GEE.normalizeBandNamesL7)
    L7 = GEE.resampleL7L8L9(L7, inputs)
    
    if inputs['Waterline']:
        wl_data=wl_data.merge(L7.map(wl_f))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L7.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L7.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L7.map(GEE.RGB))
#%%L9 ini

if useL9:
    L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA")
    L9 = L9.filterDate(start,end).filter(ee.Filter.bounds(polygon_tmp)).filter(ee.Filter.lt('CLOUD_COVER',inputs['MaxCloudCover']))
    L9 = L9.map(GEE.normalizeBandNamesL9)
    L9 = GEE.resampleL7L8L9(L9, inputs)
  
    if inputs['Waterline']:
        wl_data=wl_data.merge(L9.map(wl_f))
    if inputs['SandBar']:
        sb_data=sb_data.merge(L9.map(sb_f))
    if inputs['Vegetation']:
        vg_data.merge(L9.map(vg_f))
    if inputs['RGB']:
        rgb_data = rgb_data.merge(L9.map(GEE.RGB))

#%% IMAGE SAVING ON THE LOCAL COMPUTER/SERVER
if inputs['Waterline']:
    sizeDataset = wl_data.size()
elif inputs['SandBar']:
    sizeDataset = sb_data.size()
elif inputs['Vegetation']:
    sizeDataset = sb_data.size()



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







for i in range(length):
    
    if inputs['Waterline']:
        img = ee.Image(wl_list.get(i))
        try:
            GEE.saveImage(img, i, filepath_wl, inputs['WaterlineIndex'], polygon_geom)
        except:
            continue
    
    if inputs['SandBar']:
        img = ee.Image(sb_list.get(i))
        GEE.saveImage(img, i, filepath_sb, inputs['SandBarIndex'], polygon_geom)
    
    if inputs['Vegetation']:
        img = ee.Image(vg_list.get(i))
        GEE.saveImage(img, i, filepath_sb, inputs['VegetationIndex'], polygon_geom)
        
    if inputs['RGB']:
        img = ee.Image(rgb_list.get(i))
        GEE.saveImage(img, i, filepath_rgb, ['red','green','blue'], polygon_geom)
        

        
