import pickle
import numpy as np
import os
from osgeo import gdal
import yaml
from yaml.loader import SafeLoader
import datetime
import sys
sys.path.insert(0, r"../../../Codes")
from functions import Tools,Segmentation,toSHP
import glob
import pandas as pd
import geopandas as gpd

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('../',file_inputs),'rb'),Loader = SafeLoader)
tag_idx = inputs['WaterlineIndex']
seg_method = inputs['Contouring']

path_index = os.path.join('./',tag_idx)
ploting=False
lag_corr=True
waterline = []
lags = {'S2':5,'L5':15,'L7':7.5,'L8':7.5,'L9':7.5}

polygon=pickle.load(open('poly.json','rb'))

lonEPSG=(polygon[0][0][0]+polygon[0][2][0])/2
latEPSG=(polygon[0][0][1]+polygon[0][2][1])/2
epsg_target = Tools.convert_wgs_to_utm(lonEPSG, latEPSG)

#%%Process
if inputs['NatureROI']=='Transects':
    transects = pickle.load(open('transects.p','rb'))
    
    
    for i in transects:
        tmp1 = [transects[i]['transect'][:][0][0],transects[i]['transect'][0][1]]
        tmp2 = [transects[i]['transect'][:][1][0],transects[i]['transect'][1][1]]
        transects[i]['transect'] = np.array([tmp1,tmp2])
    
    for i in transects:
      transects[i]['transect_proj'] = Tools.convert_epsg(transects[i]['transect'][:,::-1],4326,int(epsg_target))[:][:2]
    
    if inputs['TRshp']:
        toSHP.transectToSHPfile(transects, os.getcwd(), crs=int(epsg_target))

elif inputs['CreateTransects']:
    wl_composite, wl_noproj_composite = Segmentation.composite_image(inputs,epsg=epsg_target,len_min=1000)
    transects = Tools.transectsFromComposite(wl_composite,epsg=epsg_target,hspace=inputs['TransectSpacing'])

OTSU=[]
sat = []
dates = []
list_img = os.listdir(path_index)
list_img.sort()
c=-1
count=0
for i in list_img:
  if i[-3:] == 'tif':
        c+=1
        if c%10==0 and ploting:
            plting=True
        else:
            plting=False
        print(i)
        
        year = int(i[:4])
        month = int(i[4:6])
        day = int(i[6:8])
        hour = int(i[9:11])
        minute = int(i[11:13])
        
        date = datetime.datetime(year,month,day,hour,minute)
        
        
        path_tmp_index = os.path.join(path_index,i)
        img = gdal.Open(path_tmp_index)
        try:
            INDEX = img.ReadAsArray()
        except:
            continue
        lflat = len(INDEX.flatten())
        idx_nan=len(np.arange(lflat)[(INDEX==0.0).flatten()])
        INDEX[np.logical_or((INDEX==0.0),np.isinf(INDEX))]=np.nan
        rate = idx_nan/lflat
        if rate<0.5:
            gt = img.GetGeoTransform()
            if inputs['Contouring'] == 'RefOtsu' :
              t_otsuref,t_otsu,valid = Segmentation.refinedOtsu(INDEX, tag_idx, ploting=plting, id_image=i)
              if not valid:
                print('Noisy : failed the quality check')
                count+=1
                continue
            elif inputs['Contouring'] == 'WP' :
              t_otsuref,t_otsu,valid = Segmentation.weightedPeaks(INDEX, tag_idx, ploting=plting)
              if not valid:
                print('Noisy : failed the quality check')
                count+=1
                continue
            elif inputs['Contouring'] == 'MSV' :
              t_otsuref,t_otsu,valid = Segmentation.MSV(INDEX, tag_idx, ploting=plting)
              if not valid:
                print('Noisy : failed the quality check')
                count+=1
                continue
            elif inputs['Contouring'] == 'Otsu' :
              t_otsu, valid = Segmentation.otsu(INDEX, tag_idx, ploting=plting)
              if not valid:
                print('Noisy : failed the quality check')
                count+=1
                continue
            wl_tmp, wl_noproj_tmp = Segmentation.getWaterline(INDEX,t_otsu,gt,inputs,i=i)
            if lag_corr:
                wl_tmp[:,0] += lags[i[16:18]]
                wl_tmp[:,1] -= lags[i[16:18]]
            
            if inputs['WLshp']:
                toSHP.waterlineToSHPfile(wl_tmp, tag_idx, seg_method, i, os.getcwd(), crs=int(epsg_target))

            waterline.append(wl_tmp)
            OTSU.append(t_otsu)
            dates.append(date)
            sat.append(i[16:18])

if inputs['NatureROI']=='Transects':
    for key in transects:
        transects[key][tag_idx]=dict()
        transects[key][tag_idx]['raw']=dict()
        transects[key][tag_idx]['raw']['sat_dates'] = np.array(dates)
        transects[key][tag_idx]['raw']['sat_missions'] = np.array(sat)
        transects[key][tag_idx]['raw']['index_threshold'] = np.array(OTSU)
    
    transects = Segmentation.computeIntersection(waterline,transects,sat,inputs)
    
    for i in transects:
        transects[i][tag_idx]['raw'] = Tools.removeNaN(transects[i][tag_idx]['raw'],varname='SDW_'+inputs['WaterlineIndex'])
    
    pickle.dump(transects,open('transects_test.p','wb'))

if inputs['WLshp']:
    
    folder = os.path.join(os.getcwd(),"WL_"+inputs['WaterlineIndex'],inputs['Contouring'])
    files = glob.glob(os.path.join(folder,"*.shp"))
    
    gdfs = []
    
    for shp in files:
        gdf = gpd.read_file(shp)
        name = os.path.splitext(os.path.basename(shp))[0]
        # Add a field identifying the source shapefile
        gdf["source"] = name
        gdf["month"] = name[4:6]
        gdf["year"] = name[:4]
        gdfs.append(gdf)
    
    merged = gpd.GeoDataFrame(
        pd.concat(gdfs, ignore_index=True),
        crs=gdfs[0].crs
    )
    
    merged.to_file("Merged_WL_"+inputs['WaterlineIndex']+"_"+inputs['Contouring']+".shp")





