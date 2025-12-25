import numpy as np
from osgeo import osr
import skimage.filters
from scipy.ndimage import median_filter
from shapely.geometry import Polygon, LineString, Point
import geopandas as gpd
import os 
import pandas as pd
from astropy import timeseries
from scipy import signal, integrate, interpolate
import math
import pickle
import copy
import cv2
import matplotlib.pyplot as plt
from pyproj import Geod, Transformer

import csv

#%% transect construction

def createTransects(inputs,current_path):
    """
    Computes transects in a .p file corresponding to the chosen ROI
    The ROI is a polygon defined with points (lonA,latA),(lonB,latB)... with coords in espg:4326
     """
     
    """the distance between 2 transects"""
    inter_tr = inputs['Inter_transect']
    """the length of your transecs"""
    length_tr = inputs['Length_transect']
    """the coordinates of your polygon that defines your ROI"""
    ROI = Polygon(inputs['ROI_array'])
    """the GSHHS file"""
    gdf_file_wl = gpd.read_file(current_path+inputs['Shapefile_path'])
    """the approximate angle to describe 1km at the equator"""
    km_equator = 0.0089832
    
    
    poly_list = []
    for idx_poly in range(len(gdf_file_wl['id'])):
        if gdf_file_wl['area'][idx_poly] > 10:
            x,y = gdf_file_wl['geometry'][idx_poly].exterior.coords.xy
            poly_list.append(np.transpose(np.array([x.tolist(),y.tolist()])))
        

    """Extraction of points in the ROI"""
    wl_final_list = []
    poly_in_ROI = []
    for poly in poly_list:
        c=0
        wl_points_list = []
        for idx_p in range(len(poly)):
            lon,lat = poly[idx_p][0],poly[idx_p][1]
            if ROI.contains(Point((lon,lat))):
                c+=1
                if c==1:
                    poly_in_ROI.append(Polygon(poly))
                wl_points_list.append([lon,lat])
            else:
                continue
        if c>0:
            #wl_list = [[poly[idx_points[0]-1][0], poly[idx_points[0]-1][1]]] + wl_points_list + [[poly[idx_points[-1]+1][0], poly[idx_points[-1]+1][1]]]         
            wl_final_list.append(np.array(wl_points_list))
    
    
    """Waterline densification"""
    wl_final_dens = []
    for k in range(len(wl_final_list)):
        wl = wl_final_list[k]
        wl_dens = []
        for p in range(len(wl)-1):
            d = math.dist(wl[p],wl[p+1])
            nb_inter_pts = int(np.floor(d/inter_tr/km_equator)-1)
            if nb_inter_pts>0 and d<1:
                lon_list = np.linspace(wl[p][0],wl[p+1][0],nb_inter_pts+2)[:-1]
                lat_list = np.linspace(wl[p][1],wl[p+1][1],nb_inter_pts+2)[:-1]
                for i in range(len(lon_list)):
                    wl_dens.append([lon_list[i],lat_list[i]])
        wl_final_dens.append(np.array(wl_dens))
    
    
    """Construction of transects"""
    tr = []
    for l in range(len(wl_final_dens)):
        wl = wl_final_dens[l]
        tr_prov = []
        for j in range(1,len(wl)-1):
            lon,lat = wl[j][0],wl[j][1]
            coeff = abs(np.cos(lat*np.pi/180))
            a = -1*(wl[j+1][0]-wl[j-1][0])/(wl[j+1][1]-wl[j-1][1])
            u = length_tr*km_equator/(2*np.sqrt(1+a*a))
            if u!=0:
                tr_prov.append(np.array([[lon+u/coeff,lat+u*a],[lon-u/coeff,lat-u*a]]))
            else:
                if poly_in_ROI[l].contains(Point((lon,lat+length_tr*km_equator/2))):
                    tr_prov.append(np.array([[lon,lat+length_tr*km_equator/2],[lon,lat-length_tr*km_equator/2]]))
                else:
                    tr_prov.append(np.array([[lon,lat-length_tr*km_equator/2],[lon,lat+length_tr*km_equator/2]]))
        tr.append(np.array(tr_prov))
    
    """Check points order"""
    for p in range(len(tr)):
        try:
            if not poly_in_ROI[p].contains(Point(tr[p][0][0])):
                for q in range(len(tr[p])):
                    cop = tr[p][q][0].copy()
                    tr[p][q][0] = tr[p][q][1]
                    tr[p][q][1] = cop
        except:
            continue

    tr_final = []
    for p in range(len(tr)):
        for q in range(len(tr[p])):
            tr_final.append(tr[p][q])
    
    """Shapefile creation"""
    if inputs['TRshp']:
        trsc_lines = [LineString(t) for t in tr_final]
        gdftr = gpd.GeoDataFrame(geometry=trsc_lines,crs=4326)
        gdftr.to_file(os.path.join('./Projects/'+inputs['Project'],'TR.shp'))
    
    """Storage of data"""
    TR = dict()
    for i in range(len(tr_final)):
        name_tr = 'PF'+str(i)
        TR[name_tr] = dict()
        TR[name_tr]['transect'] = tr_final[i]
        TR[name_tr]['transect_proj'] = convert_epsg(tr_final[i], 4326, 3857)
    pickle.dump(TR,open(os.path.join('./DataExample',inputs['NewPathTransect']),'wb'))

def convert_wgs_to_utm(lon, lat):
    """Based on lat and lng, return best utm epsg-code"""
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
        return epsg_code
    epsg_code = '327' + utm_band
    return epsg_code


def polyFromTransects(transect,metric = False, d_ref = 0.05):
    polygon=[]
    classi_transect=[[]]
    lon_transect=[[]]
    lat_transect=[[]]
    first_key=list(transect.keys())[0]
    #print(transect[first_key])
    lon_ref=(transect[first_key]['transect'][0][0]+transect[first_key]['transect'][1][0])/2
    lat_ref=(transect[first_key]['transect'][0][1]+transect[first_key]['transect'][1][1])/2 
    for i in transect:
      lon_tmp=(transect[i]['transect'][0][0]+transect[i]['transect'][1][0])/2
      lat_tmp=(transect[i]['transect'][0][1]+transect[i]['transect'][1][1])/2
      d=np.sqrt((lon_tmp-lon_ref)**2+(lat_tmp-lat_ref)**2)
      if d<d_ref:
        lon_transect[-1].append(transect[i]['transect'][0][0])
        lon_transect[-1].append(transect[i]['transect'][1][0])
        lat_transect[-1].append(transect[i]['transect'][0][1])
        lat_transect[-1].append(transect[i]['transect'][1][1])
        classi_transect[-1].append(i)
      else:
        lon_transect.append([transect[i]['transect'][0][0],transect[i]['transect'][1][0]])
        lat_transect.append([transect[i]['transect'][0][1],transect[i]['transect'][1][1]])
        classi_transect.append([i])
        lon_ref=lon_tmp
        lat_ref=lat_tmp
      
    indx_boxes=[]
    rate_overl=0.01
    for j in range(len(classi_transect)):
       if not(lon_transect[j] == [] or lat_transect[j] == []):
         indx_boxes.append(j)
         lon1=max(lon_transect[j])+rate_overl*(max(lon_transect[j])-min(lon_transect[j]))
         lon2=min(lon_transect[j])-rate_overl*(max(lon_transect[j])-min(lon_transect[j]))
         lat1=max(lat_transect[j])+rate_overl*(max(lat_transect[j])-min(lat_transect[j]))
         lat2=min(lat_transect[j])-rate_overl*(max(lat_transect[j])-min(lat_transect[j]))
         deltalat = d_ref-(lat1-lat2)
         deltalon = d_ref-(lon1-lon2)*np.cos(lat1*2*3.1415/360)
         print(deltalat,deltalon)
         if deltalon>0:
            lon1 += deltalon/2/np.cos(lat1*2*3.1415/360)
            lon2 -= deltalon/2/np.cos(lat1*2*3.1415/360)
         if deltalat>0:
             lat1 += deltalat/2
             lat2 -= deltalat/2
         print('idx : '+str(j)+' | lat : [ '+str(lat1)+' , '+str(lat2)+' ]')
         tmp=[[[lon1,lat1],
               [lon1,lat2],
               [lon2,lat2],
               [lon2,lat1],
               [lon1,lat1]]]
         #if abs(lat1)<60 and abs(lat2)<60:
         polygon.append(tmp)
    return polygon, classi_transect


def convert_epsg(points, epsg_in, epsg_out):
    """
    Converts from one spatial reference to another using the epsg codes
    
    KV WRL 2018

    Arguments:
    -----------
    points: np.array or list of np.ndarray
        array with 2 columns (rows first and columns second)
    epsg_in: int
        epsg code of the spatial reference in which the input is
    epsg_out: int
        epsg code of the spatial reference in which the output will be            
                
    Returns:    
    -----------
    points_converted: np.array or list of np.array 
        converted coordinates from epsg_in to epsg_out
        
    """
    
    # define input and output spatial references
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(epsg_in)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(epsg_out)
    # create a coordinates transform
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    # if list of arrays
    if type(points) is list:
        points_converted = []
        # iterate over the list
        for i, arr in enumerate(points): 
            points_converted.append(np.array(coordTransform.TransformPoints(arr)))
    # if single array
    elif type(points) is np.ndarray:
        points_converted = np.array(coordTransform.TransformPoints(points))  
    else:
        raise Exception('invalid input type')

    return points_converted

def runnanmean(X,n):
    if isinstance(X,pd.core.series.Series):
        return pd.Series(np.convolve(X, np.ones(n)/n, mode='same'),list(X.index))
    else:
        tmp = X
        for i in range(len(X)):
            a = max(0,i-n//2)
            b = min(len(X),i+n//2)
            idx = np.arange(a,b)
            try:
                X[i] = np.nanmean(tmp[idx])
            except:
                X[i] = np.nan
                continue
        return X

def runmedian(X,n):
    if isinstance(X,pd.core.series.Series):
        return pd.Series(median_filter(X[X != np.nan], n),list(X.index))
    else:
        tmp = X
        for i in range(len(X)):
            a = max(0,i-n//2)
            b = min(len(X),i+n//2)
            idx = np.arange(a,b)
            try:
                X[i] = np.nanmedian(tmp[idx])
            except:
                X[i] = np.nan
                continue
    return X



#%% POST PROCESS
def removeNaN(TR,varname,refmin = 0):
    X = copy.deepcopy(TR[varname])
    nans = np.isnan(X)
    belowzero = X<refmin #means can't have a waterline more landward than the transect origin
    cond = np.logical_or(nans,belowzero)
    idx_nan = np.arange(len(X))[cond]
    for i in TR:
        TR[i] = np.delete(TR[i],idx_nan)
    return(TR)


def Mode(TR,inputs,n=0.65,valid=True):
    
    newTR = copy.deepcopy(TR)
    X = copy.deepcopy(newTR['SDW_'+inputs['WaterlineIndex']])
    if valid:
        tmpX = X - np.nanmax([0,np.nanmedian(X)])
    else:
        tmpX = X - np.nanmedian(X)
    
    threshold = abs(skimage.filters.threshold_otsu(tmpX))
    # print(threshold)
    # print(n*np.std(tmpX))
    idx_otsu = []
    if threshold > n*np.std(tmpX):
        for i in range(len(tmpX)):
            if tmpX[i]>threshold or tmpX[i]<-threshold: #we remove data both sides, but most certainly outliers will be in the tmpX[i]<-threshold  part.
                idx_otsu.append(i)
    for i in newTR:
        newTR[i] = np.delete(newTR[i],idx_otsu)
    return(newTR)



def IQR(TR, varname,val1 = 0.25, val2 = 0.75, ratio=1.5):
    
    """
    data[PROFILE]['satellite'] = IQR(data[PROFILE]['satellite'])
    Clean a SDW timeseries by removing data deviating too much
    from the distribution, clean also the corresponding dates and others...
    """
    
    newTR = copy.deepcopy(TR)
    X = copy.deepcopy(newTR[varname])
    Q1 = np.quantile(X,val1)
    Q3 = np.quantile(X,val2)
    IQR = Q3-Q1
    valmax = Q3 + ratio*IQR
    valmin = Q1 - ratio*IQR
    
    idx_iqr=[]
    for i in range(len(X)):
        if X[i]>valmax or X[i]<valmin:
            idx_iqr.append(i)
            
    for i in newTR:
        try:
            newTR[i] = np.delete(newTR[i],idx_iqr)
        except:
            newTR[i] = np.delete(np.array(newTR[i],dtype='object'),idx_iqr)
            continue
    
    return newTR

def findNearestTimes(data_in,dates_in,dates_out,delta_t):
    
    """
    new_data = findNearestTimes(old_data,old_dates,new_dates)
    Returns the data
    """
    
    data_out = []
    tmp_in=np.array([dates_in[i].timestamp() for i in range(len(dates_in))])
    tmp_out=np.array([dates_out[i].timestamp() for i in range(len(dates_out))])

    for i in tmp_out:
        tmpdiff = abs(tmp_in-i)
        indx= np.argmin(tmpdiff)
        if min(tmpdiff)/3600/24<delta_t:
            data_out.append(data_in[indx])
        else:
            data_out.append(np.nan)
    return data_out



#%% Water level corrections

def rangeSlopes(minSlope,maxSlope,deltaSlope=0.0025):
    return np.arange(max(minSlope,deltaSlope),maxSlope+deltaSlope,deltaSlope)

def wlCorrect(Xi,wl,slopes,zref=0.0):
    'apply waterlevel correction with a range of slopes'
    Xall = []
    for i in range(len(slopes)):
        # apply tidal correction
        tmpX=Xi.copy()
        tide_correction = (-zref+wl)/slopes[i]
        tmpX += tide_correction
        Xall.append(tmpX)
    return Xall

def wlPeak(dates,wl):
    'find the high frequency peak in the tidal time-series'
    # create frequency grid
    t = np.array([_.timestamp() for _ in dates]).astype('float64')
    #days_in_year = 365.2425
    seconds_in_day = 24*3600
    time_step = 8*seconds_in_day
    freqs = getFrequencyGrid(t,time_step,50)
    # compute power spectrum
    ps_tide,_,_ = powerSpectrum(t,wl,freqs,[])
    # find peaks in spectrum
    idx_peaks,_ = signal.find_peaks(ps_tide, height=0)
    y_peaks = _['peak_heights']
    idx_peaks = idx_peaks[np.flipud(np.argsort(y_peaks))]
    # find the strongest peak at the high frequency (defined by freqs_cutoff[1])
    idx_max = idx_peaks[freqs[idx_peaks] > 1./(seconds_in_day*30)][0]
    # compute the frequencies around the max peak with some buffer (defined by buffer_coeff)
    freqs_max = [freqs[idx_max] - 1e-8, freqs[idx_max] + 1e-8]
    # make a plot of the spectrum
    return freqs_max

def getFrequencyGrid(time,time_step,n0):
    'define frequency grid for Lomb-Scargle transform'
    T = np.max(time) - np.min(time)
    fmin = 1/T
    fmax = 1/(2*time_step) # Niquist criterium
    df = 1/(n0*T)
    N = np.ceil((fmax - fmin)/df).astype(int)
    freqs = fmin + df * np.arange(N)
    return freqs

def powerSpectrum(t,y,freqs,idx_cut):
    'compute power spectrum and integrate'
    model = timeseries.LombScargle(t, y, dy=None, fit_mean=True, center_data=True, nterms=1, normalization='psd')
    ps = model.power(freqs)
    # integrate the entire power spectrum
    E = integrate.simpson(ps, x=freqs, even='avg')
    if len(idx_cut) == 0:
        idx_cut = np.ones(freqs.size).astype(bool)
    # integrate only frequencies above cut-off
    Ec = integrate.simpson(ps[idx_cut], x=freqs[idx_cut], even='avg')
    return ps, E, Ec

def integratePowerSpectrum(dates,Xall,freqMax):
    'integrate power spectrum at the frequency band of peak tidal signal'
    t = np.array([_.timestamp() for _ in dates]).astype('float64')
    seconds_in_day = 24*3600
    time_step = 8*seconds_in_day
    freqs = getFrequencyGrid(t,time_step,50)    
    beach_slopes = rangeSlopes(0.0025, 0.2, 0.0025)
    # integrate power spectrum
    idx_interval = np.logical_and(freqs >= freqMax[0], freqs <= freqMax[1]) 
    E = np.zeros(beach_slopes.size)
    for i in range(len(Xall)):
        ps, _, _ = powerSpectrum(t,Xall[i],freqs,[])
        E[i] = integrate.simpson(ps[idx_interval], x=freqs[idx_interval], even='avg')
    # calculate confidence interval
    delta = 0.0001
    prc = 0.05
    f = interpolate.interp1d(beach_slopes, E, kind='linear')
    beach_slopes_interp = rangeSlopes(0.0025,0.2-delta,delta)
    E_interp = f(beach_slopes_interp)
    # find values below minimum + 5%
    slopes_min = beach_slopes_interp[np.where(E_interp <= np.min(E)*(1+prc))[0]]
    if len(slopes_min) > 1:
        ci = [slopes_min[0],slopes_min[-1]]
    else:
        ci = [beach_slopes[np.argmin(E)],beach_slopes[np.argmin(E)]]
    
    return beach_slopes[np.argmin(E)], ci

def drawConceptualGraph(inputs,pathsave):
    
    fig,ax = plt.subplots(figsize=(20,12))
    img = np.zeros((250, 500, 3), dtype="uint8")
    img[:] = (255,255,255)
    for i in range(4):
        for j in range(2):
            x1 = 5+130*i
            y1 = 5+130*j
            x2 = 105+130*i
            y2 = 105+130*j
            # print([(x1,y1),(x2,y2)])
            cv2.rectangle(img,(x1,y1),(x2,y2),(0,0,0),(1))
    cv2.arrowedLine(img,(105,55),(134,55),(0,0,0),(1))
    cv2.arrowedLine(img,(235,55),(264,55),(0,0,0),(1))
    cv2.arrowedLine(img,(365,55),(394,55),(0,0,0),(1))
    cv2.arrowedLine(img,(445,105),(445,134),(0,0,0),(1))
    cv2.arrowedLine(img,(135,185),(106,185),(0,0,0),(1))
    cv2.arrowedLine(img,(265,185),(236,185),(0,0,0),(1))
    cv2.arrowedLine(img,(395,185),(366,185),(0,0,0),(1))
    plt.imshow(img)
    plt.text(10,22,'STUDY CASE\nDEFINITION',fontsize=15,fontweight='bold')
    plt.text(140,15,'IMAGE COLLECTION',fontsize=15,fontweight='bold')
    plt.text(270,15,'PRE-PROCESSING',fontsize=15,fontweight='bold')
    plt.text(400,15,'BAND COMBINATION',fontsize=15,fontweight='bold')
    plt.text(400,145,'THRESHOLDING',fontsize=15,fontweight='bold')
    plt.text(270,145,'CONTOURING',fontsize=15,fontweight='bold')
    plt.text(140,152,'PROJECTION OVER\nTRANSECTS',fontsize=15,fontweight='bold')
    plt.text(10,145,'POST-PROCESSING',fontsize=15,fontweight='bold')
    
    ##STUDY
    
    plt.text(10,45,'localisation : '+inputs['Project'],fontsize=13)
    plt.text(10,60,'subROI size : '+str(inputs['SizeROI'])+'°',fontsize=13)
    plt.text(10,75,'start : '+inputs['Dates'][0],fontsize=13)
    plt.text(10,90,'end : '+inputs['Dates'][1],fontsize=13)
    
    ##IMAGE COLLECTION
    missions = ''
    for i in inputs['Missions']:
        missions += i+', '
    missions = missions[:-2]
    plt.text(140,45,'missions : '+missions,fontsize=13)
    plt.text(140,60,'cloud filter : '+str(inputs['MaxCloudCover'])+'%',fontsize=13)
    
    ##PRE-PROCESSING
    prepro = inputs['Preprocessing']
    if 'Pansharpening' in prepro:
        plt.text(270,45,'pansharpening : yes',fontsize=13)
    else:
        plt.text(270,45,'pansharpening : no',fontsize=13)
    if 'Bicubic' in prepro:
        plt.text(270,75,'resampling : bicubic\n(S2 : SWIR1, SWIR2 : 10m\nLandsat : all bands : 15 m)',fontsize=13)
    elif 'Bilinear' in prepro :
        plt.text(270,75,'resampling : bilinear\n(S2 : SWIR1, SWIR2 : 10m\nLandsat : all bands : 15 m)',fontsize=13)
    
    ##BAND COMBINATION
    plt.text(400,45,'index : '+inputs['WaterlineIndex'],fontsize=13)
    if inputs['PeakNoiseFilter']:
        plt.text(400,60,'peak/noise filter : yes (ratio : '+str(inputs['PeakNoiseFilterRatio'])+')',fontsize=13)
    else:
        plt.text(400,60,'peak/noise filter : no',fontsize=13)
    
    ##THRESHOLDING
    method = inputs['Contouring']
    if method == 'RefOtsu':
        method='Local Minimum'
    if method == 'WP':
        method = 'Weighted Peaks'
    if method == 'MSV':
        method = 'Minimum Shoreline\nVariability'
    plt.text(400,175,'method : '+method,fontsize=13)
    
    ##CONTOURING
    plt.text(270,175,'method : Marching Squares',fontsize=13)
    ##POST-PROCESSING
    postpro = inputs['Postprocessing']
    if 'IQR' in postpro:
        plt.text(10,175,'outlier removal : IQR',fontsize=13)
    else:
        plt.text(10,175,'outlier removal : no',fontsize=13)
    if 'TideCorrection' in postpro:
        plt.text(10,205,'tide correction : yes',fontsize=13)
    else:
        plt.text(10,205,'tide correction : no',fontsize=13)
    # if 'Resampling' in postpro:
    #     plt.text(10,190,'resampling : '+inputs['Postprocessing']['Resampling'],fontsize=13)
    # else:
    #     plt.text(10,190,'resampling : no',fontsize=13)
    # if 'Interpolation' in postpro:
    #     plt.text(10,205,'interpolation : '+inputs['Postprocessing']['Interpolation'],fontsize=13)
    # else:
    #     plt.text(10,205,'interpolation : no',fontsize=13)
    
    fig.patch.set_visible(False)
    ax.axis('off')
    
    plt.savefig(os.path.join(pathsave,'SDWConceptualGraph.png'))

def upsert_satellite_csv(records, csv_path='img_collections.csv',UNIQUE_KEY = 'id'):
    """
    records: list of dicts with new metadata
    csv_path: path to CSV catalog
    """

    existing_rows = {}
    existing_fields = set()

    # 1. Load existing CSV if it exists
    if os.path.exists(csv_path):
        with open(csv_path, mode="r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            existing_fields = set(reader.fieldnames)

            for row in reader:
                existing_rows[row[UNIQUE_KEY]] = row

    # 2. Collect fields from new records
    new_fields = set()
    for record in records:
        new_fields.update(record.keys())

    # 3. Union of all fields (stable order)
    all_fields = list(existing_fields | new_fields)

    # Ensure UNIQUE_KEY is first column
    if UNIQUE_KEY in all_fields:
        all_fields.remove(UNIQUE_KEY)
    all_fields.insert(0, UNIQUE_KEY)

    # 4. Update or insert records
    for record in records:
        key = record[UNIQUE_KEY]

        if key in existing_rows:
            existing_rows[key].update(record)
        else:
            existing_rows[key] = {k: record.get(k, "") for k in all_fields}

    # 5. Normalize rows (fill missing fields)
    normalized_rows = []
    for row in existing_rows.values():
        normalized_rows.append({k: row.get(k, "") for k in all_fields})

    # 6. Write back CSV
    with open(csv_path, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(normalized_rows)

def generate_attributes(metadata):
    
    img_index_downloaded=0
    out = []
    
    for i in range(len(metadata)):
        tmp = metadata[i]
        name = tmp['properties']['date'][:8]+'T'+tmp['properties']['date'][8:]+'_'+tmp['properties']['satellite']
        
        if name+'.tif' in os.listdir('./'+tmp['bands'][0]['id']):
            img_index_downloaded=1
            
        tmp_out = {
            'id':name,
            'cloud_cover':tmp['properties']['cloud_cover'],
            tmp['bands'][0]['id']:img_index_downloaded}
        
        out.append(tmp_out)
        
    return out
