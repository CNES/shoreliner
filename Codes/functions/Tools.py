import numpy as np
import ee
from osgeo import osr
import skimage.filters
from skimage.measure import find_contours
from skimage.transform import AffineTransform
from scipy.ndimage import median_filter
from shapely import geometry
from shapely.geometry import Polygon, LineString, Point
import geopandas as gpd
import os 

import pandas as pd
from scipy.stats import linregress,pearsonr
from astropy import timeseries
from scipy import signal, integrate, interpolate
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import math
import pickle
import copy

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
#%%
def dates(dates):
    start = dates[0]
    end = dates[1]
    eeStart = ee.Date.fromYMD(int(start[:4]),int(start[5:7]),int(start[8:10]))
    eeEnd = ee.Date.fromYMD(int(end[:4]),int(end[5:7]),int(end[8:10]))
    return eeStart,eeEnd

def satMissions(missions):
    S2,L9,L8,L7,L5 = False,False,False,False,False
    if 'S2' in missions:
        S2=True
    if 'L9' in missions:
        L9=True
    if 'L8' in missions:
        L8=True
    if 'L7' in missions:
        L7=True
    if 'L5' in missions:
        L5=True
    return S2,L9,L8,L7,L5

def polyFromTransects(transect, d_ref = 0.05):
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
         delta = abs(deltalat-deltalon)
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
         if abs(lat1)<60 and abs(lat2)<60:
             polygon.append(tmp)
    return polygon, classi_transect


def norm(band):
    band_min, band_max = band.min(), band.max()
    return ((band - band_min)/(band_max - band_min))


def refinedOtsu(img,tag_idx,ax=[],val=256,ploting=False,id_image=0):
    img_val = img[np.logical_not(np.isnan(img))]
    if not checkHisto(img_val):
        return np.nan,np.nan,False
    
    hist, bins = np.histogram(img_val,val,density=True)
    pdf=runmedian(hist,4)

    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    
    try:
        t_otsu = skimage.filters.threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    
    
    """Minimum between 2 main peaks method"""
    peaks, properties = find_peaks(pdf, height=(None,None))
    maxs_histo = bins[peaks]

    res_max_left = maxs_histo[maxs_histo < t_otsu]
    res_max_right = maxs_histo[maxs_histo >= t_otsu]
    height_left = properties['peak_heights'][:len(res_max_left)]
    height_right = properties['peak_heights'][len(res_max_left):]

    if res_max_left.size==0:
        maxl = t_otsu 
    else:
        idx=np.argmax(height_left)
        maxl = res_max_left[idx]

    if res_max_right.size==0:
        maxr = t_otsu
    else:
        idx=np.argmax(height_right)
        maxr = res_max_right[idx]

    # select the part of the histogram between 2 peaks around the threshold
    pdf_selected = pdf[(bins>maxl)&(bins<maxr)]
    bins_selected = bins[(bins>maxl)&(bins<maxr)]
    if bins_selected.size==0:
        t_opti = t_otsu
    else:
        min_pdf = np.argmin(pdf_selected)
        t_opti = bins_selected[min_pdf]
    

    """Minimum at the bottom left of the water peak for SCoWI/AWEIsh/AWEI"""
    # peaks, properties = find_peaks(pdf, height=(None,None), distance=10)
    # maxs_histo = bins[peaks]

    # res_max_left = maxs_histo[maxs_histo < t_otsu]
    # height_left = properties['peak_heights'][:len(res_max_left)]
    # res_max_right = maxs_histo[maxs_histo >= t_otsu]
    # height_right = properties['peak_heights'][len(res_max_left):]

    # if tag_idx in ['SCoWI','AWEIsh','AWEI']:
    #     if res_max_right.size==0:
    #         maxl = maxr = t_otsu
    #     else:
    #         idx=np.argmax(height_right)
    #         maxr = res_max_right[idx]
    #         if idx!=0:
    #             maxl = res_max_right[idx-1]
    #         else:
    #             maxl = t_otsu
                
    # else:
    #     if res_max_left.size==0:
    #         maxl = t_otsu 
    #     else:
    #         idx=np.argmax(height_left)
    #         maxl = res_max_left[idx]

    #     if res_max_right.size==0:
    #         maxr = t_otsu
    #     else:
    #         idx=np.argmax(height_right)
    #         maxr = res_max_right[idx]

    # # select the part of the histogram between 2 peaks around the threshold
    # pdf_selected = pdf[(bins>maxl)&(bins<maxr)]
    # bins_selected = bins[(bins>maxl)&(bins<maxr)]
    # if bins_selected.size==0:
    #     t_opti = t_otsu
    # else:
    #     min_pdf = np.argmin(pdf_selected)
    #     t_opti = bins_selected[min_pdf]



    if ploting :
        plt.figure()
        ax=plt.subplot()
        ax.plot(bins,pdf)
        maxy = ax.get_ylim()[1]
        # ax.plot([bins[peaks[maxi1]],bins[peaks[maxi1]]],[0,maxy],'k',linewidth=2)
        # ax.plot([bins[peaks[maxi2]],bins[peaks[maxi2]]],[0,maxy],'k',linewidth=2)
        ax.plot([t_otsu,t_otsu],[0,maxy],'r',linewidth=2)
        ax.plot([t_opti,t_opti],[0,maxy],'g',linewidth=2)
        ax.grid()
        ax.set_xlabel('Index Value')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0,maxy])
        plt.title(id_image)
        
    return t_otsu,t_opti,True



def otsu(img,tag_idx,ax=[],val=256,ploting=False):
    img_val = img[np.logical_not(np.isinf(img))]
    hist=np.histogram(img_val,val,density=True)
    pdf=runmedian(hist[0],5)
    bins=hist[1]
    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    try:
        t_otsu = skimage.filters.threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    if ploting :
        ax.plot(bins,pdf)
        maxy = ax.get_ylim()[1]
        ax.plot([t_otsu,t_otsu],[0,maxy],'k',linewidth=2)
        ax.grid()
        ax.set_xlabel('Index Value')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0,maxy])
    return t_otsu


def checkHisto(X,ploting=False,val=5):
    """
    It compute the mean of the two main peaks and the median of the histogram,
    if the ratio is smaller than a certain threshold, the histogram isn't valid
    """
    
    hist=np.histogram(X,bins=256,density=True)
    pdf = hist[0]
    bins = hist[1]
    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    thresh1 = skimage.filters.threshold_otsu(X)
    
    try:
        max1 = max(pdf[bins>=thresh1])
        max2 = max(pdf[bins<=thresh1])
    except:
        return False
    meanmax = (max1+max2)/2
    idx = np.arange(np.argmax(pdf[bins<=thresh1])+1,len(pdf[bins<=thresh1])+np.argmax(pdf[bins>=thresh1])-1)
    intermed = np.median(pdf[idx])

    if ploting:
        plt.figure()
        #print('Ratio = '+str(meanmax/intermed))
        plt.plot(bins[bins>=thresh1],pdf[bins>=thresh1],'red')
        plt.plot(bins[bins<=thresh1],pdf[bins<=thresh1],'green')
        plt.plot([bins[idx[0]],bins[idx[-1]]],[meanmax,meanmax])
        plt.plot([bins[idx[0]],bins[idx[-1]]],[intermed,intermed])
        
    if meanmax/intermed>val:
        return True
    else:
        return False


def getWaterline(img,threshold,georef,transects,date,inputs,i='    ',ax=[],MIN_LENGTH_SL=0,ploting=False):
    contours=find_contours(img,threshold)
    
    contours_out = [] #non_projected contours
    for j in contours:
        for k in j:
            contours_out.append(k)
    
    aff_mat=np.array([[georef[1], georef[2], georef[0]],
                            [georef[4], georef[5], georef[3]],
                            [0, 0, 1]])
    
    tform = AffineTransform(aff_mat)
    if type(contours) is list:
        points_converted = []
        points_regular = []
        # iterate over the list
        for l, arr in enumerate(contours):
            
            tmp = arr[:,[1,0]]
            if len(tmp)>10:
                points_regular.append(tmp)
                points_converted.append(tform(tmp))
    # elif type(contours) is np.ndarray:
    #     tmp = contours[:,[1,0]]
    #     points_converted = tform(tmp)
        
        
        
        
    if ploting:
        plt.figure()
        L,l =  img.shape[0],img.shape[1]
        x = np.arange(georef[0],georef[0]+(L+1)*georef[1],georef[1])
        y = np.arange(georef[3],georef[3]+(l+1)*georef[5],georef[5])
        X,Y = np.meshgrid(x,y)
        quantiles=np.nanquantile(img.flatten(),[0.05,0.95])
        try:
            plt.pcolormesh(X,Y,img,cmap='Greys_r',vmin=quantiles[0],vmax=quantiles[1])
            for j in transects:
                lon = transects[j]['transect_proj'][:,0]
                lat = transects[j]['transect_proj'][:,1]
                plt.plot(lon,lat,'r')
                if 'situ' in transects[j]:
                    d = ((lon[0]-lon[1])**2+(lat[0]-lat[1])**2)**.5
                    tsitu = transects[j]['situ']['dates']
                    timg = np.array([date])
                    idx = findNearestTimes(range(len(tsitu)), tsitu, timg, delta_t=inputs['DeltaTimeSituSat'])[0]
                    if not(np.isnan(idx)):
                        X = transects[j]['situ']['chainage'][idx]
                        Z = np.array(transects[j]['situ']['elevation'][idx])
                        Xsitu = getCrossPos(X,Z,inputs['MSLOffset'])
                        posx = lon[0] + Xsitu/d*(lon[1] - lon[0])
                        posy = lat[0] + Xsitu/d*(lat[1] - lat[0])
                        plt.plot([posx],[posy],'x',color='orange',markersize=10)
            for j in points_converted:
                plt.plot(j[:,0],j[:,1],'.b')
        except:
            plt.pcolormesh(img,cmap='Greys_r',vmin=quantiles[0],vmax=quantiles[1])
            tform_inv = AffineTransform(np.linalg.inv(aff_mat))
            for j in transects:
                lon = tform_inv(transects[j]['transect_proj'])[:,0]
                lat = tform_inv(transects[j]['transect_proj'])[:,1]
                plt.plot(lon,lat,'r')
            for j in points_regular:
                plt.plot(j[:,0],j[:,1],'.b')

        plt.title(i)
        
    
    
    contours_coord=points_converted
    contours_long = [] #projected contours
    
    for l, wl in enumerate(contours_coord):
        coords = [(wl[k,0], wl[k,1]) for k in range(len(wl))]
        # line = geometry.LineString(coords) # shapely LineString structure
        dist = np.sqrt((coords[0][0]-coords[-1][0])**2+(coords[0][1]-coords[-1][1])**2)
        if len(coords)>30 and dist!=0:
            contours_long.append(wl)
            
    
    x_points = np.array([])
    y_points = np.array([])
    for k in range(len(contours_long)):
        x_points = np.append(x_points,contours_long[k][:,0])
        y_points = np.append(y_points,contours_long[k][:,1])
    shoreline = np.transpose(np.array([x_points,y_points]))
    contours_out = np.array(contours_out)
    return shoreline, contours_out

def computeIntersection(shorelines, transects, sat_id, inputs):
    """
    Computes the intersection between the 2D shorelines and the shore-normal.
    transects. Adapted from CoastSat compute_intersections function (Vos et al., 2019)
    """
    if 'Pansharpening' in inputs['Preprocessing']:
        along_dist = {'S2':5,'L5':15,'L7':7.5,'L8':7.5,'L9':7.5}
    else:
        along_dist = {'S2':5,'L5':15,'L7':15,'L8':15,'L9':15}
    # loop through shorelines and compute the median intersection    
    intersections = np.zeros((len(shorelines),len(transects)))
    for i in range(len(shorelines)):
        #print('shoreline n°'+str(i+1)+' over '+str(len(shorelines)))
        sl = shorelines[i]
        max_along = along_dist[sat_id[i]]
        for j,key in enumerate(list(transects.keys())): 
            
            # compute rotation matrix
            X0 = transects[key]['transect_proj'][0,0]
            Y0 = transects[key]['transect_proj'][0,1]
            X1 = transects[key]['transect_proj'][1,0]
            Y1 = transects[key]['transect_proj'][1,1]
            temp = np.array(transects[key]['transect_proj'][-1,:]) - np.array(transects[key]['transect_proj'][0,:])
            phi = np.arctan2(temp[1], temp[0])
            Mrot = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
    
            # calculate point to line distance between shoreline points and the transect
            p1 = np.array([X0,Y0])
            p2 = np.array([X1,Y1])
            d_line = np.abs(np.cross(p2-p1,sl-p1)/np.linalg.norm(p2-p1))
            # calculate the distance between shoreline points and the origin of the transect
            d_origin = np.array([np.linalg.norm(sl[k,:] - p1) for k in range(len(sl))])
            # find the shoreline points that are close to the transects and to the origin
            # the distance to the origin is hard-coded here to 1 km 
            idx_along_dist = d_line <= max_along
            idx_cross_dist = d_origin <= ((X1-X0)**2+(Y1-Y0)**2)**.5
            # find the shoreline points that are in the direction of the transect (within 90 degrees)
            temp_sl = sl - p1
            scal = np.dot(temp_sl,p2-p1)
            idx_angle = scal>=0
            # combine the transects that are close in distance and close in orientation
            idx_close = np.where(np.logical_and(idx_along_dist,idx_cross_dist,idx_angle))[0]     
            # idx_close = np.where(idx_dist)[0]
            
            # in case there are no shoreline points close to the transect 
            if len(idx_close) == 0:
                intersections[i,j] = np.nan
            else:
                # change of base to shore-normal coordinate system
                xy_close = np.array([sl[idx_close,0],sl[idx_close,1]]) - np.tile(np.array([[X0],
                                   [Y0]]), (1,len(sl[idx_close])))
                xy_rot = np.matmul(Mrot, xy_close)
                # compute the max of the intersections along the transect
                intersections[i,j] = np.nanmax(xy_rot[0,:])
                        
    for j,key in enumerate(list(transects.keys())):
        transects[key][inputs['WaterlineIndex']]['raw']['SDW_'+inputs['WaterlineIndex']] = intersections[:,j]
        # transects[key][inputs['WaterlineIndex']]['raw']['SDW_'+inputs['WaterlineIndex']] = intersections[:,j]
    return transects

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

#%% QUICK CHECKS

def quickCheck(profile,inputs,wl=0.,var = 'SDW_', index='SCoWI', ploting = True, t_d=np.array([])):
    
    """
    use it as quickCheck(transects[PROFILEYOUWANNASEE] 
                         for visual comparison situ/sat
    """
    z = wl + inputs['MSLOffset']
    Xsat = profile['satellite'][var+index]
    Xref = profile['satellite']['SDW_'+index]
    dXswash = profile['satellite']['dXswash']
    dXsetup = profile['satellite']['dXsetup']
    nwl = profile['satellite']['dZtide']
    Msat = profile['satellite']['sat']
    tsat = profile['satellite']['dates']
    tsitu = []
    # X,Z = [],[]
    X = copy.deepcopy(profile['situ']['chainage'])
    Z = copy.deepcopy(profile['situ']['elevation'])
    Xsitu = []
    idx_error = []
    
    
    for i in range(len(X)):
        try:
            Xsitu.append(getCrossPos(X[i],Z[i], z))
        except:
            idx_error.append(i)
            continue
    for i in profile['situ']:
        profile['situ'][i] = np.delete(np.array(profile['situ'][i],dtype='object'),idx_error)
    profile['situ']['Xsitu_'+var+index] = Xsitu
    # print(len(Xsitu))
    # print(len(tsitu)) 
    profile['situ'] = IQR(profile['situ'],inputs,'Xsitu_'+var+index,ratio=1.5)
    Xsitu = profile['situ']['Xsitu_'+var+index]
    tsitu = profile['situ']['dates']
    

    nXsat = findNearestTimes(Xsat, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
    nXref = findNearestTimes(Xref, tsat, tsitu ,delta_t=inputs['DeltaTimeSituSat'])
    ndXswash = findNearestTimes(dXswash, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
    ndXsetup = findNearestTimes(dXsetup, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
    nwl = findNearestTimes(nwl, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
    nMsat = findNearestTimes(Msat, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])

    idx_nan = np.arange(len(nXsat))[np.isnan(nXsat)]
    nXsat = np.delete(nXsat,idx_nan)
    nXref = np.delete(nXref,idx_nan)
    ndXswash = np.delete(ndXswash,idx_nan)
    ndXsetup = np.delete(ndXsetup,idx_nan)
    nwl = np.delete(nwl,idx_nan)
    nMsat = np.delete(nMsat,idx_nan)
    Xsitu = np.delete(Xsitu,idx_nan)
    tsitu = np.delete(tsitu,idx_nan)
    
    #print(len(nXsat),len(Xsitu))
    
    out = getStats(nXsat,Xsitu)
    maxitmp = max([max(nXsat),max(Xsitu)])
    minitmp = min([min(nXsat),min(Xsitu)])
    mini = minitmp - 0.05*(maxitmp-minitmp)
    maxi = maxitmp + 0.05*(maxitmp-minitmp)
    if ploting:
        fig,ax = plt.subplot_mosaic("AAB",figsize=((18,6)))
        ax['A'].plot(tsitu,nXsat,'.k',label='satellite_corr')
        ax['A'].plot(tsitu,nXref,'.g',label='satellite_no_corr')
        ax['A'].plot(tsitu,Xsitu,'.r',label='situ')
        ax['A'].fill_between(tsitu,nXsat-ndXswash,nXsat+ndXswash)
        ax['A'].set_xlabel('Dates')
        ax['A'].set_ylabel('Cross-shore position (m)')
        ax['A'].legend()
        
        Xsitu_S2,nXsat_S2,Xsitu_L8,nXsat_L8,Xsitu_L7,nXsat_L7,Xsitu_L5,nXsat_L5,tsitu_L7,tsitu_L5 = [],[],[],[],[],[],[],[],[],[]
        for m in range(len(nMsat)):
            if nMsat[m] == 'S2':
                Xsitu_S2.append(Xsitu[m])
                nXsat_S2.append(nXsat[m])
            if nMsat[m] == 'L8':
                Xsitu_L8.append(Xsitu[m])
                nXsat_L8.append(nXsat[m])
            if nMsat[m] == 'L7':
                Xsitu_L7.append(Xsitu[m])
                nXsat_L7.append(nXsat[m])
                tsitu_L7.append(tsitu[m])
            if nMsat[m] == 'L5':
                Xsitu_L5.append(Xsitu[m])
                nXsat_L5.append(nXsat[m])
                tsitu_L5.append(tsitu[m])
        ax['B'].plot(Xsitu_S2, nXsat_S2, '.', label='S2')
        ax['B'].plot(Xsitu_L8, nXsat_L8, '+', label='L8')
        ax['B'].plot(Xsitu_L7, nXsat_L7, 'x', label='L7')
        ax['B'].plot(Xsitu_L5, nXsat_L5, '*', label='L5')
        
        # tsitu_L7 = [7-abs(t.month-7) for t in tsitu_L7]
        # tsitu_L7 = [t.year for t in tsitu_L7]
        # ax['B'].scatter(Xsitu_L7, nXsat_L7, c=tsitu_L7, cmap='coolwarm')
        # tsitu_L5 = [7-abs(t.month-7) for t in tsitu_L5]
        # tsitu_L5 = [t.year for t in tsitu_L5]
        # ax['B'].scatter(Xsitu_L5, nXsat_L5, c=tsitu_L5, cmap='coolwarm')
        
        # ax['B'].scatter(Xsitu, nXsat, c=nwl, cmap='coolwarm', vmin=-1.5, vmax=1.5)
        ax['B'].set_xlabel('Situ (m)')
        ax['B'].set_ylabel('Satellite (m)')
        ax['B'].set_xlim([mini,maxi])
        ax['B'].set_ylim([mini,maxi])
        ax['B'].plot([mini,maxi],[mini,maxi],'--b')
        # ax['B'].legend(loc=4)
        
        tRMSE = 'RMSE = '+str(round(out['RMSE'],2))+' m'
        tR2 = 'R² = '+str(round(out['R2'],2))
        tbias = 'bias = '+str(round(out['bias'],2))+' m'
        tstd = r'$\sigma$ = '+str(round(out['STD'],2))+' m'
        tdatasize = 'datasize = '+str(out['datasize'])+' points'
        allt = tR2+'\n'+tRMSE+'\n'+tbias+'\n'+tstd+'\n'+tdatasize
        plt.text(minitmp,maxitmp - 0.15*(maxitmp-minitmp),allt)
    return(out,nXsat,Xsitu,tsitu,nwl,nXref)


def StatsAndPlots(TR, i, inputs, steps):
    
    tsteps = ''
    for s in steps:
        tsteps = tsteps+s[-5:]
        Xsat = copy.deepcopy(TR[i][s]['SDW_'+inputs['WaterlineIndex']])
        tsat = copy.deepcopy(TR[i][s]['dates'])
        # if (s[-5:] in ['tcorr','wcorr']) and inputs['TideCorrectionWithSitu']:
        #     Xsitu = copy.deepcopy(TR[i]['situ_tcorr_used']['Xsitu'])
        #     tsitu = copy.deepcopy(TR[i]['situ_tcorr_used']['tsitu'])
        # else:
        Xsitu = copy.deepcopy(TR[i]['situ_used']['Xsitu'])
        tsitu = copy.deepcopy(TR[i]['situ_used']['tsitu'])
            
        
        nXsat,Xsitu,tsitu = findNearestTimes2(Xsat,tsat,Xsitu,tsitu,delta_t=inputs['DeltaTimeSituSat'])
        
        idx_nan = np.arange(len(Xsitu))[np.isnan(Xsitu)]
        nXsat = np.delete(nXsat,idx_nan)
        Xsitu = np.delete(Xsitu,idx_nan)
        tsitu = np.delete(tsitu,idx_nan)

        
        if len(nXsat)>10 and len(nXsat)==len(Xsitu): #we want points for stats
            if not 'stats' in TR[i]:
                TR[i]['stats'] = dict()
            out = TR[i]['stats'][s[-5:]] = getStats(nXsat, Xsitu)
            if s[-5:] == 'wcorr':
                dXswash = copy.deepcopy(TR[i][s]['swash_uncert'])
                dXswash = findNearestTimes(dXswash, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
                match = np.logical_and(np.array(nXsat)-np.array(dXswash)<Xsitu,Xsitu<np.array(nXsat)+np.array(dXswash))
                c=0
                for m in match:
                    if m == True:
                        c+=1 
                TR[i]['stats'][s[-5:]]['PartInSwashUncert'] = c/len(Xsitu)
            
            if inputs['Ploting'] and i in inputs['TransectsToPlot']:
                maxitmp = max([max(nXsat),max(Xsitu)])
                minitmp = min([min(nXsat),min(Xsitu)])
                mini = minitmp - 0.05*(maxitmp-minitmp)
                maxi = maxitmp + 0.05*(maxitmp-minitmp)
                
                fig,ax = plt.subplot_mosaic("AAB",figsize=((18,6)))
                ax['A'].set_title(inputs['Project']+' -- transect '+i+' -- '+inputs['WaterlineIndex']+' -- '+tsteps)
                ax['A'].plot(tsitu,nXsat,'k',label='satellite')
                ax['A'].plot(tsitu,Xsitu,'r',label='situ')
                if s[-5:] == 'wcorr':
                    ax['A'].fill_between(tsitu,np.array(nXsat)-np.array(dXswash),np.array(nXsat)+np.array(dXswash))
                ax['A'].set_xlabel('Dates')
                ax['A'].set_ylabel('Cross-shore position (m)')
                ax['A'].legend()
                
                if inputs['SatelliteSymbol']:
                    Msat = copy.deepcopy(TR[i][s]['sat'])
                    Msat = findNearestTimes(Msat, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
                
                    Xsitu_S2,nXsat_S2,Xsitu_L8,nXsat_L8,Xsitu_L7,nXsat_L7,Xsitu_L5,nXsat_L5 = [],[],[],[],[],[],[],[]
                    for m in range(len(Msat)):
                        if Msat[m] == 'S2':
                            Xsitu_S2.append(Xsitu[m])
                            nXsat_S2.append(nXsat[m])
                        if Msat[m] == 'L8':
                            Xsitu_L8.append(Xsitu[m])
                            nXsat_L8.append(nXsat[m])
                        if Msat[m] == 'L7':
                            Xsitu_L7.append(Xsitu[m])
                            nXsat_L7.append(nXsat[m])
                        if Msat[m] == 'L5':
                            Xsitu_L5.append(Xsitu[m])
                            nXsat_L5.append(nXsat[m])
                    ax['B'].plot(Xsitu_S2, nXsat_S2, '.', label='S2')
                    ax['B'].plot(Xsitu_L8, nXsat_L8, '+', label='L8')
                    ax['B'].plot(Xsitu_L7, nXsat_L7, 'x', label='L7')
                    ax['B'].plot(Xsitu_L5, nXsat_L5, '*', label='L5')
                    ax['B'].legend(loc=4)
                
                else:
                    if s[-5:] in ['tcorr','wcorr']:
                        wl = copy.deepcopy(TR[i][inputs['WaterlineIndex']+'_satellite_tcorr']['waterlevel'])
                        tsat = copy.deepcopy(TR[i][inputs['WaterlineIndex']+'_satellite_tcorr']['dates'])
                        wl = findNearestTimes(wl, tsat, tsitu, delta_t=inputs['DeltaTimeSituSat'])
                        ax['B'].scatter(Xsitu, nXsat, c=wl, cmap='coolwarm')
                    else:
                        ax['B'].plot(Xsitu, nXsat, '.g')
                
                ax['B'].set_xlabel('Situ (m)')
                ax['B'].set_ylabel('Satellite (m)')
                ax['B'].set_xlim([mini,maxi])
                ax['B'].set_ylim([mini,maxi])
                ax['B'].plot([mini,maxi],[mini,maxi],'--k')
                tRMSE = 'RMSE = '+str(round(out['RMSE'],2))+' m'
                tR2 = 'R² = '+str(round(out['R2'],2))
                tbias = 'bias = '+str(round(out['bias'],2))+' m'
                tstd = r'$\sigma$ = '+str(round(out['STD'],2))+' m'
                tdatasize = 'datasize = '+str(out['datasize'])
                if s[-5:] == 'wcorr':
                    tPartInSwashUncert = 'PartInSwashUncert = '+str(round(c/len(tsitu)*100,2))+' %'
                    allt = tR2+'\n'+tRMSE+'\n'+tbias+'\n'+tstd+'\n'+tPartInSwashUncert+'\n'+tdatasize
                else:
                    allt = tR2+'\n'+tRMSE+'\n'+tbias+'\n'+tstd+'\n'+tdatasize
                plt.text(minitmp,maxitmp - 0.2*(maxitmp-minitmp),allt)
                
            tsteps = tsteps+' + '
        else:
            continue



def transectToSHPfile(transects, current_path, crs=3857):
    """
    Create a .shp file of the transects to visualize it on QGIS
    """
    trsc_lines = [LineString(transects[tr]['transect_proj']) for tr in transects]
    gdftr = gpd.GeoDataFrame(geometry=trsc_lines,crs=crs)
    gdftr.to_file(os.path.join(current_path,'TRSC.shp'))

def waterlineToSHPfile(wl_tmp, tag_idx, image_name, current_path, crs=3857):
    """
    Create a .shp file of the waterline to visualize it on QGIS
    """
    if not('WL_'+tag_idx in os.listdir(current_path)):
        os.mkdir('WL_'+tag_idx)
    points = [Point(p) for p in wl_tmp]
    gdf = gpd.GeoDataFrame(geometry=points,crs=crs)
    gdf.to_file(os.path.join(os.path.join(current_path,'WL_'+tag_idx),image_name[:-4]+'_WL.shp'))


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


def findNearestTimes2(data_1, data_1_bis, data_1_ter, t_1, data_2, data_2_bis, t_2, delta_t):
    """
    It return datas and dates of each list. Each date of t_out_1 corresponds to a unique date in t_out_2
    so the process is a bijection between the two datasets
    """
    t_out_1=[]
    t_out_2=[]
    data_out_1=[]
    data_out_1_bis=[]
    data_out_1_ter=[]
    data_out_2=[]
    data_out_2_bis=[]
    tmp_1=np.array([t_1[i].timestamp() for i in range(len(t_1))])
    tmp_2=np.array([t_2[i].timestamp() for i in range(len(t_2))])
    tuplelist1=[]
    tuplelist2=[]
    
    for t1 in tmp_1:
        tmpdiff1 = abs(tmp_2-t1)
        if min(tmpdiff1)/3600/24<delta_t:
            tuplelist1.append((tmp_1.tolist().index(t1),tmpdiff1.tolist().index(min(tmpdiff1))))
    for t2 in tmp_2:
        tmpdiff2 = abs(tmp_1-t2)
        if min(tmpdiff2)/3600/24<delta_t:
            tuplelist2.append((tmpdiff2.tolist().index(min(tmpdiff2)),tmp_2.tolist().index(t2)))
    for tup in tuplelist1:
        if tup in tuplelist2:
            data_out_1.append(data_1[tup[0]])
            data_out_1_bis.append(data_1_bis[tup[0]])
            data_out_1_ter.append(data_1_ter[tup[0]])
            data_out_2.append(data_2[tup[1]])
            data_out_2_bis.append(data_2_bis[tup[1]])
            t_out_1.append(t_1[tup[0]])
            t_out_2.append(t_2[tup[1]])
    return data_out_1, data_out_1_bis, data_out_1_ter, t_out_1, data_out_2, data_out_2_bis, t_out_2


def getStats(X,Y):
    
    """
    STATS = getStats(X,Y)
    returns the common validations stats when comparing 2 lists of the same length
    """
    X = np.array(X)
    Y = np.array(Y)
    
    out = dict()
    R = out['R'] = np.corrcoef(X,Y)[0][1]
    out['R2'] = R**2
    MSE = np.mean((X - Y) ** 2)
    out['RMSE'] = np.sqrt(MSE)
    bias = out['bias'] = np.mean(X - Y)
    out['STD'] = (MSE - bias**2)**.5
    out['datasize'] = len(X)
    if len(X)>1:
        out['pvalue'] = pearsonr(X,Y)[1]
    return out

def slopeFromProfile(X,Z,zref=0,window=1.0):
    
    """
    slope = slopeFromProfile(X,Z)
    returns the slope around the elevation Z=zref
    """
    if Z[0]<Z[-1]:
        Z = Z[::-1]
        X = X[::-1]
    if not Z[0]>zref>Z[-1]:
        return np.nan
    Z=np.array(Z)
    X=np.array(X)
    idx = np.logical_and(Z>zref-window,Z<zref+window)
    if sum(bool(x) for x in idx)<2:
        z = [Z[np.where(Z>zref)[0][-1]],Z[np.where(Z<zref)[0][0]]]
        x = [X[np.where(Z>zref)[0][-1]],X[np.where(Z<zref)[0][0]]]
    else:
        newidx=[]
        i=list(idx).index(True)
        while idx[i]==True:
            newidx.append(i)
            i+=1
            if i==len(idx):
                break
        x=X[newidx]
        z=Z[newidx]
    slope = linregress(x,z)[0]
    return -slope



def slopeFromProfile2(X,Z,z1,z2):
    
    x1 = getCrossPos(X, Z, z1)
    x2 = getCrossPos(X, Z, z2)
    return (z2-z1)/(x1-x2)
    

def getCrossPos(X,Z,z):
    
    """
    x = getCrossPos(X,Z,z)
    returns the cross-shore (x) position of the intersection between
    the elevation z and the profile(X,Z)
    """
    if Z[0]<Z[-1]:
        Z = Z[::-1]
        X = X[::-1]
    if X[0]<=0:
        offset=X[0]
        for j in range(len(X)):
            X[j]-=offset

    X = np.array(X)
    Z = np.array(Z)
    idx = np.argmin(abs(Z-z))
    if Z[idx]<z and not idx==0:
        z1 = abs(Z[idx-1]-z)
        z2 = abs(Z[idx]-z)
        dX = X[idx]-X[idx-1]
        return(X[idx-1]+z1/(z1+z2)*dX)
    elif Z[idx]>z and not idx==len(Z)-1:
        z1 = abs(Z[idx]-z)
        z2 = abs(Z[idx+1]-z)
        dX = X[idx+1]-X[idx]
        return(X[idx]+z1/(z1+z2)*dX)
    elif Z[idx]==z:
        return X[idx]
    else:
        return np.nan
    
    


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

