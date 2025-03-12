import numpy as np
import skimage.filters
from skimage.measure import find_contours
from skimage.transform import AffineTransform
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, r"../../../Codes")
from functions import Tools

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

def refinedOtsu(img,tag_idx,ax=[],val=256,ploting=False,id_image=0):
    img_val = img[np.logical_not(np.logical_or(np.isnan(img),np.isinf(img)))]
    if not checkHisto(img_val):
        return np.nan,np.nan,False
    
    hist, bins = np.histogram(img_val,val,density=True)

    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    
    try:
        t_otsu = skimage.filters.threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    
    
    """Minimum between 2 main peaks method"""
    peaks, properties = find_peaks(hist, height=(None,None))
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
    pdf_selected = hist[(bins>maxl)&(bins<maxr)]
    bins_selected = bins[(bins>maxl)&(bins<maxr)]
    if bins_selected.size==0:
        t_opti = t_otsu
    else:
        min_pdf = np.argmin(pdf_selected)
        t_opti = bins_selected[min_pdf]

    if ploting :
        plt.figure()
        ax=plt.subplot()
        ax.plot(bins,hist)
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
    img_val = img[np.logical_not(np.logical_or(np.isnan(img),np.isinf(img)))]
    if not checkHisto(img_val):
        return np.nan,False
    hist,bins=np.histogram(img_val,val,density=True)
    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    try:
        t_otsu = skimage.filters.threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    if ploting :
        ax.plot(bins,hist)
        maxy = ax.get_ylim()[1]
        ax.plot([t_otsu,t_otsu],[0,maxy],'k',linewidth=2)
        ax.grid()
        ax.set_xlabel('Index Value')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0,maxy])
    return t_otsu, True

def weightedPeaks(img,tag_idx,ax=[],val=256,ploting=False,id_image=0):
    img_val = img[np.logical_not(np.logical_or(np.isnan(img),np.isinf(img)))]
    if not checkHisto(img_val):
        return np.nan,np.nan,False
    
    hist, bins = np.histogram(img_val,val,density=True)

    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    t_otsu=0
    
    try:
        t_otsu = skimage.filters.threshold_otsu(img_val,val)
    except:
        print('fail_normal')
    
    
    """Minimum between 2 main peaks method"""
    peaks, properties = find_peaks(hist, height=(None,None))
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
    if tag_idx == 'SCoWI' or tag_idx == 'AWEIsh' or tag_idx == 'AWEIns':
        t_opti = maxr + 0.7*(maxl-maxr)
    else:
        t_opti = maxl + 0.7*(maxr-maxl)
    return t_otsu, t_opti, True

def MSV(img,tag_idx,ax=[],val=256,ploting=False,valcheckHisto=5,val_iter=100):
    img_val = img.flatten()[np.logical_not(np.isnan(img.flatten()))]
    img_val = img_val[np.logical_not(np.isinf(img_val))]
    if not checkHisto(img_val,valcheckHisto):
        return np.nan, np.nan, False
    
    hist, bins = np.histogram(img_val,val,density=True)

    bins = np.array([bins[i]+(bins[i+1]-bins[i])/2 for i in range(len(bins)-1)])
    
    try:
        thresh1 = skimage.filters.threshold_otsu(img_val)
        idmax2 = np.argmax(hist[bins>=thresh1])
        idmax1 = np.argmax(hist[bins<=thresh1])
        max2 = bins[bins>=thresh1][idmax2]
        max1 = bins[bins<=thresh1][idmax1]
    except:
        print('fail_normal')
    
    #return([thresh1,max1,max2])
    ranges = np.linspace(max1,max2,val_iter)
    length=[0]
    for value in ranges:
        contours = find_contours(img,value)
        l= np.max([len(i) for i in contours])
        length.append(l)


    ndelta = Tools.runnanmean(np.diff(np.diff(length)/length[1:]),5)
    ndelta = np.diff(length)/length[1:]
    #ndelta = np.diff(np.diff(length)/length[1:])
    ranges=np.flipud(ranges)
    ndelta = np.flipud(ndelta)

    return thresh1,ranges[1:][np.argmin(abs(ndelta))],True#,pd.Series(ndelta,index=ranges)