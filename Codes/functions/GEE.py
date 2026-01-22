import ee
import os
from urllib.request import urlretrieve
import shutil
import zipfile
from osgeo import gdal
global inputs
import datetime as dt


#%% General

def index_f(name):
    if name == 'SCoWI':
        return(SCoWI)
    elif name == 'MNDWI':
        return(MNDWI)
    elif name == 'NDWI':
        return(NDWI)
    elif name == 'AWEIsh':
        return(AWEIsh)
    elif name == 'AWEIns':
        return(AWEIns)
    elif name == 'NDVI':
        return(NDVI)
    elif name == 'SRWI':
        return(SRWI)
    elif name == 'SBI':
        return(SBI)
    else:
        raise('Error : Unknown index. Index  known by default : SCoWI, MNDWI, NDWI, NDVI, AWEIns and AWEIsh')

def resampleL5S2(satlist,inputs,val=32629):
    process = inputs['Preprocessing']
    if 'Bicubic' in process:
        satlist = satlist.map(bicubicResampleAll)
    if 'Bilinear' in process:
        satlist = satlist.map(bilinearResampleAll)
    return satlist

def resampleL7L8L9(satlist,inputs):
    process = inputs['Preprocessing']
    if 'Pansharpening' in process:
        satlist = satlist.map(pansharpening)
        if 'Bicubic' in process:
            satlist = satlist.map(bicubicResamplePan)
        if 'Bilinear' in process:
            satlist = satlist.map(bilinearResamplePan)
    else:
        if 'Bicubic' in process:
            satlist = satlist.map(bicubicResampleAll)
        if 'Bilinear' in process:
            satlist = satlist.map(bilinearResampleAll)
    return satlist

def preprocess(satlist,inputs):
    process = inputs['Preprocessing']
    if 'Pansharpening' in process:
        satlist = satlist.map(pansharpening)
    return satlist

def dates(dates):
    start = dates[0]
    end = dates[1]
    eeStart = ee.Date.fromYMD(int(start[:4]),int(start[5:7]),int(start[8:10]))
    eeEnd = ee.Date.fromYMD(int(end[:4]),int(end[5:7]),int(end[8:10]))
    return eeStart,eeEnd

def add_simple_date(image):
    date = ee.String(image.get('date')).slice(0, 8)
    return image.set('simple_date', date)

def filterNonFullImages(coll,inputs):
    pix_max = ee.Number(coll.aggregate_max('valid_pixels'))
    return coll.filter(ee.Filter.gte('valid_pixels', pix_max.multiply(inputs['FillingThreshold'])))

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

#%% Preprocess
def bicubicResampleAll(image):
    return image.resample('bicubic').reproject(
        **{
      'crs':image.select('blue').projection(),
      'scale': image.get('res')}).copyProperties(image,['res','date','satellite'])

def bilinearResampleAll(image):
    return image.resample('bilinear').reproject(
        **{
      'crs':image.select('blue').projection(),
      'scale': image.get('res')}).copyProperties(image,['res','date','satellite'])

def bicubicResamplePan(image):
    red = image.select(['red'],['red'])
    blue = image.select(['blue'],['blue'])
    green =image.select(['green'],['green'])
    nir = image.select(['nir'],['nir']).resample('bicubic')
    swir1 = image.select(['swir1'],['swir1']).resample('bicubic')
    swir2 = image.select(['swir2'],['swir2']).resample('bicubic')
    panchromatic = image.select(['panchromatic'],['panchromatic'])
    newImage = ee.Image.cat([red,green,blue,nir,swir1,swir2,panchromatic])
    return newImage.reproject(
        **{
      'crs':image.select('blue').projection(),
      'scale': image.get('res')}).copyProperties(image,['res','date','satellite'])

def bilinearResamplePan(image):
    red = image.select(['red'],['red'])
    blue = image.select(['blue'],['blue'])
    green =image.select(['green'],['green'])
    nir = image.select(['nir'],['nir']).resample('bilinear')
    swir1 = image.select(['swir1'],['swir1']).resample('bilinear')
    swir2 = image.select(['swir2'],['swir2']).resample('bilinear')
    panchromatic = image.select(['panchromatic'],['panchromatic'])
    newImage = ee.Image.cat([red,green,blue,nir,swir1,swir2,panchromatic])
    return newImage.reproject(
        **{
      'crs':image.select('blue').projection(),
      'scale': image.get('res')}).copyProperties(image,['res','date','satellite'])

def pansharpening(image):
    hsv = image.select(['red', 'green', 'blue']).reproject(
        **{
      'crs':image.select('blue').projection(),
      'scale': image.get('res')}).rgbToHsv();
    
    sharpened = ee.Image.cat([
            hsv.select('hue'), hsv.select('saturation'), image.select('panchromatic')
            ]).hsvToRgb();

    blue = sharpened.select(['blue'],['blue'])
    green = sharpened.select(['green'],['green'])
    red = sharpened.select(['red'],['red'])
    
    nir = image.select(['nir'],['nir'])
    swir1 = image.select(['swir1'],['swir1'])
    swir2 = image.select(['swir2'],['swir2'])
    panchromatic = image.select(['panchromatic'],['panchromatic'])
    
    newImage = ee.Image.cat([red,green,blue,nir,swir1,swir2,panchromatic])
    return newImage.copyProperties(image,['date','satellite'])

def normalizeBandNamesS2(image):
    newImage = image.select(['B2','B3','B4','B8','B11','B12'],['blue','green','red','nir','swir1','swir2'])
    return newImage.set({'res':10,'satellite':'S2', 'date':image.date().format('YYYYMMddHHmmss'),'cloud_cover':image.get("CLOUDY_PIXEL_PERCENTAGE")}).float()
    
def normalizeBandNamesL9(image):
    newImage = image.select(['B2','B3','B4','B5','B6','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'res':15,'satellite':'L9', 'date':image.date().format('YYYYMMddHHmmss'),'cloud_cover':image.get('CLOUD_COVER')}).float()

def normalizeBandNamesL8(image):
    newImage = image.select(['B2','B3','B4','B5','B6','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'res':15,'satellite':'L8', 'date':image.date().format('YYYYMMddHHmmss'),'cloud_cover':image.get('CLOUD_COVER')}).float()

def normalizeBandNamesL7(image):
    newImage = image.select(['B1','B2','B3','B4','B5','B7','B8'],['blue','green','red','nir','swir1','swir2','panchromatic'])
    return newImage.set({'res':15,'satellite':'L7', 'date':image.date().format('YYYYMMddHHmmss'),'cloud_cover':image.get('CLOUD_COVER')}).float()

def normalizeBandNamesL5(image):
    newImage = image.select(['B1','B2','B3','B4','B5','B7'],['blue','green','red','nir','swir1','swir2'])
    return newImage.set({'res':15,'satellite':'L5', 'date':image.date().format('YYYYMMddHHmmss'),'cloud_cover':image.get('CLOUD_COVER')}).float()

def maskLandsat(image):
  # Bits 3 and 5 are cloud shadow and cloud, respectively.
  cloudShadowBitMask = 1 << 3;
  cloudsBitMask = 1 << 5;
  # Get the pixel QA band.
  qa = image.select('QA_PIXEL');
  #  Both flags should be set to zero, indicating clear conditions.
  maskShadow = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
  maskCloud = qa.bitwiseAnd(cloudsBitMask).eq(0);
  #Return the masked image, scaled to reflectance, without the QA bands.
  return image.updateMask(maskShadow).updateMask(maskCloud).copyProperties(image, ['SCENE_CENTER_TIME'])

def get_s2_sr_cld_col(aoi, start_date, end_date,cloud_cover_treshold = 50):
    # Import and filter S2 SR.
    s2_sr_col = (ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_cover_treshold)))

    # Import and filter s2cloudless.
    s2_cloudless_col = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        .filterBounds(aoi)
        .filterDate(start_date, end_date))

    # Join the filtered s2cloudless collection to the SR collection by the 'system:index' property.
    return s2_sr_col.combine(s2_cloudless_col)


def add_cloud_maskS2(img):
    # Get s2cloudless image, subset the probability band.
    cld_prb = img.select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(75).reproject(**{'crs': 'EPSG:3857', 'scale': 20})
    
    no_cloud_img = is_cloud.select('probability').Not()
    # Add the cloud probability layer and cloud mask as image bands.
    return img.select('B.*').updateMask(no_cloud_img).addBands(ee.Image([cld_prb, is_cloud]))
#%% Indexes
def SCoWI(image):
    return image.expression('(BLUE + 2 * GREEN - 2 * NIR - 0.75 * SWIR1 - 0.5 * SWIR2)',#
    {
      'BLUE': image.select('blue'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
      'SWIR2': image.select('swir2')
    }).rename('SCoWI').copyProperties(image,['date','satellite','cloud_cover'])

def MNDWI(image):
    return image.expression('(SWIR1-GREEN)/(SWIR1+GREEN)',#
    {
      'GREEN': image.select('green'),
      'SWIR1': image.select('swir1')
    }).rename('MNDWI').copyProperties(image,['date','satellite','cloud_cover'])

def NDVI(image):
    return image.expression('(NIR-RED)/(NIR+GREEN)',#
    {
      'RED': image.select('red'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir')
    }).rename('NDVI').copyProperties(image,['date','satellite','cloud_cover'])

def SBI(image):
    return image.expression('BLUE+2*GREEN-1.75*NIR',#
    {
      'BLUE': image.select('blue'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir')
    }).rename('SBI').copyProperties(image,['date','satellite','cloud_cover'])


def NDWI(image):
    return image.expression('(NIR-GREEN)/(NIR+GREEN)',#
    {
      'GREEN': image.select('green'),
      'NIR': image.select('nir')
    }).rename('NDWI').copyProperties(image,['date','satellite','cloud_cover'])


def AWEIsh(image):
    return image.expression('(BLUE + 2.5 * GREEN - 1.5 * NIR - 1.5 * SWIR1 - 0.25 * SWIR2)',#
    {
      'BLUE': image.select('blue'),
      'GREEN': image.select('green'),
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
      'SWIR2': image.select('swir2')
    }).rename('AWEIsh').copyProperties(image,['date','satellite','cloud_cover'])

def AWEIns(image):
    return image.expression('(4 * GREEN - 4 * SWIR1 - 0.25 * NIR - 2.75 * SWIR2)',#
    {
      'GREEN': image.select('green'),
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
      'SWIR2': image.select('swir2')
    }).rename('AWEIns').copyProperties(image,['date','satellite','cloud_cover'])

def SRWI(image):
   return image.expression('( GREEN + BLUE - NIR - SWIR1)/(GREEN + BLUE + NIR + SWIR1)',#
    {
      'GREEN': image.select('green'),
      'NIR': image.select('nir'),
      'SWIR1': image.select('swir1'),
      'BLUE': image.select('blue')
    }).rename('SRWI').copyProperties(image,['date','satellite','cloud_cover'])

def RGB(image):
        #image = image.resample('bicubic')
        
        red = image.select(['red'])
        blue = image.select(['blue'])
        green = image.select(['green'])
        
        new_image = ee.Image.cat([red,green,blue])
        return new_image.copyProperties(image,['date','satellite','cloud_cover'])
    
def ALLBANDS(image):
    
    blue = image.select(['blue'])
    green = image.select(['green'])
    red = image.select(['red'])
    nir = image.select(['nir'])
    swir1 = image.select(['swir1'])
    swir2 = image.select(['swir2'])
    
    new_image = ee.Image.cat([blue,green,red,nir,swir1,swir2])
    return new_image.copyProperties(image,['date','satellite','cloud_cover'])
    
#%% Data/Folder 

def folderCreation(inputs):
    
    if inputs['Waterline']:
        try:
            os.mkdir(inputs['WaterlineIndex'])
        except:
            pass
    if inputs['SandBar']:
        try:
            os.mkdir(inputs['SandBarIndex'])
        except:
            pass
    if inputs['Vegetation']:
        try:
            os.mkdir(inputs['VegetationIndex'])
        except:
            pass

def saveImage(imggee, i, length, filepath, index, poly):
    
    img = imggee.getInfo()
    name=img['properties']['date'][:8]+'T'+img['properties']['date'][8:]+'_'+img['properties']['satellite']
    print(name)
    cond1=not(True in [name in j for j in os.listdir(filepath)])
    if (cond1):
        print('Downloading... ('+str(i+1)+'/'+str(length)+')')
        url = ee.data.makeDownloadUrl(ee.data.getDownloadId({
            'image': imggee,
            'bands': index,
            'region': poly,
            'name': name
            }))
    
        try:
            local_zip, headers = urlretrieve(url)
            # move zipfile from temp folder to data folder
            dest_file = os.path.join(filepath, 'imagezip')
            shutil.move(local_zip,dest_file)
            # unzip file
            with zipfile.ZipFile(dest_file) as local_zipfile:
                for fn in local_zipfile.namelist():
                    local_zipfile.extract(fn, filepath)
            # filepath + filename to single bands
                fn_tifs = [os.path.join(filepath,_) for _ in local_zipfile.namelist()]
        # stack bands into single .tif
            outds = gdal.BuildVRT(os.path.join(filepath,'stacked'+str(i)+'.vrt'), fn_tifs, separate=True)
            outds = gdal.Translate(os.path.join(filepath,name+'.tif'), outds)
            for fn in fn_tifs: os.remove(fn)
        # delete .vrt file
            os.remove(os.path.join(filepath,'stacked'+str(i)+'.vrt'))
        # delete zipfile
            os.remove(dest_file)
        # delete data.tif.aux (not sure why this is created)
            if os.path.exists(os.path.join(filepath,'data.tif.aux')):
                os.remove(os.path.join(filepath,'data.tif.aux'))
        except:
            pass
    # else:
    #     name=name+'_1'
    #     cond2 = not(True in [name in j for j in os.listdir(filepath)])
    #     if (cond2):
    #         print('Downloading... ('+str(i+1)+'/'+str(length)+')')
    #         url = ee.data.makeDownloadUrl(ee.data.getDownloadId({
    #             'image': imggee,
    #             'bands': index,
    #             'region': poly,
    #             'name': name
    #             }))
        
    #         try:
    #             local_zip, headers = urlretrieve(url)
    #             # move zipfile from temp folder to data folder
    #             dest_file = os.path.join(filepath, 'imagezip')
    #             shutil.move(local_zip,dest_file)
    #             # unzip file
    #             with zipfile.ZipFile(dest_file) as local_zipfile:
    #                 for fn in local_zipfile.namelist():
    #                     local_zipfile.extract(fn, filepath)
    #             # filepath + filename to single bands
    #                 fn_tifs = [os.path.join(filepath,_) for _ in local_zipfile.namelist()]
    #         # stack bands into single .tif
    #             outds = gdal.BuildVRT(os.path.join(filepath,'stacked'+str(i)+'.vrt'), fn_tifs, separate=True)
    #             outds = gdal.Translate(os.path.join(filepath,name+'.tif'), outds)
    #         # delete single-band files
    #             for fn in fn_tifs: os.remove(fn)
    #         # delete .vrt file
    #             os.remove(os.path.join(filepath,'stacked'+str(i)+'.vrt'))
    #         # delete zipfile
    #             os.remove(dest_file)
    #         # delete data.tif.aux (not sure why this is created)
    #             if os.path.exists(os.path.join(filepath,'data.tif.aux')):
    #                 os.remove(os.path.join(filepath,'data.tif.aux'))
    #         except:
    #             pass

