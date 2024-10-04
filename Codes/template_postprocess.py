import numpy as np
import pickle
import sys
sys.path.insert(0, r"../../Codes")
from functions import Tools
import os
import yaml
from yaml.loader import SafeLoader
import copy
import xarray as xr


file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('./',file_inputs),'rb'),Loader = SafeLoader)
tag_idx = inputs['WaterlineIndex']



if 'transects_post.p' in os.listdir():
    transects = pickle.load(open('transects_post.p','rb'))
    for i in sorted(os.listdir()):
        if i[:5]=='poly_':
            pathtmp = os.path.join(os.getcwd(),i, 'transects.p')
            tmptransect = pickle.load(open(os.path.join(pathtmp),'rb'))
            for j in tmptransect:
                if 'raw' in tmptransect[j][tag_idx]:
                    transects[j][tag_idx] = tmptransect[j][tag_idx]
    
else:
    transects = dict()
    for i in sorted(os.listdir()):
        if i[:5]=='poly_':
            pathtmp = os.path.join(os.getcwd(),i, 'transects.p')
            tmptransect = pickle.load(open(os.path.join(pathtmp),'rb'))
            for j in tmptransect:
                if 'raw' in tmptransect[j][tag_idx]:
                    transects[j] = tmptransect[j]


"""Raw"""
for i in transects:

    tsat = copy.deepcopy(transects[i][tag_idx]['raw']['sat_dates'])
    Msat = copy.deepcopy(transects[i][tag_idx]['raw']['sat_missions'])
    for d in range(1,len(tsat)):
        if tsat[d].month==tsat[d-1].month and tsat[d].day==tsat[d-1].day and Msat[d]==Msat[d-1]:
            transects[i][tag_idx]['raw']['SDW_'+tag_idx][d] = np.nan
    transects[i][tag_idx]['raw'] = Tools.removeNaN(transects[i][tag_idx]['raw'], inputs, varname='SDW_'+tag_idx)
    


"""IQR"""
if inputs['IQR']:
    for i in transects:
        transects[i][tag_idx]['IQR'] = Tools.IQR(transects[i][tag_idx]['raw'],inputs,varname='SDW_'+tag_idx)



"""TideCorrection"""
if inputs['TideCorrection']:

    for i in transects:
        transects[i][tag_idx]['tcorr'] = copy.deepcopy(transects[i][tag_idx]['IQR'])    


    print('Computing tide predictions')
    carrier_file = 'latLonTime_maree.nc'
    time = []
    lon = []
    lat = []
    n_t = []
    for t in transects : 
      tmp_time = copy.deepcopy(transects[t][tag_idx]['tcorr']['sat_dates'])

      tmp_lon = transects[t]['transect'][0][0] + 5*(transects[t]['transect'][1][0] - transects[t]['transect'][0][0])
      tmp_lat = transects[t]['transect'][0][1] + 5*(transects[t]['transect'][1][1] - transects[t]['transect'][0][1])
      n_t.append(len(tmp_time))

      for j in range(n_t[-1]):
        time.append(tmp_time[j])
        lon.append(tmp_lon)
        lat.append(tmp_lat)
          
    VAR = {'time':(['t'],time),'longitude':(['t'],lon),'latitude':(['t'],lat)}
    COORD = {'t':(["t"], np.arange(len(time)))}
    ds_FES = xr.Dataset(data_vars = VAR,coords=COORD)
    ds_FES.to_netcdf(carrier_file)
      
    model = 'FES2022'
    ext = '-FES2022.nc'
    d = "/work/EOLAB/tools/FES/FES2022/FES2022"
    waves_files = os.path.join(d, "WAVE{}".format(ext))
    
    list_files = [x for x in os.listdir(d) if x.endswith(ext)]
    f = [s.replace(ext, '') for s in list_files]
    f = " ".join(f)
    
    prediction_file = carrier_file[:-3] + "-predictions2.nc"
    prediction_file_path = os.path.join(os.getcwd(),prediction_file)
    mesh = os.path.join(d,"mesh"+ ext)
    command = "predictor -p {} -a {} -g {} -o {} -w {}".format("".join(carrier_file), waves_files, mesh, prediction_file, f)
    #LOGGER.info("Launch predictor: {:s}".format(command))
    #print(command)
    os.system(command)
    
    ds_prediction = xr.load_dataset(prediction_file, engine='netcdf4')
    ds_prediction = xr.merge([ds_FES, ds_prediction.prediction])
    
    df_prediction = ds_prediction.to_dataframe()
    
    pred = df_prediction.prediction.values
    c=0
    i=0
    for t in transects:
        transects[t][tag_idx]['tcorr']['waterlevel']=pred[c:c+n_t[i]]
        c += n_t[i]
        i += 1


    
    c=-1
    for i in transects:
        c+=1
        print(i + '   ' + str(int(c/len(transects))))

        """Sat_slope_calculation"""
        X = copy.deepcopy(transects[i][tag_idx]['tcorr']['SDW_'+tag_idx])
        t = copy.deepcopy(transects[i][tag_idx]['tcorr']['sat_dates'])
        wl = copy.deepcopy(transects[i][tag_idx]['tcorr']['waterlevel'])
        
        if not True in np.isnan(wl):
            slopes = Tools.rangeSlopes(0.0, 0.2)
            Xall = Tools.wlCorrect(X,wl,slopes)
            freqMax = Tools.wlPeak(t,wl)
            slope ,uncSlope = Tools.integratePowerSpectrum(t,Xall,freqMax)
            transects[i]['post'] = dict()
            transects[i]['post']['slope'+tag_idx]=slope
            transects[i]['post']['uncertaintySlope'+tag_idx]=uncSlope
            slope = np.array(len(t)*[np.nanmean(slope)])

            transects[i][tag_idx]['tcorr']['slope'] = np.array(slope)
            transects[i][tag_idx]['tcorr'] = Tools.removeNaN(transects[i][tag_idx]['tcorr'], inputs, varname='slope')

            slope = copy.deepcopy(transects[i][tag_idx]['tcorr']['slope'])
            wl = copy.deepcopy(transects[i][tag_idx]['tcorr']['waterlevel'])
            transects[i][tag_idx]['tcorr']['SDW_'+tag_idx] += wl/slope
        
        else:
            transects[i][tag_idx].pop('tcorr')
            print('No tide data for transects '+str(i))

pickle.dump(transects,open('transects_post.p','wb'))



