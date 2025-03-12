import os
import pickle
import shutil
#
import osgeo
from shapely.geometry import Polygon
import geopandas as gpd
import yaml
from yaml.loader import SafeLoader
import sys
sys.path.insert(0, r".\Codes")
from Codes.functions import Tools

force=True

path_codes = './Codes'
path_work = os.getcwd()
path_inputs  = './Inputs'
file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join(path_inputs,file_inputs),'rb'),Loader = SafeLoader)

if not 'Projects' in os.listdir():
	os.mkdir('Projects')

path_project = os.path.join('./Projects',inputs['Project'])

try :
    os.mkdir(path_project)
    print('Project '+inputs['Project']+' has been created in Projects')
except:
    if force:
        pass
    else:
        print('The project : '+inputs['Project']+' already exists, if you wanna continue put force=True in this script.')
        raise
    
shutil.copy(os.path.join(path_inputs,'config.yaml'),os.path.join(path_project,'config.yaml'))
shutil.copy(os.path.join(path_codes,'template_postprocess.py'),os.path.join(path_project,'postprocess.py'))
shutil.copy('job_check.py',os.path.join(path_project,'job_check.py'))
shutil.copy('job_post.slurm',os.path.join(path_project,'job_post.slurm'))

if inputs['NatureROI']=='Transects':
    if inputs['CreateTransects']:
        Tools.createTransects(inputs,os.getcwd())
        transects = pickle.load(open(os.path.join('./DataExample',inputs['NewPathTransect']),'rb'))
    else:
        transects = pickle.load(open(os.path.join('./DataExample',inputs['PathTransect']),'rb'))
    
    ROIs,TRs = Tools.polyFromTransects(transects,d_ref = inputs['SizeROI'])

elif inputs['NatureROI']=='Polygon':
    ROIs = [[inputs['ROI_array']]]
    TRs=[[0]]
    transects=TRs

if inputs['SaveShpROI']:
    polygons = [Polygon(coords[0]) for coords in ROIs]
#    gdf = gpd.GeoDataFrame(geometry=polygons,crs = 4326)
#    gdf.to_file(os.path.join(path_project,'ROIs.shp'))
    
c=0
for i in ROIs:
    TR_out = dict()
    for j in TRs[c]:
        TR_out[j] = transects[j]
    tmp = 'poly_'+str(c)
    path_tmp = os.path.join(path_project,tmp)
    if not(tmp in os.listdir(path_project)):
        os.mkdir(path_tmp)
        pickle.dump(i,open(os.path.join(path_tmp,'poly.json'),'wb'))
        pickle.dump(TR_out,open(os.path.join(path_tmp,'transects.p'),'wb'))
    
        shutil.copy(os.path.join(path_codes,'template_download.py'),os.path.join(path_tmp,'download.py'))
        shutil.copy(os.path.join(path_codes,'template_main_process.py'),os.path.join(path_tmp,'main_process.py'))
        shutil.copy('job.slurm',os.path.join(path_tmp,'job.slurm'))
    # L=["python download.py \n", "python main_process.py \n"]
    # job = open(os.path.join(path_tmp,"job.slurm"),"a")
    # job.writelines(L)
    # job.close()
    if inputs['typeUse']=='HPC':
        os.chdir(path_tmp)
        os.system('sbatch job.slurm')
        os.chdir(path_work)
    c+=1
